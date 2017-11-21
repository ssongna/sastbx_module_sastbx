from libtbx import easy_pickle
from libtbx.test_utils import run_command
from stdlib import math as smath
from scitbx.array_family import flex
import os, sys
from sastbx.zernike_model import hcluster,model_interface
from scitbx import math, matrix
from iotbx import pdb, ccp4_map

import iotbx.phil
from libtbx.utils import Sorry
from scitbx.math import zernike_align_fft as fft_align
from iotbx.xplor import map as xplor_map
from cctbx import uctbx, sgtbx
from scitbx.golden_section_search import gss
import time

from sastbx.interface import get_input
from sastbx.zernike_model import error_model, model_consistency
from sastbx.zernike_model import build_pymol_script, generate_html
from sastbx.intensity.sas_I import write_json
from sastbx.data_reduction import curves

#### fXS part ####
from sastbx.fXS import fxs_tools

master_params = iotbx.phil.parse("""\
query{
  target = None
  .type=path
  .help="the experimental intensity profile"

  nmax = 10
  .type=int
  .help="maximum order of zernike polynomial: FIXED for the existing database"

  lmax = None
  .type=int
  .help="maximum order of Legendre polynomial: FIXED for the existing database"

  pdb_files=None
  .multiple=True
  .type=path
  .help="If provided, align this structure to the models"

  qmax = 0.20
  .type=float
  .help="maximum q value where data beyond are disgarded"

  q_level = 0.01
  .type=float
  .help="ratio between I_stop and I_max"

  q_background=None
  .type = float
  .help = "the intensity beyond q-background is treated as background"

  rmax = None
  .type=float
  .help="estimated rmax of the molecule"

  scan = True
  .type = bool
  .help = "scan for different rmax?"

  prefix="query"
  .type=path
  .help="the output prefix"

  dbpath=None
  .type=path
  .help="the directory of database file, i.e., the pickle files"

  db_choice = *pisa piqsi allpdb user
  .type = choice
  .help = "Data base name"

  db_user_prefix="mydb"
  .type=path
  .help="the prefix of database filename"

  buildmap=True
  .type=bool
  .help="align the top models and generate xplor files"

  calc_cc=True
  .type=bool
  .help="calculate Correlation Coefficient or just Coefficient distance"

  smear=True
  .type=bool
  .help="smear the calculated data to remove the spikes (fits better to expt data)"

  weight=*i s
  .type=choice
  .help="the weights to be used in chi-score calculation"

  delta_q = None
  .type=float
  .help="linear smearing distance, default is set to q_step*0.1"

  ntop = 10
  .type=int
  .help="number of top hits returned per search"

  fraction = 0.9
  .type=float
  .help="fraction in zernike moments calculation on 1-D axis: This is FIXED, unless the database is changed"

  scale_power = 4.0
  .type = float
  .help = "Parameter controlling the scale factor calculation. Default should be good."


}

""")

banner = "-------------------Searching the protein DATABASE for similar shapes-------------------"


def build_3d_grid( np_on_grid, rmax_over_fraction ):
  grid = flex.vec3_double()
  one_d = range( -np_on_grid, np_on_grid+1 )
  one_d_x = flex.double(one_d)/np_on_grid*rmax_over_fraction
  for x in one_d_x:
    for y in one_d_x:
      for z in one_d_x:
        grid.append([x,y,z])
  return grid

def get_moments_for_scaled_model( map, np_on_grid, grid, nmax, rmax, external_rmax ):
  ### default params for moments calculation ###
  splat_range=1
  uniform=True
  fix_dx=False
  default_dx=0.7
  fraction=0.9
  ### end of default params for moments calculation ###

  threshold = flex.max( map )/3.0
  select = flex.bool( map.as_1d()>=threshold )
  xyz = grid.select( select )
  xyz_norms=xyz.norms()
  select = flex.bool ( xyz_norms<=rmax )
  xyz = xyz.select( select )
  density = flex.double(xyz.size(), 1.0)
  vox_obj = math.sphere_voxel( np_on_grid, splat_range, uniform, fix_dx, external_rmax, default_dx, fraction, xyz, density )
  grid_obj = math.sphere_grid( np_on_grid, nmax )
  grid_obj.clean_space(vox_obj, False)
  grid_obj.construct_space_sum()
  moment_obj = math.zernike_moments( grid_obj, nmax)
  nlm_array=moment_obj.moments()
  return nlm_array.coefs()


def smear_data( I_s, q_s, delta_q ):
  tmp_q1=q_s - delta_q
  tmp_q2=q_s + delta_q
  tmp_i1 = flex.linear_interpolation( q_s, I_s, tmp_q1[1:-2] )
  tmp_i2 = flex.linear_interpolation( q_s, I_s, tmp_q2[1:-2] )
  new_i = (tmp_i1 + tmp_i2)/2.0
  return new_i

def gauss_smear_data( I_s, w_s, np_total, span, f2c_ratio ):
  new_I = flex.double(np_total, 0)
  for ii in range( np_total ):
    jj = (ii+span)*f2c_ratio
    new_I[ii] = flex.sum( I_s[jj-span:jj+span+1]*w_s )
  return new_I


def read_codes(path, dbprefix):
  code_file = dbprefix+'.codes'
  if not os.path.isfile(path+code_file):
    raise Sorry("Datbase file %s is not present. Please check input parameters"%code_file)

  codes = easy_pickle.load( path+code_file )
  return codes

def read_nlm(path, dbprefix):
  nlm_file = path+dbprefix + '.nlm'
  if( os.path.exists(nlm_file) ):
    nlm = easy_pickle.load( nlm_file )
    return nlm
  else:
    return None

class intoshape(object):
  def __init__(self, data, rg,io, nmax=20, lmax=10, rmax=None, scan=True, fraction=0.9, smear=True, prefix=None, weight='i', delta_q=None, scale_power=2.0):
    self.data=data
    blq0 = self.data.blq[0]
    self.data.blq = self.data.blq/blq0
    self.setup_weighting_scheme( weight )
    self.smear=smear
    self.rg = rg
    self.io = io
    self.int_scale =  self.io
    self.nmax = nmax
    self.lmax = lmax
    if( rmax is None):
      rmax = int( self.rg * 3.0 / 2.0)
    self.rmax = rmax
    self.scan = scan
    self.prefix=prefix
    self.fraction = fraction
    self.nlm_array = math.nlm_array(nmax)
    self.nlm = self.nlm_array.nlm()
    self.nlm_total = self.nlm_array.coefs().size()
    self.scale_power  = scale_power


  def setup_weighting_scheme( self, weight ):
    if(weight == 'i'):
      self.weights = self.data.blq
    elif(weight == 's'):
      self.weights=flex.sqrt( self.data.blq )


  def lookup( self, coefs, codes, ntop):
    self.out = open(self.prefix+'.sum', 'w')
    self.coefs = []
    for c in coefs:
      self.coefs.append(c[0:self.nlm_total])
    self.codes = codes
    self.ntop = ntop
    self.mean_ws=flex.sqrt( 1.0/flex.double( range(1,ntop+1) ) )

    if(self.scan):

      self.rmax_list=flex.double()
      self.top_hits = []
      self.scores = []
      self.scales = []
      self.ave_scores = flex.double()
      #self.rmax_max = 3.14/self.data.q[0]
      self.rmax_max = self.rmax*2.5
      self.rmax_min = max( self.rmax/2.0, 1)
      print "   Search range of rmax  :   %5.2f A  ----  %5.2f A"%(self.rmax_max, self.rmax_min)
      gss( self.score_at_rmax, self.rmax_min, self.rmax_max, eps=0.5, N=30,monitor_progress=True )
      rmax_indx = flex.min_index( self.ave_scores )
      self.best_rmax = self.rmax_list[ rmax_indx ]
      self.best_models = self.top_hits[rmax_indx]
      print "   Best rmax found       :   %5.2f A"%self.best_rmax

      print >>self.out, "Best Result from Golden Section Search:", self.best_rmax
      self.show_result(self.best_rmax, self.best_models, self.scores[rmax_indx], codes, self.out)
      self.plot_intensity(self.best_rmax, self.best_models, self.scores[rmax_indx], self.coefs, codes, qmax=None,scales=self.scales[rmax_indx])
      self.print_rmax_profile( self.rmax_list, self.ave_scores )
      self.summary( self.top_hits, self.scores, comment="----Statistics from Golden Section Scan----" )

      self.local_scan( int( self.best_rmax + 0.5) , local=5 )

    else:
      print
      print "   Not performing search in rmax. Using fixed value of rmax=%5.3f"%self.rmax
      print
      self.score_at_rmax( self.rmax )
      self.plot_intensity(self.best_rmax, self.best_models, self.scores, self.coefs, codes, qmax=None)

    self.out.close()

  def local_scan( self, rmax, local=5 ):
    print >>self.out, "Results from Local Scan around the optimal rmax:"
    for r in range(rmax-local, rmax+local+1):
      self.score_at_rmax( r )
    self.summary( self.top_hits[-local*2:], self.scores[-local*2:], comment="----Statistics from Local Scan----" )

  def summary(self, top_models, top_scores,  comment="---"):
    model_dict={}
    score_dict={}
    for models, scores in zip(top_models, top_scores):
      for m, i in zip(models, range(self.ntop) ):
        if( model_dict.__contains__(m)):
          model_dict[m] = model_dict[m] + 1
          score_dict[m] = score_dict[m] + scores[i]
        else:
          model_dict[m] = 1
          score_dict[m] = scores[i]
    model_dict.items().sort()
    out=open(self.prefix+'.sta', 'a')
    keys = model_dict.keys()
    counts = flex.double( model_dict.values() )
    scores = flex.double( score_dict.values() ) / counts

    print>>out, comment
    print>>out, "Sorted by Appearances:"
    order = flex.sort_permutation( counts )
    for o in order:
      print>>out, self.codes[keys[o]], counts[o], scores[o]

    print>>out, "Sorted by Average Score:"
    order = flex.sort_permutation( scores )
    for o in order:
      print>>out, self.codes[keys[o]], counts[o], scores[o]
    out.close()


  def find_average_score(self, top_scores):
    diff = top_scores[1:] - top_scores[0:self.ntop-1]
    for ii in range(1, self.ntop ):
      if( diff[ii]<diff[ii-1] ):
        break
    nsample=ii+1
    ave_score = flex.mean( top_scores[0:nsample] )
    return ave_score

  def weight_mean( self, top_scores ):
    return flex.mean( top_scores*self.mean_ws )


  def score_at_rmax(self, rmax):
    top_hits, scores, scales = self.fixed_rmax( rmax )
    if( self.scan ):
      self.rmax_list.append( rmax )
      self.top_hits.append( top_hits )
      self.scores.append( scores )
      self.scales.append( scales )
      ave_score = flex.mean( flex.sqrt(scores) )
      self.ave_scores.append( ave_score )
      return ave_score
    else:
      self.scores = scores
      self.scales = scales
      self.best_models = top_hits
      self.best_rmax = rmax
      return


  def fixed_rmax(self, rmax):
    z_model = fxs_tools.znk_blq(self.nlm_array, self.data.q, rmax/self.fraction, self.nmax, self.lmax )
    top_hits, scores, scales = self.search( z_model, self.coefs, self.ntop, rmax )
    self.show_result( rmax, top_hits, scales, self.codes, self.out)
    return top_hits, scores, scales


  def get_rg(self, data, out=None):
    msga = guinier_analyses.multi_step_rg_engine( data, out)
    return msga.median_rg, msga.median_io

  def score_cc(self,obs_data, calc_data, weights):
    weights = weights*weights*weights
    obj = math.weighted_covariance(obs_data, calc_data, weights)
    return obj.correlation


  def approx_scale( self, obs, calc, w ):
    rats = flex.sum( calc * w )/flex.sum( obs *w  )
    return rats

  def search(self, z_model, coefs, ntop, rmax ):
    scores = flex.double()
    scales = flex.double()
    w = 1.0 ## temporary use
    for coef in coefs:
      self.nlm_array.load_coefs( self.nlm, coef )
      calc_blq = z_model.get_all_blq( self.nlm_array )
      scale = calc_blq[0]
      calc_blq = calc_blq/scale
      scale = self.approx_scale(self.data.blq, calc_blq,w)
      score = flex.sum_sq( ( calc_blq - self.data.blq*scale )/(scale*w) )

      scores.append( score )
      scales.append( scale )

    top_hits = self.tops( scores, ntop )
    top_scores = flex.double()
    top_scales = flex.double()
    for tt in top_hits:
      top_scores.append( scores[tt] )
      top_scales.append( scales[tt] )

    return top_hits, top_scores, top_scales


  def show_result( self, rmax, top_hits, scores, codes, out):
    print>>out, "RMAX:: ", rmax
    for ii in range( top_hits.size() ):
      indx = top_hits[ii]
      print>>out,  "The top %s match is "%ii, codes[indx], scores[ii]

  def print_rmax_profile( self, rmax_list, ave_scores ):
    out=open(self.prefix+".rmax_profile", 'w')
    for r, s in zip( rmax_list, ave_scores):
      print>>out, r, s

  def tops( self, x, n ):
    order = flex.sort_permutation( x )
    return order[0:n]

  def plot_intensity( self, rmax, top_hits, scores, coefs, codes, qmax=None, scales=None ):
    prefix = self.prefix
    if qmax is None:
      q_array = self.data.q
    else:
      q_array = flex.double( range(int((qmax-self.data.q[0])*100)) )/100.0 + self.data.q[0]
    z_model = fxs_tools.znk_blq(self.nlm_array, q_array, rmax/self.fraction, self.nmax, self.lmax )

    for nn in range( top_hits.size() ):
      oo = top_hits[nn]
      self.nlm_array.load_coefs(self.nlm, self.coefs[oo] )
      cal_blq = z_model.get_all_blq( self.nlm_array )
      cal_blq = cal_blq/cal_blq[0]
      filename = prefix+"_"+str(nn+1)+".blq"
      file = open(filename, 'w')
      print>>file, "# RMAX: "+str(rmax)+"  TOP: "+str(nn)+"  SCORE: "+str( scores[nn] )+"  PDB CODE: "+str(codes[oo])
      self.data.print_out( cal_blq, out=file )
      file.close()


  def pair_align(self, nlm_coefs, calc_cc=True):
    self.cc_array = []
    for ii in range( self.ntop ):
      self.cc_array.append( flex.double( self.ntop, 1 ) )

    if( nlm_coefs is not None  and calc_cc):
      comment = "# Correlation Coefficient <rho_1(r)*rho_2(r)>"
      fix_nlm_array = math.nlm_array(self.nmax)
      mov_nlm_array = math.nlm_array(self.nmax)
      nlm = fix_nlm_array.nlm()
      nlm_total = fix_nlm_array.coefs().size()
      top_coefs = []
      mean=flex.double()
      sig=flex.double()

      for ii in range(self.ntop):
        ff = self.best_models[ii]
        fix=nlm_coefs[ff][0:nlm_total]
        top_coefs.append( fix )
        fix_nlm_array.load_coefs( nlm, fix )
        m,s = get_mean_sigma( fix_nlm_array )
        mean.append( m )
        sig.append( s )

      for ii in range(self.ntop):
        fix=top_coefs[ii]
        fix_nlm_array.load_coefs( nlm, fix )
        for jj in range(ii):
          mov=top_coefs[jj]
          mov_nlm_array.load_coefs( nlm, mov )
          cc = fft_align.align( fix_nlm_array, mov_nlm_array, nmax=self.nmax, refine=True ).best_score
          cc = (cc-mean[ii]*mean[jj])/(sig[ii]*sig[jj])
          self.cc_array[ii][jj]=cc
          self.cc_array[jj][ii]=cc

    outfile = self.prefix+".cc"
    out=open(outfile, 'w')
    print>>out, comment
    for ii in range(1,self.ntop+1):
      print>>out,"%6d"%ii,
    print>>out, "   average"

    for ii in range(self.ntop):
      for jj in range(self.ntop):
        print>>out,"%6.3f"%self.cc_array[ii][jj],
      print>>out, "%6.3f"%flex.mean( self.cc_array[ii] )

    clusters = hcluster.hcluster(self.cc_array, 0.8)
    clusters.print_hclust()
    out.close()
    tree_dot_file = self.prefix+".tree"
    clusters.print_neato(tree_dot_file)
    # generate image using neato
#    run_command('/sw/bin/neato -Tpng -o cluster.png '+tree_dot_file )
    #clusters.print_neato()
    self.clusters = clusters
    return

def get_rg(data, out=None):
    msga = guinier_analyses.multi_step_rg_engine( data, out)
    return msga.median_rg, msga.median_io


def process( pdb_files, nmax, rmax=50.0, fraction=0.9 ):
  pdb_models = []
  rmax_over_fraction = rmax/fraction
  shift = (rmax_over_fraction,rmax_over_fraction, rmax_over_fraction)
  for file in pdb_files:
    model = model_interface.container( pdbfile=file, rmax=rmax, nmax=nmax )
    if(model is not None):
      if(len(pdb_models)==0):
        ref_nlm_array = model.nlm_array
        pdb_models.append( pdb_model( ref_nlm_array.coefs(), file, model.rmax, model) )
        ea=(0,0,0)
        outname = model.id+'_sa.pdb'
        model.write_pdb( rmax=rmax_over_fraction,rotation=ea, filename=outname )
      else:
        mov_nlm_array = model.nlm_array
        align_obj = fft_align.align(ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True)
        pdb_models.append( pdb_model(align_obj.moving_nlm.coefs(), file, model.rmax, model) )
        ea=align_obj.best_ea
        outname = model.id+'_sa.pdb'
        model.write_pdb( rmax=rmax_over_fraction,rotation=ea, filename=outname )
  return pdb_models

class pdb_model(object):
  def __init__(self, nlm_coef, filename, rmax, pdb_model):
    self.nlm_coef = nlm_coef
    self.filename = filename
    self.rmax=rmax
    self.pdb_model = pdb_model


def build_map( nmax, rmax, coefs, codes, model_indx, pdb_models, clusters=None, fract=0.9, type='.ccp4'):
  rmax_over_fraction = rmax/fract
  np_on_grid=30
  zga = math.zernike_grid(np_on_grid, nmax, False)

  ref_nlm_array = math.nlm_array(nmax)
  mov_nlm_array = math.nlm_array(nmax)

  nlm = ref_nlm_array.nlm()
  nlm_total = ref_nlm_array.coefs().size()

  top_cc = flex.double()
  top_ids = []
  ntop = model_indx.size()
  rank = 0
  ave_c = flex.complex_double(nlm_total, 0)

  aligned_coefs = []
  map_files = []
  levels = []
  scale_model=False
  if( (pdb_models is not None) and (pdb_models[0].rmax > rmax ) ):
    scale_model=True
    external_rmax=pdb_models[0].rmax
    grid=build_3d_grid(np_on_grid, rmax_over_fraction)

  print "Rank PDB_code cc (to the given model or the first model):"
  if( pdb_models is not None):
    ref_nlm_array.load_coefs( nlm, pdb_models[0].nlm_coef )
    fraction = 0
  else:
    c = coefs[model_indx[0] ][0:nlm_total]
    ave_c = c.deep_copy()
    ref_nlm_array.load_coefs( nlm, c )
    rank = 1
    filename = "m"+str(rank)+"_"+codes[ model_indx[0] ]
    map_files.append( filename+type )
    fraction,map=write_map( zga, nlm, c,np_on_grid, rmax_over_fraction, filename, type=type )
    level = flex.max(map)/3.0
    levels.append(level)
    top_cc.append( 1.0 )  # refering to itself
    if( scale_model ):
      c = get_moments_for_scaled_model( map, np_on_grid, grid, nmax, rmax, external_rmax )
    print rank, codes[ model_indx[0] ]
    aligned_coefs.append( c.deep_copy() )# save the aligned nlm coefs

  mean_frac = fraction
  mean_sqr_frac = fraction*fraction
  for ii in model_indx[rank:]:
    rank = rank + 1
    c = coefs[ ii ][0:nlm_total]
    mov_nlm_array.load_coefs( nlm, c )
    align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    new_c = align_obj.moving_nlm.coefs()
    filename = "m"+str(rank)+"_"+codes[ii]
    map_files.append( filename+type )
    fraction,map = write_map( zga, nlm, new_c,np_on_grid, rmax_over_fraction, filename, type=type )
    if( scale_model ):
      c = get_moments_for_scaled_model( map, np_on_grid, grid, nmax, rmax, external_rmax )
      mov_nlm_array.load_coefs( nlm, c )
      align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
      new_c = align_obj.moving_nlm.coefs()
      fraction,map = write_map( zga, nlm, new_c,np_on_grid, rmax_over_fraction, filename, type=type )
    level = flex.max(map)/3.0
    levels.append(level)
    cc =  align_obj.get_cc()
    print "%2d  %5s  %5.3f"%(rank, codes[ ii ], cc )
    top_cc.append( cc )
    top_ids.append( codes[ii] )
    ave_c = ave_c + new_c
    aligned_coefs.append( new_c.deep_copy() )  # save the aligned nlm coefs
    mean_frac = mean_frac + fraction
    mean_sqr_frac = mean_sqr_frac + fraction*fraction

  sphere_volume = rmax_over_fraction**3.0*8.0  # cube with d=2.0*r
  mean_frac = mean_frac / ntop
  sigma_frac = smath.sqrt( mean_sqr_frac/ntop - mean_frac*mean_frac )
  print "Volume is ", mean_frac * sphere_volume, "+/-", sigma_frac*sphere_volume, "(A^3)"
  #### Write average map ####
  ave_maps = []
  ave_levels=[]
  ave_cc=[]
  cluster_ids=[1]*len(model_indx)  # initialize cluster_ids
#  filename = "ave_map"
#  map_files.append( filename+type )
#  fraction, map=write_map( zga, nlm, ave_c/ntop, np_on_grid, rmax_over_fraction, filename, type=type )
#  levels.append( flex.max(map)/3.0 )
#  if( len(clusters.nodes) == 1): return top_cc, top_ids, map_files, levels, [1]*len(map_files)
  cluster_id = 1
  print "cc. between Cluster average and PDB model"
  for node in clusters.nodes:
    ave_c = ave_c*0
    coefs_list = []
    for ii in node.leaf_eles:
      ave_c = ave_c + aligned_coefs[ii]
      cluster_ids[ii]=cluster_id
      coefs_list.append( aligned_coefs[ii] )
    ave_c = ave_c/len(node.leaf_eles)
    level_n = model_consistency.good2n( nmax, coefs_list, ave_c )
    print "consistency level to order n: %d"%level_n
    mov_nlm_array.load_coefs( nlm, ave_c )
    align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    cc = align_obj.get_cc()
    ave_cc.append(cc)
    print "cluster  # ", cluster_id, "cc=", cc
    filename = "ave_"+str(cluster_id)
    ave_maps.append( filename+type )
    fraction, map = write_map( zga, nlm, ave_c, np_on_grid, rmax_over_fraction, filename, type=type )
    ave_levels.append( flex.max(map)/3.0 )
    print "Averaged Model #%d Volume is %f (A^3)"%(cluster_id, fraction * sphere_volume)
    cluster_id = cluster_id+1

  return top_cc, top_ids, map_files, levels, cluster_ids, ave_maps, ave_levels, ave_cc

def get_mean_sigma( nlm_array ):
  coef = nlm_array.coefs()
  mean = abs( coef[0] )
  var = flex.sum( flex.norm(coef) )
  sigma = smath.sqrt( var-mean*mean )
  return mean, sigma

def write_map( zga, nlm, coef, np_on_grid, rmax, filename, type='.ccp4'):
  zga.load_coefs( nlm, coef )
  map = flex.abs( zga.f() )
  if( type=='.ccp4'):
    ccp4_map_type( map, np_on_grid, rmax, filename+'.ccp4')
  else:
    xplor_map_type( map, np_on_grid, rmax, filename+'.xplor')
  return volume_weight(map), map

def ccp4_map_type(map, N, radius,file_name='map.ccp4'):
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  map.reshape( grid )
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90"),
      space_group=sgtbx.space_group_info("P1").group(),
      gridding_first=(0,0,0),
      gridding_last=(N*2, N*2, N*2),
      map_data=map,
      labels=flex.std_string(["generated from zernike moments"]))


def xplor_map_type(m,N,radius,file_name='map.xplor'):
  gridding = xplor_map.gridding( [N*2+1]*3, [0]*3, [2*N]*3)
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  m.reshape( grid )
  uc = uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90")
  xplor_map.writer( file_name, ['no title lines'],uc, gridding,m)


def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\nUsage: \n"
  print >> out, "   sastbx.shapeup <target=target.iq> [rmax=rmax nmax=nmax scan=True*/False buildmap=True*/False pdb=pdbfile path=database_path]\n"
  print >> out, "   The intensity profile is the only required input file  (in theory)\n"
  print >> out, "   Optional control parameters:"
  print >> out, "     rmax     : radius of the molecule (default: guessed from Rg)"
  print >> out, "     nmax     : maximum order of the zernike polynomial expansion (<=20 for precomputed database; 10 is the default)"
  print >> out, "     qmax     : maximum q value, beyond which the intensity profile will not be considered (default 0.20)"
  print >> out, "     path     : path to the database (this MUST be correct to execute the searching)"
  print >> out, "     buildmap : build electron density map in xplor format, all the map will be aligned"
  print >> out, "     pdb      : any pdb model to be compared, and the maps will be aligned to the first pdb file"
  print >> out, "     prefix   : the output prefix\n\n"

def set_default_db_path():
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/database/'
  print "\nATTENTION: Database path was set to : >>%s<<"%path
  return path

def volume_weight(map, threshold=0.3):
  #threshold = threshold * flex.max( map )
  map_1d = map.as_1d()
  mean_and_var = flex.mean_and_variance(map_1d )
  threshold = mean_and_var.mean() - mean_and_var.unweighted_standard_error_of_mean()
  molecule=flex.bool( map_1d >= threshold )
  molecule_tot = flex.sum( molecule.as_double() )
  fraction = molecule_tot/map.size()
  return fraction


def run(args):
  t1 = time.time()
  params = get_input(args, master_params, "query", banner, help)
  if( params is None):
    exit()
  target_file = params.query.target
  rmax = params.query.rmax
  nmax = params.query.nmax
  lmax = params.query.lmax
  if( lmax is None ): lmax=nmax
  smear = params.query.smear
  dbpath = params.query.dbpath
  pdb_files = params.query.pdb_files
  db_choice = params.query.db_choice
  weight = params.query.weight
  delta_q = params.query.delta_q
  if( db_choice == "user" ):
    dbprefix = params.query.db_user_prefix
  else:
    dbprefix = db_choice

  if (dbpath is None):
    dbpath = set_default_db_path()
  ntop = params.query.ntop
  scan = params.query.scan
  fraction = params.query.fraction
  scale_power = params.query.scale_power
  q_step = 1.0/100.0

  data = fxs_tools.read_blq( target_file, lmax=lmax )
  saxs_i = data.blq[0:-1:int(lmax/2)+1]
  saxs_i = flex.sqrt( flex.abs(saxs_i) )
  saxs_data = curves.simple_saxs_data(data.q, saxs_i, saxs_i)
  try:
    rg, io = get_rg(saxs_data)
  except:
    print "Guinier analysis failed, R_max is required"
    print "ATTENTION: dummy values for Rg and Io set"
    rg=50
    io=1

  qmax = params.query.qmax
  q_background = params.query.q_background

  print " ==== Reading in shape database ==== "
  begin_time = time.time()
  nlm_coefs = read_nlm(dbpath, dbprefix)
  codes = read_codes(dbpath, dbprefix)
  ready_time = time.time()
  delta_time = ready_time-begin_time
  print
  print "   Done reading database with %i entries in %5.4e seconds"%(len(codes),delta_time)
  print

  print " ==== Shape retrieval ==== "
  print "   Constructing shape retrieval object"
  shapes = intoshape( data, rg=rg, io=io, nmax=nmax, lmax=lmax, rmax=rmax, scan=scan, fraction=fraction, smear=smear, prefix=params.query.prefix, weight=weight, delta_q=delta_q,scale_power=scale_power)
  print "   Shape search  .... "
  shapes.lookup(nlm_coefs, codes, ntop)

  shapes.pair_align( nlm_coefs, params.query.calc_cc )

  pdb_models=None
  if( len(pdb_files) > 0 ):
    pdb_models = process(pdb_files,nmax, rmax=shapes.best_rmax, fraction=fraction)

  if(params.query.buildmap):
    top_cc, top_ids, map_files, levels, cluster_ids, ave_maps, ave_levels, ave_cc = build_map( nmax, shapes.best_rmax, nlm_coefs, codes, shapes.best_models, pdb_models, clusters=shapes.clusters,fract=fraction )
                   # need to use rmax/fraction to get right size of box
    build_pymol_script.write_pymol_scripts(map_files,levels)
    pdb_out_name=None
    if( pdb_models is not None ):
      pdb_out_name = pdb_files[0].split('.')[0]+'_sa.pdb'
    generate_html.generate_jmol_html(ave_maps, ave_cc, ave_levels, map_files, top_cc, levels, cluster_ids, 'models.html', pdb=pdb_out_name)
    if(len(pdb_files) > 0 ):
      out=open(params.query.prefix+"_cc2pdb.dat", 'w')
      print >> out, "Correlation coefficients of retrieved shapes vs input model"
      for cc,id in zip(top_cc,top_ids):
        print>>out,"Code: %5s    CC: %5.1f  "%(id,100*cc)

      print>>out, "mean: %8.5f"%flex.mean(top_cc)
      print "Compared to the PDB model (%s)"%pdb_models[0].filename
      print "mean cc: %8.5f"%flex.mean(top_cc)
      print "first cc: %8.5f"%top_cc[0]
      print "best cc: %8.5f"%flex.max(top_cc)
      print "worst cc: %8.5f"%flex.min(top_cc)
      out.close()
      print "Rmax: estimated vs PDB", shapes.best_rmax, pdb_models[0].rmax

  t2 = time.time()
  print "total time used: ", t2-t1, "(seconds)"


if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
