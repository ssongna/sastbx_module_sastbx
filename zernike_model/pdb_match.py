from libtbx import easy_pickle
from stdlib import math as smath
from scitbx.array_family import flex
import os, sys
from sastbx.data_reduction import saxs_read_write
from sastbx.zernike_model import pdb2zernike, zernike_model, hcluster
from sastbx.basic_analysis import guinier_analyses
from scitbx import math, matrix
from iotbx import pdb

import iotbx.phil
from sastbx.pr.pregxs_tools import get_q_array_uniform_body
from scitbx.math import zernike_align_fft as fft_align
from iotbx.xplor import map as xplor_map
from cctbx import uctbx
from scitbx.golden_section_search import gss
import time

from sastbx.interface import get_input

master_params = iotbx.phil.parse("""\
query{
  target = None
  .type=path
  .help="the experimental intensity profile"

  nmax = 10
  .type=int
  .help="maximum order of zernike polynomial: FIXED for the existing database"

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

}

""")

banner = "-------------------Searching the protein DATABASE for similar shapes-------------------"


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


def read_pickle(path, dbprefix):
  nn_file = dbprefix+'.nn'
  code_file = dbprefix+'.codes'
  rmax_file = dbprefix+'.rmax'
  nn = easy_pickle.load( path+nn_file )
  codes = easy_pickle.load( path+code_file )
  rmaxs = easy_pickle.load( path+rmax_file )
  return nn, codes, rmaxs

def read_nlm(path, dbprefix):
  nlm_file = path+dbprefix + '.nlm'
  if( os.path.exists(nlm_file) ):
    nlm = easy_pickle.load( nlm_file )
    return nlm
  else:
    return None

class intoshape(object):
  def __init__(self, data, nmax=20, rmax=None, scan=True, fraction=0.9, smear=True, prefix=None, weight='i', delta_q=None):
    self.data=data
    self.smear=smear
    self.delta_q=delta_q
    self.setup_weighting_scheme(weight)
    if(self.smear):
      self.set_up_smear()
    self.int_scale = flex.max( self.data.i )
    self.nmax = nmax
    if( rmax is None):
      rmax = int( get_rg(data) * 3.0 / 2.0)
    self.rmax = rmax
    self.scan = scan
    self.prefix=prefix
    self.fraction = fraction
    self.nlm_array = math.nlm_array(nmax)
    self.nlm_total = self.nlm_array.coefs().size()
    self.nn_array = math.nl_array(nmax)
    self.nn = self.nn_array.nl()
    self.nn_total = self.nn_array.coefs().size()

  def setup_weighting_scheme( self, weight ):
    self.expt_var = flex.pow2( self.data.s )
    if(weight == 'i'):
      self.weights = self.data.i
    elif(weight == 's'):
      self.weights=self.data.s

  def set_up_smear(self):
    if(self.delta_q is None):
      self.delta_q = (self.data.q[1]-self.data.q[0])/10.0

  def lookup( self, coefs, codes, ntop):
    self.best_rmax = flex.double()
    self.best_codes = []
    self.coefs = []
    for c in coefs:
      self.coefs.append(c[0:self.nn_total])
    self.codes = codes
    self.ntop = ntop
    self.mean_ws=flex.sqrt( 1.0/flex.double( range(1,ntop+1) ) )

    if(self.scan):
      self.rmax_list=flex.double()
      self.top_hits = []
      self.scores = flex.double()

      self.rmax_max = self.rmax*2.0
      self.rmax_min = max( self.rmax/2.0, 1)
      for coef in self.coefs:
        self.this_coef = coef
        gss( self.search, self.rmax_min, self.rmax_max, eps=0.5, N=30 )
        self.rmax_list.append( self.this_rmax )
        self.scores.append( self.this_score )
      top_indices = self.tops( self.scores, self.ntop )
      for ii in top_indices:
        self.best_rmax.append( self.rmax_list[ii] )
        self.best_codes.append( codes[ii] )
      self.best_indices = top_indices

      self.plot_intensity(self.rmax_list, top_indices, self.scores, self.coefs, codes, qmax=0.5)


  def search(self, rmax ):
    z_model = zernike_model(self.nlm_array, self.data.q, rmax/self.fraction, self.nmax )
    self.nn_array.load_coefs( self.nn, self.this_coef )
    i_cal = z_model.calc_intensity( self.nn_array )
    scale = i_cal[0]
    i_cal = i_cal/(scale/self.int_scale)
    if(self.smear):
      new_i = smear_data( i_cal, self.data.q, self.delta_q )
      self.this_score = flex.sum_sq( ( new_i - self.data.i[1:-2] )/self.weights[1:-2])
    else:
      self.this_score = flex.sum_sq( ( i_cal - self.data.i )/self.weights )
    self.this_rmax = rmax
    return self.this_score



  def tops( self, x, n ):
    order = flex.sort_permutation( x )
    return order[0:n]

  def plot_intensity( self, rmax_list, indices, scores, coefs, codes, qmax=None ):
    prefix = self.prefix
    if qmax is None:
      q_array = self.data.q
    else:
      q_array = flex.double( range(int((qmax-self.data.q[0])*100)) )/100.0 + self.data.q[0]

    for nn in range( indices.size() ):
      ii = indices[ nn ]
      z_model = zernike_model(self.nlm_array, q_array, rmax_list[ii]/self.fraction, self.nmax )

      oo = indices[nn]
      self.nn_array.load_coefs(self.nn, self.coefs[oo] )
      i_cal = z_model.calc_intensity( self.nn_array )
      if(self.smear):
        new_i = smear_data( i_cal, q_array, self.delta_q )
        for ii in range(1, new_i.size()+1):
          i_cal[ii] = new_i[ii-1]
      i_cal = i_cal / i_cal[0] * self.int_scale
      filename = prefix+"_"+str(nn)+".iq"
      file = open(filename, 'w')
      print>>file, "# RMAX: "+str(self.best_rmax[nn])+"  TOP: "+str(nn)+"  SCORE: "+str( self.scores[oo] )+"  PDB CODE: "+str(codes[oo])
      for qq, ii in zip( q_array, i_cal ):
        print>>file, qq, ii
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
      nlm_total = fix_nlm_array.nlm().size()
      top_coefs = []
      mean=flex.double()
      sig=flex.double()

      for ii in self.best_indices:
        fix = nlm_coefs[ii][0:nlm_total]
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

    else: # There is no nlm coefs loaded
      comment = "# Coefficient distance, similar to the eq. (12) in L. Mak et al, JMGM.26 (2008) P.1035"
      all_nn_coefs = []
      for ii in range(self.ntop):
        nn_i = self.coefs[ self.best_models[ii] ].deep_copy()
        nn_i = nn_i/nn_i[0]
        all_nn_coefs.append( nn_i )
      for ii in range(self.ntop):
        for jj in range(ii+1):
          cc = (all_nn_coefs[ii]-all_nn_coefs[jj]).norm()
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
    clusters.print_dot(tree_dot_file)
    clusters.print_neato()

def get_rg(data, out=None):
    msga = guinier_analyses.multi_step_rg_engine( data, out)
    return msga.median_rg

def reduce_raw_data(raw_data, qmax, bandwidth, level=0.05, q_background=None):
    print "delta_q is ", bandwidth
    if qmax > raw_data.q[-1]:
      qmax = raw_data.q[-1]
    ### Get rid of noisy signall at very low q range ###
    qmin_indx = flex.max_index( raw_data.i )
    qmin = raw_data.q[qmin_indx]
    new_data = get_q_array_uniform_body(raw_data, q_min=qmin, q_max=qmax, level=level )
    qmax = new_data.q[-1]
    print "LEVEL=%f"%level, "and Q_MAX=%f"%qmax
    raw_q = raw_data.q[qmin_indx:]
    raw_i = raw_data.i[qmin_indx:]
    raw_s = raw_data.s[qmin_indx:]
    ### Take care of the background (set zero at very high q) ###
    if( q_background is not None):
      cutoff = flex.bool( raw_q > q_background )
      q_bk_indx = flex.last_index( cutoff, False )
      if( q_bk_indx < raw_q.size() ):
        bkgrd = flex.mean( raw_i[q_bk_indx:] )
        print "Background correction: I=I-background, where background=", bkgrd
        raw_i = flex.abs( raw_i -bkgrd )

    q = flex.double( range(int( (qmax-qmin)/bandwidth )+1) )*bandwidth + qmin
    raw_data.i = flex.linear_interpolation( raw_q, raw_i, q )
    raw_data.s = flex.linear_interpolation( raw_q, raw_s, q )
    raw_data.q = q

    return raw_data

def process( pdb_files, nmax, rmax ):
  pdb_models = []
  shift = (rmax, rmax, rmax)
  for file in pdb_files:
    mom_obj, vox_obj, ipdb = pdb2zernike.zernike_moments( file, nmax=nmax, coef_out=False, calc_intensity=False)
    if(mom_obj is not None):
      if(len(pdb_models)==0):
        ref_nlm_array = mom_obj.moments()
        pdb_models.append( pdb_model( mom_obj.moments().coefs(), file, vox_obj.rmax(), ipdb ) )
        ea=(0,0,0)
        write_pdb( file, vox_obj, ea, shift, ipdb )
      else:
        mov_nlm_array = mom_obj.moments()
        align_obj = fft_align.align(ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True)
        pdb_models.append( pdb_model(align_obj.moving_nlm.coefs(), file, vox_obj.rmax(), ipdb) )
        ea=align_obj.best_ea
        write_pdb( file, vox_obj, ea, shift, ipdb )
  return pdb_models

class pdb_model(object):
  def __init__(self, nlm_coef, filename, rmax, pdb_inp):
    self.nlm_coef = nlm_coef
    self.filename = filename
    self.rmax=rmax
    self.pdb_inp = pdb_inp

def write_pdb(filename, voxel_object, euler_angle, shift, pdbinput):
  base=filename.split('.')[0]
  base = base.split('/')
  n=len(base)
  base = base[n-1]
  print base
  # we need a more clever way tp write out pdb files
  out_pdb_name = base+"_sa.pdb"
  aligned_xyz = voxel_object.rotate((-euler_angle[0],euler_angle[1], -euler_angle[2]),False)+shift
  for a,xyz in zip( pdbinput.hierarchy.atoms(), aligned_xyz):
    a.set_xyz( new_xyz=xyz)
  pdbinput.hierarchy.write_pdb_file( file_name=out_pdb_name, open_append=False)


def build_map( nmax, shapes, coefs, codes, pdb_models):
  np_on_grid=30
  zga = math.zernike_grid(np_on_grid, nmax, False)

  ref_nlm_array = math.nlm_array(nmax)
  mov_nlm_array = math.nlm_array(nmax)

  nlm = ref_nlm_array.nlm()
  nlm_total = ref_nlm_array.coefs().size()

  top_cc = flex.double()
  ntop = shapes.ntop
  rank = 0
  ave_c = flex.complex_double(nlm_total, 0)
  if( pdb_models is not None):
    ref_nlm_array.load_coefs( nlm, pdb_models[0].nlm_coef )
    fraction = 0
  else:
    c = coefs[ shapes.best_indices[0] ][0:nlm_total]
    ave_c = c.deep_copy()
    ref_nlm_array.load_coefs( nlm, c )
    rank = 1
    filename = "m"+str(rank)+"_"+shapes.best_codes[ 0 ]+".xplor"
    fraction=write_xplor( zga, nlm, c,np_on_grid, shapes.best_rmax[0], filename )
    print rank, shapes.best_codes[0], fraction

  ref_mean, ref_s = get_mean_sigma( ref_nlm_array )

  mean_frac = fraction
  mean_sqr_frac = fraction*fraction
  for ii, code in zip(shapes.best_indices[rank:], shapes.best_codes[rank:]):
    coef = coefs[ii][0:nlm_total]
    mov_nlm_array.load_coefs( nlm, coef )
    mov_mean, mov_s = get_mean_sigma( mov_nlm_array )
    align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    new_c = align_obj.moving_nlm.coefs()
    cc =  align_obj.get_cc()
    top_cc.append( cc )
    ave_c = ave_c + new_c
    filename = "m"+str(rank+1)+"_"+code+".xplor"
    fraction = write_xplor( zga, nlm, new_c,np_on_grid, shapes.best_rmax[rank], filename )
    rank = rank + 1
    print rank, code, fraction
    mean_frac = mean_frac + fraction
    mean_sqr_frac = mean_sqr_frac + fraction*fraction

  #sphere_volume = 4.0/3.0*smath.pi*rmax*rmax*rmax
  rmax=shapes.best_rmax[0]
  sphere_volume = rmax*rmax*rmax*8.0
  mean_frac = mean_frac / ntop
  sigma_frac = smath.sqrt( mean_sqr_frac/ntop - mean_frac*mean_frac )
  print "Volume is ", mean_frac * sphere_volume, "+/-", sigma_frac*sphere_volume
  #### Write average map ####
  filename = "ave_map.xplor"
  write_xplor( zga, nlm, ave_c/ntop, np_on_grid, rmax, filename )
  return top_cc

def get_mean_sigma( nlm_array ):
  coef = nlm_array.coefs()
  mean = abs( coef[0] )
  var = flex.sum( flex.norm(coef) )
  sigma = smath.sqrt( var-mean*mean )
  return mean, sigma

def write_xplor( zga, nlm, coef, np_on_grid, rmax, filename):
  zga.load_coefs( nlm, coef )
  map = flex.abs( zga.f() )
  xplor_map_type( map, np_on_grid, rmax, filename)
  return volume_weight(map)

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
  print >> out, "   sastbx.shapeup <target=target.iq> [nmax=nmax scan=True*/False buildmap=True*/False pdb=pdbfile path=database_path]\n"
  print >> out, "   The intensity profile is the only required input file  (in theory)\n"
  print >> out, "   Optional control parameters:"
  print >> out, "     nmax     : maximum order of the zernike polynomial expansion (<=20 for precomputed database; 10 is the default)"
  print >> out, "     qmax     : maximum q value, beyond which the intensity profile will not be considered (default 0.20)"
  print >> out, "     path     : path to the database (this MUST be correct to execute the searching)"
  print >> out, "     buildmap : build electron density map in xplor format, all the map will be aligned"
  print >> out, "     pdb      : any pdb model to be compared, and the maps will be aligned to the first pdb file"
  print >> out, "     prefix   : the output prefix\n\n"

def set_default_db_path():
  import libtbx
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/database/'
  print "\nATTENT: database path was set to default:"
  print ">>>>  dbpath = ", path, "  <<<<"
  return path

def volume_weight(map, threshold=0.3):
  threshold = threshold * flex.max( map )
  molecule=flex.bool( map.as_1d() >= threshold )
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
  q_step = 1.0/100.0

  data = saxs_read_write.read_standard_ascii_qis( target_file )
  if(rmax is None):
    rmax = get_rg( data ) * 3.0 /2.0

  qmax = params.query.qmax
  q_background = params.query.q_background
  #qmax = 0.44*smath.exp( -0.00023*rmax*rmax )
  ######### Interpolation ##########
  bandwidth = min( q_step, data.q[2]/2.0) # smath.pi/2.0/rmax )
  data = reduce_raw_data( data, qmax, bandwidth, q_background=q_background,level=params.query.q_level )
  #saxs_read_write.write_standard_ascii_qis(data, 'reduced'+target_file )
  ###### END of Interpolation ##########

  nn_coefs, codes, rmaxs = read_pickle(dbpath, dbprefix)

  shapes = intoshape( data, nmax=nmax, rmax=rmax, scan=scan, fraction=fraction, smear=smear, prefix=params.query.prefix, weight=weight, delta_q=delta_q)
  shapes.lookup(nn_coefs, codes, ntop)

  pdb_models = None
  if( len(pdb_files) > 0 ):
    pdb_models = process(pdb_files,nmax, shapes.best_rmax[0]/fraction)

  nlm_coefs = None
  if(params.query.buildmap):
    nlm_coefs = read_nlm(dbpath, dbprefix)
    top_cc = build_map( nmax, shapes, nlm_coefs, codes, pdb_models )
                   # need to use rmax/fraction to get right size of box
    if(len(pdb_files) > 0 ):
      out=open(params.query.prefix+"_cc2pdb.dat", 'w')
      for cc in top_cc:
        print>>out, cc
      print>>out, "mean: %8.5f"%flex.mean(top_cc)
      print "mean cc: %8.5f"%flex.mean(top_cc)
      print "first cc: %8.5f"%top_cc[0]
      print "best cc: %8.5f"%flex.max(top_cc)
      print "worst cc: %8.5f"%flex.min(top_cc)
      out.close()
      print "Rmax: estimated vs PDB", shapes.best_rmax[0], pdb_models[0].rmax

  shapes.pair_align( nlm_coefs, params.query.calc_cc )

  t2 = time.time()
  print "total time used: ", t2-t1


if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
