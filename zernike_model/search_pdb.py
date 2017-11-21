from libtbx import easy_pickle
from libtbx.test_utils import run_command
from stdlib import math as smath
from scitbx.array_family import flex
import os, sys
from sastbx.data_reduction import saxs_read_write
from sastbx.zernike_model import pdb2zernike, zernike_model, hcluster
from sastbx.basic_analysis import guinier_analyses
from scitbx import math, matrix
from iotbx import pdb, ccp4_map

import iotbx.phil
from libtbx.utils import Sorry
from sastbx.pr.pregxs_tools import get_q_array_uniform_body
from scitbx.math import zernike_align_fft as fft_align
from iotbx.xplor import map as xplor_map
from cctbx import uctbx, sgtbx
from scitbx.golden_section_search import gss
import time

from sastbx.interface import get_input
from sastbx.zernike_model import error_model, model_consistency
from sastbx.zernike_model import build_pymol_script, generate_html
from sastbx.intensity.sas_I import write_json



base_path = os.path.split(sys.path[0])[0]
global stdfile
global outfilelog
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

  #threshold = flex.max( map )/3.0
  threshold = map.standard_deviation_of_the_sample()
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


def read_pickle(path, dbprefix):
  nn_file = dbprefix+'.nn'
  code_file = dbprefix+'.codes'
  rmax_file = dbprefix+'.rmax'
  if not os.path.isfile(path+nn_file):
    raise Sorry("Datbase file %s is not present. Please check input parameters"%nn_file)
  if not os.path.isfile(path+code_file):
    raise Sorry("Datbase file %s is not present. Please check input parameters"%code_file)
  if not os.path.isfile(path+rmax_file):
    raise Sorry("Datbase file %s is not present. Please check input parameters"%rmax_file)

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
  def __init__(self, data, rg,io, nmax=20, rmax=None, scan=True, fraction=0.9, smear=True, prefix=None, weight='i', delta_q=None, scale_power=2.0):
    global stdfile
    global outfilelog
    self.stdfile = stdfile
    self.outfilelog = outfilelog
    self.data=data
    self.smear=smear
    self.delta_q=delta_q
    self.setup_weighting_scheme(weight)
    if(self.smear):
      self.set_up_smear()


    self.rg = rg
    self.io = io
    self.int_scale =  self.io
    self.nmax = nmax
    if( rmax is None):
      rmax = int( self.rg * 3.0 / 2.0)
    self.rmax = rmax
    self.scan = scan
    self.prefix=prefix
    self.fraction = fraction
    self.nlm_array = math.nlm_array(nmax)
    self.nlm_total = self.nlm_array.coefs().size()
    self.nn_array = math.nl_array(nmax)
    self.nn = self.nn_array.nl()
    self.nn_total = self.nn_array.coefs().size()
    self.scale_power  = scale_power


  def setup_weighting_scheme( self, weight ):
    self.expt_var = flex.pow2( self.data.s )
    if(weight == 'i'):
      self.weights = self.data.i
    elif(weight == 's'):
      self.weights=self.data.s

  def update_model_error(self, rmax, a=-19.93, b=18.95,c=0.52, d=-1.13):
    self.ratio = self.data.q*rmax
    self.ratio = flex.pow( self.ratio, c)*d
    self.ratio = a+b/(flex.exp( self.ratio ) + 1.0 )
    self.ratio = flex.exp( self.ratio )

  def set_up_smear(self):
    if(self.delta_q is None):
      self.delta_q = (self.data.q[1]-self.data.q[0])/10.0

  def variance_est(self,q,i,rmax,a=1.0,b=0.14,c=55.0,ref_rmax=50.0):
    frac_sig = a/( 1+ flex.exp( (b-q)*c*rmax/ref_rmax))
    sig = frac_sig*i
    return sig


  def lookup( self, coefs, codes, ntop):
    global stdfile
    global outfilelog
    with open(outfilelog,"a") as filelog:
      print >> filelog, self.prefix+'.sum'

    self.out = open(self.prefix+'.sum', 'w')
    self.coefs = []
    for c in coefs:
      self.coefs.append(c[0:self.nn_total])
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
      with open(stdfile,"a") as log:
        print >>log,  "   Search range of rmax  :   %5.2f A  ----  %5.2f A"%(self.rmax_max, self.rmax_min)
      print   "   Search range of rmax  :   %5.2f A  ----  %5.2f A"%(self.rmax_max, self.rmax_min)

      gss( self.score_at_rmax, self.rmax_min, self.rmax_max, eps=0.5, N=30,out=self.stdfile, monitor_progress=True )
      
     

      rmax_indx = flex.min_index( self.ave_scores )
      self.best_rmax = self.rmax_list[ rmax_indx ]
      self.best_models = self.top_hits[rmax_indx]
      with  open(stdfile,"a") as log:
        print >>log, "   Best rmax found       :   %5.2f A"%self.best_rmax
      print  "   Best rmax found       :   %5.2f A"%self.best_rmax

      print >>self.out, "Best Result from Golden Section Search:", self.best_rmax
      with open(stdfile,"a") as log:
        self.show_result(self.best_rmax, self.best_models, self.scores[rmax_indx], codes, self.out, log)
      #self.show_result2(self.best_rmax, self.best_models, self.scores[rmax_indx], codes,self.stdfile)
      # import threading
      # t = []
      # t.append(threading.Thread(target=self.show_result()))
      # t.append(threading.Thread(target=self.show_result2()))
      # for t1 in t:
      #   t1.setDaemon(True)
      #   t1.start()
      #   t1.join()
      self.plot_intensity(self.best_rmax, self.best_models, self.scores[rmax_indx], self.coefs, codes, qmax=None,scales=self.scales[rmax_indx])
      #self.plot_intensity(self.best_rmax, self.best_models, self.scores[rmax_indx], self.coefs, codes, qmax=0.5,scales=self.scales[rmax_indx])
      self.print_rmax_profile( self.rmax_list, self.ave_scores )
      with open(self.outfilelog,"a") as f:
        print >> f, self.prefix+'.sta'
      self.summary( self.top_hits, self.scores, comment="----Statistics from Golden Section Scan----" )

      self.local_scan( int( self.best_rmax + 0.5) , local=5 )

    else:
      print
      log = open(stdfile,"a")
      print >>log, "   Not performing search in rmax. Using fixed value of rmax=%5.3f"%self.rmax
      log.close()
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
    self.update_model_error( rmax )
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
    z_model = zernike_model(self.nlm_array, self.data.q, rmax/self.fraction, self.nmax )
    top_hits, scores, scales = self.search( z_model, self.coefs, self.ntop, rmax )
    with open(self.stdfile,"a") as log:
      self.show_result( rmax, top_hits, scales, self.codes, self.out, log)
    return top_hits, scores, scales


  def get_rg(self, data, out=None):
    out = self.file
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
    fancy_weight_obj = error_model.two_stage_basic(self.data.q,self.data.i,self.data.s,0.02,0.75)
    w = fancy_weight_obj.error_model_smear(rmax,3)
    w = flex.sqrt(w)

    for coef in coefs:
      self.nn_array.load_coefs( self.nn, coef )
      i_cal = z_model.calc_intensity( self.nn_array )
      scale = i_cal[0]
      i_cal = i_cal/(scale/self.int_scale)
      scale = self.approx_scale(self.data.i,i_cal,w)
      if(self.smear):
        new_i = smear_data( i_cal, self.data.q, self.delta_q )
        score = flex.sum_sq( ( new_i - scale*self.data.i[1:-2] )/(scale*w[1:-2]) )
        #score = 1.0 - self.score_cc( new_i,self.data.i[1:-2], self.weights[1:-2] )
      else:
        score = flex.sum_sq( ( i_cal - self.data.i*scale )/(scale*w) )
        #score = 1.0 - self.score_cc( new_i,self.data.i,self.weights )

      scores.append( score )
      scales.append( scale )

    top_hits = self.tops( scores, ntop )
    top_scores = flex.double()
    top_scales = flex.double()
    for tt in top_hits:
      top_scores.append( scores[tt] )
      top_scales.append( scales[tt] )

    return top_hits, top_scores, top_scales


  def show_result( self, rmax, top_hits, scores, codes, out, out2): 
    print>>out, "RMAX:: ", rmax
    print>>out2, "RMAX: ", rmax
    for ii in range( top_hits.size() ):
      indx = top_hits[ii]
      print>>out,  "The top %s match is %s %.3f" %(str(ii+1), codes[indx], scores[ii])
      print>>out2,  "The top %s match is %s %.3f" %(str(ii+1), codes[indx], scores[ii])

      #print>>out,  "The top %s match is "%ii, codes[indx], scores[ii]
      #print>>out2,  "The top %s match is "%ii, codes[indx], scores[ii]

  # def show_result2( self, rmax, top_hits, scores, codes, outfile): 
  #   print outfile
  #   with open(outfile) as out:
  #     print>>out, "RMAX:: ", rmax
  #     for ii in range( top_hits.size() ):
  #       indx = top_hits[ii]
  #       print>>out,  "The top %s match is "%ii, codes[indx], scores[ii]

  def print_rmax_profile( self, rmax_list, ave_scores ):
    with open(self.outfilelog,"a") as f:
      print >>f ,self.prefix+".rmax_profile"
    out=open(self.prefix+".rmax_profile", 'w')
    for r, s in zip( rmax_list, ave_scores):
      print>>out, r, s
    out.close()

  def tops( self, x, n ):
    order = flex.sort_permutation( x )
    return order[0:n]

  def plot_intensity( self, rmax, top_hits, scores, coefs, codes, qmax=None, scales=None ):
    prefix = self.prefix
    if qmax is None:
      q_array = self.data.q
    else:
      q_array = flex.double( range(int((qmax-self.data.q[0])*100)) )/100.0 + self.data.q[0]
    z_model = zernike_model(self.nlm_array, q_array, rmax/self.fraction, self.nmax )

    for nn in range( top_hits.size() ):
      scale  = 1.0
      if scales is not None:
        scale = 1.0/scales[nn]
      oo = top_hits[nn]
      self.nn_array.load_coefs(self.nn, self.coefs[oo] )
      i_cal = z_model.calc_intensity( self.nn_array )
      if(self.smear):
        new_i = smear_data( i_cal, q_array, self.delta_q )
        for ii in range(1, new_i.size()+1):
          i_cal[ii] = new_i[ii-1]
      i_cal = i_cal / i_cal[0] * self.int_scale * scale
      filename = prefix+"_"+str(nn+1)+".iq"
      file = open(filename, 'w')
      print>>file, "# RMAX: "+str(rmax)+"  TOP: "+str(nn)+"  SCORE: "+str( scores[nn] )+"  PDB CODE: "+str(codes[oo])
      for qq, ii in zip( q_array, i_cal ):
        print>>file, qq, ii
      file.close()
      write_json(filename+".json", q_array, i_cal, self.data.i)


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

    clusters = hcluster.hcluster(self.cc_array, 0.8,self.stdfile)
    clusters.print_hclust(self.stdfile)
    out.close()
    tree_dot_file = self.prefix+".tree"
    clusters.print_neato(tree_dot_file)
    # generate image using neato
#    run_command('/sw/bin/neato -Tpng -o cluster.png '+tree_dot_file )
    #clusters.print_neato()
    self.clusters = clusters
    return

def get_rg(data, out=None):
    global stdfile
    out = stdfile
    msga = guinier_analyses.multi_step_rg_engine( data, out)
    return msga.median_rg, msga.median_io

def reduce_raw_data(raw_data, qmax, bandwidth, level=0.05, q_background=None):
    global stdfile
    
    log = open(stdfile,"a")
    print
    print >>log, " ====  Data reduction ==== "
    print
    print >>log, "  Preprocessing of data increases efficiency of shape retrieval procedure."
    print
    print >>log, "   -  Interpolation stepsize                           :  %4.3e"%bandwidth
    print >>log, "   -  Uniform density criteria:  level is set to       :  %4.3e"%level
    print >>log,"                                 maximum q to consider :  %4.3e"%qmax
    log.close()

    print  " ====  Data reduction ==== "
    print
    print  "  Preprocessing of data increases efficiency of shape retrieval procedure."
    print
    print  "   -  Interpolation stepsize                           :  %4.3e"%bandwidth
    print  "   -  Uniform density criteria:  level is set to       :  %4.3e"%level
    print "                                 maximum q to consider :  %4.3e"%qmax
  
    qmin_indx = flex.max_index( raw_data.i )
    qmin = raw_data.q[qmin_indx]
    new_data = get_q_array_uniform_body(raw_data, q_min=qmin, q_max=qmax, level=level )
    qmax = new_data.q[-1]
    if qmax > raw_data.q[-1]:
      qmax = raw_data.q[-1]
    log = open(stdfile,"a")
    print >>log,"      Resulting q range to use in  search:   q start   :  %4.3e"%qmin
    print >>log, "                                             q stop    :  %4.3e"%qmax
    log.close()
    print "      Resulting q range to use in  search:   q start   :  %4.3e"%qmin
    print  "                                             q stop    :  %4.3e"%qmax
  
    raw_q = raw_data.q[qmin_indx:]
    raw_i = raw_data.i[qmin_indx:]
    raw_s = raw_data.s[qmin_indx:]
    ### Take care of the background (set zero at very high q) ###
    if( q_background is not None):
      cutoff = flex.bool( raw_q > q_background )
      q_bk_indx = flex.last_index( cutoff, False )
      if( q_bk_indx < raw_q.size() ):
        bkgrd = flex.mean( raw_i[q_bk_indx:] )
        log = open(stdfile,"a")
        print >>log, "Background correction: I=I-background, where background=", bkgrd
        log.close()
        print  "Background correction: I=I-background, where background=", bkgrd

        raw_i = flex.abs( raw_i -bkgrd )

    q = flex.double( range(int( (qmax-qmin)/bandwidth )+1) )*bandwidth + qmin
    raw_data.i = flex.linear_interpolation( raw_q, raw_i, q )
    raw_data.s = flex.linear_interpolation( raw_q, raw_s, q )
    raw_data.q = q

    return raw_data

def process( pdb_files, nmax, rmax=50.0, fraction=0.9 ):
  pdb_models = []
  rmax_over_fraction = rmax/fraction
  shift = (rmax_over_fraction,rmax_over_fraction, rmax_over_fraction)
  for file in pdb_files:
    mom_obj, vox_obj, ipdb = pdb2zernike.zernike_moments( file, nmax=nmax, coef_out=False, calc_intensity=False)
    pdb_rmax = vox_obj.rmax()
    if( vox_obj.rmax() < rmax ):
      mom_obj, vox_obj, ipdb = pdb2zernike.zernike_moments( file, nmax=nmax, external_rmax=rmax, coef_out=False, calc_intensity=False)
    if(mom_obj is not None):
      if(len(pdb_models)==0):
        ref_nlm_array = mom_obj.moments()
        pdb_models.append( pdb_model( mom_obj.moments().coefs(), file, pdb_rmax, ipdb ) )
        ea=(0,0,0)
        write_pdb( file, vox_obj, ea, shift, ipdb )
      else:
        mov_nlm_array = mom_obj.moments()
        align_obj = fft_align.align(ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True)
        pdb_models.append( pdb_model(align_obj.moving_nlm.coefs(), file, pdb_rmax, ipdb) )
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
  global stdfile
  base=filename.split('.')[0]
  # base = base.split('/')
  # n=len(base)
  # base = base[n-1]
  log = open(stdfile,"a")
  print >>log , base
  log.close()
  print base
  # we need a more clever way tp write out pdb files
  out_pdb_name = base+"_sa.pdb"
  aligned_xyz = pdbinput.hierarchy.atoms().extract_xyz()
  center_xyz = aligned_xyz.mean()
  aligned_xyz = aligned_xyz - center_xyz + shift
  pdbinput.hierarchy.atoms().set_xyz( aligned_xyz )
  pdbinput.hierarchy.write_pdb_file( file_name=out_pdb_name, open_append=False)
  return out_pdb_name


def build_map( nmax, rmax, coefs, codes, model_indx, pdb_models,prefix, clusters=None, fract=0.9, type='.ccp4'):
  global stdfile
  global outfilelog
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

  with open(stdfile,"a") as log:
    print >>log, "Rank PDB_code cc (to the given model or the first model):"
  print  "Rank PDB_code cc (to the given model or the first model):"

  if( pdb_models is not None):
    ref_nlm_array.load_coefs( nlm, pdb_models[0].nlm_coef )
    fraction = 0
  else:
    c = coefs[model_indx[0] ][0:nlm_total]
    ave_c = c.deep_copy()
    ref_nlm_array.load_coefs( nlm, c )
    rank = 1
    if prefix!=None:
      filename = prefix+"m"+str(rank)+"_"+codes[ model_indx[0] ]
      #filename = os.path.join(prefix,"m"+str(rank)+"_"+codes[ model_indx[0] ])
    else:
      filename = "m"+str(rank)+"_"+codes[ model_indx[0] ]
      

    map_files.append( filename+type )
    fraction,map=write_map( zga, nlm, c,np_on_grid, rmax_over_fraction, filename, type=type )
    
    #level = flex.max(map)/3.0
    print "map in search pdb.py: \n"
    print map
    level = map.standard_deviation_of_the_sample()
    print "level in search pdb:  ",  level
    levels.append(level)
    top_cc.append( 1.0 )  # refering to itself
    if( scale_model ):
      c = get_moments_for_scaled_model( map, np_on_grid, grid, nmax, rmax, external_rmax )
    with open(stdfile,"a") as log:
      print >>log,rank, codes[ model_indx[0] ]
    print  rank, codes[ model_indx[0] ]

    aligned_coefs.append( c.deep_copy() )# save the aligned nlm coefs

  mean_frac = fraction
  mean_sqr_frac = fraction*fraction
  for ii in model_indx[rank:]:
    rank = rank + 1
    c = coefs[ ii ][0:nlm_total]
    mov_nlm_array.load_coefs( nlm, c )
    align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    new_c = align_obj.moving_nlm.coefs()
    if prefix!=None:
      #filename = os.path.join(prefix,"m"+str(rank)+"_"+codes[ii])
      filename = prefix+"m"+str(rank)+"_"+codes[ii]
      print "**********************************"
      print "outfilelog: ",outfilelog
      with open(outfilelog,"a") as f:
        print >> f, filename+".ccp4"
    else:
      filename = "m"+str(rank)+"_"+codes[ii]
      with open(outfilelog,"a") as f:
        print >> f,filename+".ccp4"
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
    with open(stdfile,"a") as log:
      print >>log, "%2d  %5s  %5.3f"%(rank, codes[ ii ], cc )
    print  "%2d  %5s  %5.3f"%(rank, codes[ ii ], cc )

    top_cc.append( cc )
    top_ids.append( codes[ii] )
    ave_c = ave_c + new_c
    aligned_coefs.append( new_c.deep_copy() )  # save the aligned nlm coefs
    mean_frac = mean_frac + fraction
    mean_sqr_frac = mean_sqr_frac + fraction*fraction

  sphere_volume = rmax_over_fraction**3.0*8.0  # cube with d=2.0*r
  mean_frac = mean_frac / ntop
  sigma_frac = smath.sqrt( mean_sqr_frac/ntop - mean_frac*mean_frac )
  
  with open(stdfile,"a") as log:
    print >>log, "Volume is ", mean_frac * sphere_volume, "+/-", sigma_frac*sphere_volume, "(A^3)"
  print  "Volume is ", mean_frac * sphere_volume, "+/-", sigma_frac*sphere_volume, "(A^3)"

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

  with open(stdfile,"a") as log:
    print >>log,  "cc. between Cluster average and PDB model"

  print   "cc. between Cluster average and PDB model"

  for node in clusters.nodes:
    ave_c = ave_c*0
    coefs_list = []
    for ii in node.leaf_eles:
      ave_c = ave_c + aligned_coefs[ii]
      cluster_ids[ii]=cluster_id
      coefs_list.append( aligned_coefs[ii] )
    ave_c = ave_c/len(node.leaf_eles)
    level_n = model_consistency.good2n( nmax, coefs_list, ave_c,threshold=0.90,outfile=stdfile)
    
    with open(stdfile,"a") as log:
      print >>log,"consistency level to order n: %d"%level_n

    print "consistency level to order n: %d"%level_n


    mov_nlm_array.load_coefs( nlm, ave_c )
    align_obj=fft_align.align( ref_nlm_array, mov_nlm_array, nmax=nmax, refine=True )
    cc = align_obj.get_cc()
    ave_cc.append(cc)

    with open(stdfile,"a") as log:
      print >>log, "cluster  # ", cluster_id, "cc=", cc
    print  "cluster  # ", cluster_id, "cc=", cc

    
    if prefix==None:
      filename = "ave_"+str(cluster_id)
      with open(outfilelog,"a") as f:
        print >> f, filename
    else:
      filename = prefix+"ave_"+str(cluster_id)

      with open(outfilelog, "a") as f:
        print >> f, filename+".ccp4"
      

    ave_maps.append( filename+type )
    fraction, map = write_map( zga, nlm, ave_c, np_on_grid, rmax_over_fraction, filename, type=type )
    ave_levels.append( flex.max(map)/3.0 )
    with open(stdfile,"a") as log:
      print >>log,"Averaged Model #%d Volume is %f (A^3)"%(cluster_id, fraction * sphere_volume)
    
    print "Averaged Model #%d Volume is %f (A^3)"%(cluster_id, fraction * sphere_volume)

    cluster_id = cluster_id+1
  
  # with open(stdfile,"a") as log:
  #   log.write("__END__")
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
  global stdfile
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/database/'
  with open(stdfile,"a") as log:
    print >> log, "\nATTENTION: Database path was set to : >>%s<<"%path
  
  print  "\nATTENTION: Database path was set to : >>%s<<"%path

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
  global stdfile
  global outfilelog
  
  targetpath_fromGUI = ''
  targetpath_fromGUI_file =  os.path.join(base_path,"targetpath_GUI.txt")
  if os.path.isfile(targetpath_fromGUI_file) and (os.stat(targetpath_fromGUI_file).st_size>0):
    with open(targetpath_fromGUI_file,"r") as f:
      targetpath_fromGUI = f.read().strip()

  if targetpath_fromGUI == '':
    stddir = "maps"
  else:
    tempfile = os.path.join(targetpath_fromGUI,"Shape_Search_Engine")
    stddir = os.path.join(tempfile,"maps")
  #stdfile = os.path.join(tempfile,"temp.txt")
  stdfile = os.path.join(os.path.split(sys.path[0])[0],"shapeup.txt")
  with open(stdfile,"w") as  f:
    f.truncate()

  outfilelog = os.path.join(os.path.split(sys.path[0])[0],"outfilelog_shapeup.txt")
  with open(outfilelog,"w") as f:
    f.truncate()
  
  t1 = time.time()
  with open(stdfile,"a") as outfile:
    params = get_input(args, master_params, "query", banner, help, outfile)
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
  scale_power = params.query.scale_power
  q_step = 1.0/200.0

  data = saxs_read_write.read_standard_ascii_qis( target_file )
  try:
    rg, io = get_rg(data)
  except:
    with open(stdfile,"a") as log:
      print >> log, "Guinier analysis failed, R_max is required"
      print >>log, "ATTENTION: dummy values for Rg and Io set"

      print  "Guinier analysis failed, R_max is required"
      print  "ATTENTION: dummy values for Rg and Io set"
    rg=50
    io=1

  qmax = params.query.qmax
  q_background = params.query.q_background
  #qmax = 0.44*smath.exp( -0.00023*rmax*rmax )
  ######### Interpolation ##########
  if (rmax is None):
    rmax=50
  bandwidth = min( q_step, smath.pi/2.0/rmax,data.q[1]-data.q[0] )
  data = reduce_raw_data( data, qmax, bandwidth, q_background=q_background,level=params.query.q_level )
  ###### END of Interpolation ##########
  with open(stdfile,"a") as log:
    print >>log, " ==== Reading in shape database ==== "

  print  " ==== Reading in shape database ==== "
  begin_time = time.time()
  nn_coefs, codes, rmaxs = read_pickle(dbpath, dbprefix)
  ready_time = time.time()
  delta_time = ready_time-begin_time
  print
  with open(stdfile,"a") as log:
    print >>log, "   Done reading database with %i entries in %5.4e seconds"%(len(codes),delta_time)

    print >>log, " ==== Shape retrieval ==== "
    print >>log, "   Constructing shape retrieval object"

  print  "   Done reading database with %i entries in %5.4e seconds"%(len(codes),delta_time)

  print  " ==== Shape retrieval ==== "
  print  "   Constructing shape retrieval object"
  shapes = intoshape( data, rg=rg, io=io, nmax=nmax, rmax=rmax, scan=scan, fraction=fraction, smear=smear, prefix=params.query.prefix, weight=weight, delta_q=delta_q,scale_power=scale_power)
  with open(stdfile,"a") as log:
    print >>log, "   Shape search  .... "

  print  "   Shape search  .... "
  shapes.lookup(nn_coefs, codes, ntop)

  nlm_coefs = read_nlm(dbpath, dbprefix)
  shapes.pair_align( nlm_coefs, params.query.calc_cc )

  pdb_models=None
  if( len(pdb_files) > 0 ):
    pdb_models = process(pdb_files,nmax, rmax=shapes.best_rmax, fraction=fraction)

  if(params.query.buildmap):
    top_cc, top_ids, map_files, levels, cluster_ids, ave_maps, ave_levels, ave_cc = build_map( nmax, shapes.best_rmax, nlm_coefs, codes, shapes.best_models, pdb_models, clusters=shapes.clusters,fract=fraction,prefix=params.query.prefix)
                   # need to use rmax/fraction to get right size of box
    #build_pymol_script.write_pymol_scripts(maps=map_files,levels=levels,root_name=stddir)

    build_pymol_script.write_pymol_shapeup(maps=map_files,root_name=stddir)
    pdb_out_name=None
    if( pdb_models is not None ):
      pdb_out_name = pdb_files[0].split('.')[0]+'_sa.pdb'
    #generate_html.generate_jmol_html(ave_maps, ave_cc, ave_levels, map_files, top_cc, levels, cluster_ids, 'models.html', pdb=pdb_out_name)
    if(len(pdb_files) > 0 ):
      with open(params.query.prefix+"_cc2pdb.dat", 'w') as out:
        print >> out, "Correlation coefficients of retrieved shapes vs input model"
        for cc,id in zip(top_cc,top_ids):
          print>>out,"Code: %5s    CC: %5.1f  "%(id,100*cc)
        print>>out, "mean: %8.5f"%flex.mean(top_cc)

      with open(stdfile,"a") as log:
        print >>log, "Compared to the PDB model (%s)"%pdb_models[0].filename
        print >>log, "mean cc: %8.5f"%flex.mean(top_cc)
        print >>log, "first cc: %8.5f"%top_cc[0]
        print >>log, "best cc: %8.5f"%flex.max(top_cc)
        print >>log, "worst cc: %8.5f"%flex.min(top_cc)
      

      print  "Compared to the PDB model (%s)"%pdb_models[0].filename
      print  "mean cc: %8.5f"%flex.mean(top_cc)
      print  "first cc: %8.5f"%top_cc[0]
      print  "best cc: %8.5f"%flex.max(top_cc)
      print  "worst cc: %8.5f"%flex.min(top_cc)

      with open(stdfile,"a") as log:
        print >>log, "Rmax: estimated vs PDB", shapes.best_rmax, pdb_models[0].rmax
      
      print  "Rmax: estimated vs PDB", shapes.best_rmax, pdb_models[0].rmax

  t2 = time.time()
  with open(stdfile,"a") as log:
    print >>log,"total time used: ", t2-t1, "(seconds)"
  print  "total time used", t2-t1, "(seconds)"

  with open(stdfile,"a") as log:
    log.write("__END__")

if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
