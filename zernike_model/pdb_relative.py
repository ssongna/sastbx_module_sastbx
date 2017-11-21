from libtbx import easy_pickle
from scitbx.array_family import flex
from scitbx import math
from scitbx.math import zernike_align_fft as fft_align
from sastbx.zernike_model import zernike_model
from scitbx.golden_section_search import gss

import os, sys

def set_default_db_path():
  import libtbx
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/database/'
  print "\nATTENTION: Database path was set to : >>%s<<"%path
  return path


prefix=set_default_db_path()

trivial_nlm_array = math.nlm_array(10)

class optimize_r(object):
  def __init__(self, nn_array, ref_int, q_array, nmax=10):
    self.nn_array = nn_array
    self.ref_int = ref_int
    self.q_array = q_array
    self.nmax=nmax

  def target(self, r):
    z_model=zernike_model( trivial_nlm_array, self.q_array, r, self.nmax )
    calc_i = z_model.calc_intensity( self.nn_array )
    calc_i = calc_i/calc_i[0]
    return flex.sum_sq(self.ref_int-calc_i)

def load_data():
  codes=easy_pickle.load(prefix+'pisa.codes')
  codes = flex.std_string( codes )
  moments=easy_pickle.load(prefix+'pisa.nlm')
  return codes, moments

def put_intensity( z_model, q_array, nn_array, out_name,ref=None):
  calc_i = z_model.calc_intensity( nn_array )
  calc_i = calc_i/calc_i[0]
  if ref is None:
    #write_file(q_array, calc_i, out_name)
    return calc_i
  else:
    #write_file(q_array, calc_i/ref, out_name)
    return calc_i/ref


def write_file(q_array, i_array, out_name):
  output=open(out_name, 'w')
  for qq,ii in zip( q_array, i_array):
    print>>output, qq,ii
  output.close()


def find_relatives( ids, cc_min, cc_max, rmax, codes, moments, nmax=10 ):
  indices = flex.int()
  idlist = open('id_list.txt','r')
  for id in idlist:
    id = id[0:4]
    indices.append( flex.first_index( codes, id ) )
  r_max = easy_pickle.load(prefix+'pisa.rmax')
  nns = easy_pickle.load(prefix+'pisa.nn')
  nn_array = math.nl_array(nmax)
  nn_indx = nn_array.nl()
  nn_total = nn_indx.size()
  q_array = flex.double( range(501) )/2000.0

  ref_nlm_array    = math.nlm_array(nmax)
  target_nlm_array = math.nlm_array(nmax)
  nlm = ref_nlm_array.nlm()
  coef_size=nlm.size()
  all_indices = range( codes.size() )

  small_q_array=flex.double(range(51))/300.0
  mean = []
  sig = []
  for indx in indices:
    print indx
    #rmax = 50.0 #r_max[indx]
    ref_coef = moments[indx]
    ref_nlm_array.load_coefs( nlm, ref_coef[0:coef_size] )
    z_model = zernike_model( ref_nlm_array, q_array, rmax, nmax )
    out_name=codes[indx]+"_.qi"
    nn_array.load_coefs( nn_indx, nns[indx][0:nn_total] )
    ref_int = put_intensity( z_model, q_array, nn_array, out_name)
    mean_r = ref_int*0.0
    sig_r = ref_int*0.0
    small_z_model = zernike_model( ref_nlm_array, small_q_array, rmax, nmax )
    small_ref_int = small_z_model.calc_intensity( nn_array )
    small_ref_int = small_ref_int/small_ref_int[0]
    N = 0.0
    for coef, ii in zip(moments, all_indices) :
      if N > 25: break
      target_nlm_array.load_coefs( nlm, coef[0:coef_size] )
      align_obj = fft_align.align( ref_nlm_array, target_nlm_array, nmax=nmax, topn=10 ,refine=False)
      cc = align_obj.get_cc()
      if (cc>=cc_min and cc<=cc_max):
        N+=1
        nn_array.load_coefs( nn_indx, nns[ii][0:nn_total] )
        opt_r_obj = optimize_r(nn_array, small_ref_int, small_q_array, nmax)
        opt_r=gss( opt_r_obj.target, rmax*0.8, rmax*1.2)
        z_model = zernike_model( ref_nlm_array, q_array, opt_r, nmax )
        out_name=codes[indx]+"_"+codes[ii]+".qi.rel"
        mod_int = put_intensity( z_model, q_array, nn_array, out_name,ref_int)
        out_name=codes[indx]+"_"+codes[ii]+".qi"
        put_intensity( z_model, q_array, nn_array, out_name)
        mod_int = mod_int -  1.0
        mean_r += mod_int
        sig_r += mod_int*mod_int
        print ii, cc, codes[ii], opt_r
    if N>3:
      mean_r /= N
      sig_r = sig_r/N - mean_r*mean_r
      mean.append(mean_r)
      sig.append(sig_r)

  N = len(mean)
  if N> 0:
    mean_r = mean[0]*0.0
    s_r    = mean[0]*0.0
    for uu in range(N):
      mean_r+=mean[uu]
      s_r += sig[uu]
    mean_r /= N
    s_r /= N
    s_r = flex.sqrt(s_r)
    f = open('q_m_s_%s.dat'%rmax,'w')
    for q,m,s in zip(q_array,mean_r,s_r):
      print >> f,q,m,s


def run(cc_min, cc_max, rmax, ids):
  codes, moments = load_data()
  find_relatives( ids, cc_min, cc_max, rmax, codes, moments)

if __name__ == "__main__":
  args = sys.argv[1:]
  if(len(args)<3):
    print "\n  Usages:\n"
    print "    find_relative.py min_cc max_cc pdb_ids \n"
    exit()

  cc_min = float(args[0])
  cc_max = float(args[1])
  rmax = float(args[2])
  ids    = args[3:]
  run(cc_min, cc_max, rmax, ids)
