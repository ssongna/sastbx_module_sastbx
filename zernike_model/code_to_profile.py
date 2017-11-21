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
  qmax = 0.20
  .type=float
  .help="maximum q value"

  delta_q = 0.002
  .type=float
  .help = "step in q"

  fraction = 0.9
  .type=float

  rmax = None
  .type=float
  .help="rmax of the molecule"

  nmax = 10
  .type=int
  .help="Expansion order"

  dbpath=None
  .type=path
  .help="the directory of database file, i.e., the pickle files"

  db_choice = *pisa piqsi allpdb user
  .type = choice
  .help = "Data base name"

  db_user_prefix="mydb"
  .type=path
  .help="the prefix of database filename"

  code=None
  .type=str
  .help="PID ID from database"

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

def get_nn_array_for_code(code, nns,codes,rmaxs,nmax=10):
  this_one = codes.index(code)
  coefs =  nns[this_one]
  nn_array = math.nl_array(nmax)
  nn = nn_array.nl()
  nn_size = nn.size()
  coefs = coefs[0:nn_size]
  nn_array.load_coefs( nn, coefs )
  print code,rmaxs[this_one]
  for c in coefs:
    print c
  return nn_array



def get_profile(nn_array,rmax,fraction=0.9,nmax=10, qmax=0.25,qstep=0.002):
  nlm_array =  math.nlm_array(nmax)
  n = int(qmax/qstep+0.5)
  q_array = flex.double( range(n+1) ) / float(n) * qmax
  z_model = zernike_model(nlm_array, q_array, rmax/fraction, nmax )
  result = z_model.calc_intensity( nn_array )
  return q_array,result



def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, " NO HELP YET"

def set_default_db_path():
  import libtbx
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/database/'
  print "\nATTENTION: database path was set to default:"
  print ">>>>  dbpath = ", path, "  <<<<"
  return path



def run(args):
  t1 = time.time()
  params = get_input(args, master_params, "query", banner, help)
  if( params is None):
    exit()
  rmax = params.query.rmax
  nmax = params.query.nmax
  code = params.query.code
  dbpath = params.query.dbpath
  db_choice = params.query.db_choice
  delta_q = params.query.delta_q
  if( db_choice == "user" ):
    dbprefix = params.query.db_user_prefix
  else:
    dbprefix = db_choice

  if (dbpath is None):
    dbpath = set_default_db_path()
  fraction = params.query.fraction


  nn_coefs, codes, rmaxs = read_pickle(dbpath, dbprefix)
  nn_array = get_nn_array_for_code(code, nn_coefs,codes,rmaxs, nmax=nmax)
  intensity = get_profile(nn_array,rmax,fraction=0.9,nmax=nmax)


if __name__=="__main__":
  args = sys.argv[1:]
  run(args)
