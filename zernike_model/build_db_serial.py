import os, sys, time
from sastbx.zernike_model import pdb2zernike
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx.option_parser import option_parser
from libtbx import easy_pickle

from sastbx.interface import get_input

master_params = iotbx.phil.parse("""\
db{
  path = None
  .type=path
  .help="the experimental intensity profile"

  qmax = 0.3
  .type=float
  .help="maximum q value, for which the intensity to be evaluated"

  nmax=20
  .type=int
  .help="maximum order of zernike expansion"

  fix_dx=True
  .type=bool
  .help="Whether keeping default dx=0.7A or not"

  np=50
  .type=int
  .help="number of point covering [0,1]"

  prefix="myDB"
  .type=path
  .help="the prefix of pickle file names"

}
""")

banner = "--------------Build Zernike Moment DataBase----------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\n  Usage: libtbx.python build_db.py path=path nmax=nmax np=np_of_point fix_dx=True/False"
  print >> out, "\n  Attention: \n  The path should be a directory that contains ONLY PDB files\n\n"


def read(path='./'):
  files = os.listdir(path)
  return files

def run(args):
  params = get_input(args, master_params, "db", banner, help)
  if params is None:
    exit()
  path = params.db.path+"/"
  nmax=params.db.nmax
  np = params.db.np
  fix_dx = params.db.fix_dx
  prefix = params.db.prefix


  files = read(path)
  nlm_coefs = []
  nn_coefs = []
  codes = []
  rmax = []
  for file in files:
    code = file.split('\n')[0].split('.')[0]
    file = path+file
    mom_obj, vox_obj, pdb = pdb2zernike.zernike_moments( file, nmax=nmax, np=np, fix_dx=fix_dx, coef_out=False, calc_intensity=False )
    if(mom_obj is None):
      print code, "NOT processed, please check the file"
      continue
    codes.append( code )
    rmax.append( vox_obj.rmax() )
    nlm_coefs.append( mom_obj.moments().coefs().deep_copy() )
    nn_coefs.append( mom_obj.fnn().coefs().deep_copy() )
    print code, "processed."

  easy_pickle.dump(prefix+".nlm", nlm_coefs)
  easy_pickle.dump(prefix+".nn", nn_coefs)
  easy_pickle.dump(prefix+".rmax", rmax)
  easy_pickle.dump(prefix+".codes", codes)


if __name__ == "__main__":
  t1 = time.time()
  args = sys.argv[1:]
  run(args)
  t2 = time.time()
  print "time used: ", t2-t1
