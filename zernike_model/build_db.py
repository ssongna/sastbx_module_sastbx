import os, sys, time
from sastbx.zernike_model import pdb2zernike
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx.option_parser import option_parser
from libtbx import easy_pickle

from multiprocessing import Pool

from sastbx.interface import get_input

master_params = iotbx.phil.parse("""\
db{
  path = None
  .type=path
  .help="path of pdb files"

  nmax=20
  .type=int
  .help="maximum order of zernike expansion"

  fix_dx=True
  .type=bool
  .help="Whether keeping default dx=0.7A or not"

  np=50
  .type=int
  .help="number of point covering [0,1]"

  nprocess=4
  .type=int
  .help="number of processes"

  prefix="myDB"
  .type=path
  .help="the prefix of pickle file names"

}
""")

banner = "--------------Build Zernike Moment DataBase----------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\n  Usage: libtbx.python build_db.py path=path nmax=nmax fix_dx=True/False"
  print >> out, "\n  Attention: \n  The path should be a directory that contains ONLY PDB files\n\n"


def read(path='./'):
  files = os.listdir(path)
  return files

def get_results( args ): # args = [path, filename, nmax, np, fix_dx]
  code = args[1].split('.')[0]
  file = args[0]+args[1]
  nmax=args[2]
  np = args[3]
  fix_dx = args[4]

  mom_obj, vox_obj, ipdb = pdb2zernike.zernike_moments( file, nmax=nmax, np=np, fix_dx=fix_dx, coef_out=False, calc_intensity=False )
  if( mom_obj is None):
    print code, "not processed, check the file"
    return None
  print code, "processed."
  return (mom_obj.moments().coefs(), mom_obj.fnn().coefs(), code, vox_obj.rmax() )


def run(args):
  params = get_input(args, master_params, "db", banner, help)
  if params is None:
    exit()
  path = params.db.path+"/"
  nmax=params.db.nmax
  np = params.db.np
  fix_dx = params.db.fix_dx
  prefix = params.db.prefix
  NPROCESS = params.db.nprocess

  files = read(path)
  nlm_coefs = []
  nn_coefs = []
  codes = []
  rmaxs = []
  inputs = []
  for file in files:
    inputs.append( (path, file, nmax, np, fix_dx) )



  # print 'NPROCESS', NPROCESS
  # p = Pool( NPROCESS )
  # CHUNKSIZE = len(inputs)/2/NPROCESS
  # print "CHUNKSIZE", CHUNKSIZE
  # results = p.map( get_results, inputs )
  results = []
  count =0
  for i in inputs:
    count = count+1
    print "count",count
    print i
    temp = get_results(i)
    results.append(temp)

  for result in results:
    if( result is None ):
      continue
    nlm = result[0]
    nn = result[1]
    code = result[2]
    rmax = result[3]
    codes.append( code )
    rmaxs.append( rmax )
    nlm_coefs.append( nlm.deep_copy() )
    nn_coefs.append( nn.deep_copy() )

  easy_pickle.dump(prefix+".nlm", nlm_coefs)
  easy_pickle.dump(prefix+".nn", nn_coefs)
  easy_pickle.dump(prefix+".rmax", rmaxs)
  easy_pickle.dump(prefix+".codes", codes)


if __name__ == "__main__":
  t1 = time.time()
  args = sys.argv[1:]
  run(args)
  t2 = time.time()
  print "time used: ", t2-t1
