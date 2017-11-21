from sastbx.zernike_model import pdb2zernike
from scitbx.array_family import flex
from sastbx.fXS import zernike_moment_variants
import os, sys

def tst(filename):
  nmax=20
  np=30
  fix_dx=True
  prefix=tst

  q_array = flex.double(range(50))/100.0
  mom_obj, vox_obj, pdb = pdb2zernike.zernike_moments( filename, nmax=nmax, np=np, fix_dx=fix_dx, coef_out=False, calc_intensity=True)
  c_nlm = mom_obj.moments()
  rmax = vox_obj.rmax()/0.9
  znk_mom_variants = zernike_moment_variants( c_nlm, q_array, rmax, nmax )

  this_blq=flex.double()
  for ii in range(50):
    real_coef = znk_mom_variants.get_real_coef(ii)
    print q_array[ii],
    for cc in real_coef:
      print cc,
      this_blq.append(cc)
    print
#  all_blq = znk_mom_variants.get_all_blq()
#  for t0,t1 in zip(this_blq,all_blq):
#    assert abs(t0-t1)<1e-6
  

def print_help():
  print "Usage: libtbx.python test_znk_fXS.py pdb_file.name"
  print "output: q, B0, B1, B2, B3, ....."

if __name__ == "__main__":
  args = sys.argv
  if(len(args) <2):
    print_help()
    exit()
  tst(sys.argv[1])
