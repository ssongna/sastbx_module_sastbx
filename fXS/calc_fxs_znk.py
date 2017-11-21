from sastbx.zernike_model import pdb2zernike
from sastbx.fXS import fxs_tools
from scitbx.array_family import flex
from sastbx.fXS import zernike_moment_variants
from sastbx.interface import get_input
import iotbx.phil
import os, sys

master_params = iotbx.phil.parse("""\
fxs_znk{
  nmax = 20
  .type=int
  .help="maximum epxansion order"

  lmax = None
  .type=int
  .help="maximum epxansion order"

  pdb = None
  .type = path
  .help = "pdb format model"

  np_on_grid = 30
  .type=int
  .help = "number of points on the grid covering [0,1]"

  q_start = 0.0
  .type = float
  .help = "start of the Q-array"

  q_stop = 0.5
  .type = float
  .help = "end of the Q-array"

  n_step = 51
  .type = int
  .help = "number of bins"

  output= "blq.dat"
  .type = path
  .help = "output filename"

  fix_dx = True
  .type = bool
  .help = "Using fixed dx when calculating zernike moments, where dx=0.7A"
}
""")

banner = """
============================================================================
      fXS auto-correlation lendedre expansion coefs calculation
                using Zernike Polynomial model
============================================================================
 """

def run(args, log=sys.stdout):
  params = get_input( args, master_params, "fxs_znk", banner, print_help )
  if (params is None):
    exit()
  nmax=params.fxs_znk.nmax
  lmax=params.fxs_znk.lmax
  if (lmax is None): lmax=nmax
  np_on_grid=params.fxs_znk.np_on_grid
  filename =params.fxs_znk.pdb
  output = params.fxs_znk.output
  fix_dx=params.fxs_znk.fix_dx
  q_array = None
  if q_array is None:
     q_array = params.fxs_znk.q_start +  \
               (params.fxs_znk.q_stop-params.fxs_znk.q_start) \
               *flex.double( range(params.fxs_znk.n_step) )/(params.fxs_znk.n_step-1)

  mom_obj, vox_obj, pdb = pdb2zernike.zernike_moments( filename, nmax=nmax, \
                          np=np_on_grid, fix_dx=fix_dx, coef_out=False, calc_intensity=True)
  c_nlm = mom_obj.moments()
  rmax = vox_obj.rmax()/0.9
  znk_mom_variants = zernike_moment_variants( c_nlm, q_array, rmax, nmax, lmax )

  out = open(output,'w')
  this_blq = znk_mom_variants.get_all_blq2()
  blq_data = fxs_tools.blq_data(q_array,this_blq,lmax )
  blq_data.print_out(out=out)
  out.close()


def print_help(out):
  print>>out, "Usage: libtbx.python calc_fxs_znk.py pdb=pdb_file.name"
  print>>out, "output: q, B0, B2, B4, B6, ..... (odd orders are Zero)"

if __name__ == "__main__":
  run(sys.argv[1:])
