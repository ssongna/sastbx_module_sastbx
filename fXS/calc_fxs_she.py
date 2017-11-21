from sastbx.intensity import she
from sastbx.fXS import fxs_tools
from scitbx.array_family import flex
from sastbx.interface import get_input
import iotbx.phil
import os, sys

master_params = iotbx.phil.parse("""\
fxs_she{
  lmax = 15
  .type=int
  .help="maximum epxansion order"

  pdb = None
  .type = path
  .help = "pdb format model"

  q_start = 0.0
  .type = float
  .help = "start of the Q-array"

  q_stop = 0.5
  .type = float
  .help = "end of the Q-array"

  n_step = 51
  .type = int
  .help = "number of bins"

  prefix = "output"
  .type = path
  .help = "output filename"

  print_c2 = false
  .type = bool
  .help = "output the auto-correlation curve"

}
""")

banner = """
============================================================================
      fXS auto-correlation lendedre expansion coefs calculation
                using Spherical Harmonics Expansion model
============================================================================
 """

def unpack( correlation, n_phi, n_q ):
  new_co = []
  index = 0
  for pp in range(n_phi):
    tmp_co = []
    for qq in range( n_q ):
      tmp_co_q = []
      for qqq in range(qq+1):
        tmp_co_q.append( correlation[index] )
        index = index + 1
      tmp_co.append( tmp_co_q )
    new_co.append( tmp_co )
  return new_co

def print_correlation( out_correlation, correlation, q_array, N_phi):
  print >>out_correlation, "# phi, q, correlation"
  two_pi = smath.pi*2.0
  ii = 0
  for qq in q_array:
    for pp in range(N_phi):
      phi = pp*two_pi/N_phi
      print >>out_correlation, qq*smath.cos(phi),qq*smath.sin(phi), correlation[pp][ii][ii]
    ii=ii+1

def run(args, log=sys.stdout):
  params = get_input( args, master_params, "fxs_she", banner, print_help )
  if (params is None):
    exit()
  lmax=params.fxs_she.lmax
  pdb_file =params.fxs_she.pdb
  prefix = params.fxs_she.prefix

  q_array = None
  if q_array is None:
     q_array = params.fxs_she.q_start +  \
               (params.fxs_she.q_stop-params.fxs_she.q_start) \
               *flex.double( range(params.fxs_she.n_step) )/(params.fxs_she.n_step-1)

  she_obj = she.she(pdb_file, q_array=q_array, max_L=lmax)
  N_phi = 51
  phi_array = flex.double( range(N_phi) ) / 25.0 * 3.1416
  she_obj.engine.calc_spatial_correlation(phi_array)  # expansion coefficient calculation is built in this call
  if( params.fxs_she.print_c2 ):
    correlation = she_obj.engine.get_spatial_correlation()
    correlation = unpack( correlation, N_phi, q_array.size() )
    out_correlation = open(prefix+'.cor','w')
    print_correlation( out_correlation, correlation, q_array, N_phi )
    close(out_correlation)

  out = open(prefix+'.blq','w')
  this_blq = she_obj.get_all_blq()
  blq_data = fxs_tools.blq_data(q_array,this_blq,lmax)
  blq_data.print_out(out=out)
  out.close()


def print_help(out):
  print>>out, "Usage: libtbx.python calc_fxs_she.py pdb=pdb_file.name"
  print>>out, "output: q, B0, B1, B2, B3, ....."

if __name__ == "__main__":
  run(sys.argv[1:])
