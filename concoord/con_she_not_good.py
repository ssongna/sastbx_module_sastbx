from scitbx.array_family import flex
from mmtbx.monomer_library import server, pdb_interpretation
from iotbx import pdb
import iotbx.phil
from sastbx.intensity import she
from sastbx import intensity
from sastbx.interface import get_input
from sastbx.data_reduction import saxs_read_write
from stdlib import math as smath
import libtbx.env_config
from libtbx.test_utils import run_command
import os, sys, time, random
from multiprocessing import Pool

master_params = iotbx.phil.parse("""\
concoord{
  start_pdb = None
  .type=path
  .help="the structure to be refined"

  target_I = None
  .type=path
  .help="the scattering profile"

  n_struct= 50
  .type=int
  .help="number of structures generated each run of disco"

  prefix = "prefix"
  .type=str
  .help="prefix of the output pdb files"

  n_process= 1
  .type=int
  .help="number of available processors"
}
""")

banner = """
============================================================================
                  Concoord-based Structure Refinement
============================================================================
"""

def print_help(out):
  print """\nThe program utilizes CONCOORD to generate diverse conformations,
        which are subsequently analyzed using SHE to get I(q), the Chi-score
        is compared to the supplied experimental data"""

class concoord(object):
  def __init__(self, pdb_file, target, nstruct=500,
               np=50,max_np=100, prefix='prefix'):
    self.pdb_file = pdb_file
    self.obs = saxs_read_write.read_standard_ascii_qis(target)
    if(self.obs.q.size() > max_np):
      self.obs= self.reduction(self.obs) # reduce number_of_point in q-array
    self.she_obj = she.she(pdb_file,self.obs.q)
      # More options are available, see line #10 for class she definition ##
    self.run_concoord(nstruct, prefix=prefix)
    self.files, self.scores = self.compute_score(nstruct, prefix=prefix)
    self.min_indx = flex.min_index( self.scores )
    self.min_file = self.files[ self.min_indx ]
    self.min_score= self.scores[ self.min_indx ]

  def get_results(self):
    return self.files, self.scores

  def reduction(self, data, np=50):
    qmin=max( data.q[0]-0.01, 0.0)
    qmax=data.q[-1]+0.01
    new_q = flex.double(range(np))/float(np)*(qmax-qmin)+qmin
    data.i = flex.linear_interpolation( data.q, data.i, new_q )
    data.s = flex.linear_interpolation( data.q, data.s, new_q )
    data.q = new_q
    return data

  def run_concoord(self, nstruct, prefix='prefix'):
    env = libtbx.env_config.unpickle()
    self.dist= env.find_dist_path('sastbx')+'/concoord/mydist.csh'
    self.disco= env.find_dist_path('sastbx')+'/concoord/mydisco.csh'
    # initialize distance constraints
    op_file = prefix+".pdb"
    og_file = prefix+".gro"
    od_file = prefix+".dat"
    run_command(self.dist+" -p %s -op %s -og %s -od %s >& %s.log0"%(self.pdb_file, op_file, og_file, od_file, prefix ) )
    # run disco to generate nstruct conformations
    #print "running:\n", self.disco+" -p %s -n %d -op %s>%s.log &"%(self.pdb_file,nstruct,prefix,prefix)
    random_number = random.randint(1,741265)
    run_command(self.disco+" -p %s -d %s -n %d -op %s -s %d >& %s.log &"%(op_file,od_file,nstruct,prefix,random_number,prefix) )


  def compute_score(self, nstruct, prefix):
    new_indx = flex.int(self.she_obj.natom)
    files = []
    results = flex.double()
    for ii in range(1,nstruct+1):
      pdb_file_name = prefix+str(ii)+'.pdb'
      while( not os.path.exists(pdb_file_name) ): time.sleep(0.02)
      time.sleep(0.05)
      new_xyz = extract_xyz(pdb_file_name)
      self.she_obj.engine.update_coord(new_xyz,new_indx)
      new_i = self.she_obj.engine.I()
      s,o = she.linear_fit(new_i, self.obs.i, self.obs.s)
      chi2= flex.mean( flex.pow2( (self.obs.i-(s*new_i+o))  /self.obs.s ))
      files.append( pdb_file_name )
      results.append( chi2 )
    return files, results

def extract_xyz(pdb_file_name):
    pdbi = pdb.hierarchy.input(file_name=pdb_file_name)
    xyz = pdbi.hierarchy.atoms().extract_xyz()
    return xyz

def get_concoord(para):
  pdb_file = para[0]
  target   = para[1]
  n_struct = para[2]
  prefix   = para[3]
  con_obj = concoord(pdb_file, target, nstruct=n_struct,
               np=50,max_np=100, prefix=prefix)
  return [con_obj.get_results(), con_obj.min_indx]

def run(args,log=sys.stdout):
  params = get_input( args, master_params, "concoord", banner, print_help)
  if (params is None): exit()
  pdb_file      = params.concoord.start_pdb
  target        = params.concoord.target_I
  n_struct      = params.concoord.n_struct
  prefix        = params.concoord.prefix
  nprocess      = params.concoord.n_process

  p = Pool(nprocess)

  inputs = []
  for ii in range(nprocess):
    sub_prefix = prefix + '_' + str(ii) + '_'
    inputs.append([ pdb_file, target, n_struct, sub_prefix] )

  results = p.map( get_concoord, inputs)
  p.close()
  p.join()

  for r in results:
    min_indx = r[1]
    print r[0][0][min_indx], r[0][1][min_indx]



if __name__ == "__main__":
  run(sys.argv[1:])
