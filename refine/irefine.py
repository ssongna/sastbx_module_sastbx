import sys,os
from scitbx import simplex
from scitbx.array_family import flex
from sastbx import refine
from sastbx.refine.PDB import PDB
from math import exp,log,sqrt
import random,time
from iotbx import pdb
from mmtbx.command_line import geometry_minimization
from mmtbx.utils import process_pdb_file_srv
from sastbx.intensity import she
from sastbx.data_reduction import saxs_read_write
from sastbx.interface import get_input
import iotbx.phil


master_params = iotbx.phil.parse("""\
pr_ref{
  start_pdb = None
  .type=path
  .help="the structure to be refined"

  target_I = None
  .type=path
  .help="the pr derived from scattering profile"

  method = *rtb ca
  .type = choice
  .help = "Method to calculate Normal Modes"

  n_modes = 20
  .type=int
  .help="number of normal modes used to refine the structure"

  max_rmsd = 0.2
  .type=float
  .help="maximum rmsd allowed for each iteration"

  backbone_force = 10.0
  .type=float
  .help="force constant of backbone (relative to the rest)"

  weight=i* s
  .type=choice
  .help="weighting scheme for the target function"

  prefix="I_ref"
  .type=str
  .help="prefix of the output pdb file"

}

""")

banner = """
============================================================================
       Iterative Normal Mode Perturbation for Structure Refinement
============================================================================
"""



class nmref(object):
  def __init__(self,start_pdb, target_I, ntotal,nmodes,max_rmsd,backbone_scale,prefix, weight='i', method='rtb', log='tmp.log'):
    self.counter = 0
    self.nmode_init = ntotal
    self.nmodes = 3 #nmodes
    self.method=method
    self.topn=3
    self.Niter=0
    self.modes=flex.int(range(self.nmode_init))+7
    self.cutoff=12
    self.weighted = True
    self.log = open(log, 'w')
    pdb_inp = pdb.input(file_name=start_pdb)
    crystal_symmetry = pdb_inp.xray_structure_simple().\
        cubic_unit_cell_around_centered_scatterers(
        buffer_size = 10).crystal_symmetry()
#    uc=cctbx.uctbx.unit_cell("300,300,300,90,90,90")
#    crystal_symmetry=cctbx.crystal.symmetry(uc, 'P1')
    self.pdb_processor = process_pdb_file_srv( crystal_symmetry=crystal_symmetry )

    self.expt = saxs_read_write.read_standard_ascii_qis(target_I)
    self.q = self.expt.q
    self.expt_I = self.expt.i
    self.expt_s = self.expt.s
    if(self.q.size() > 50):
      self.q = self.interpolation(self.q, n_pts=50)
      self.expt_I = flex.linear_interpolation( self.expt.q, self.expt.i, self.q)
      self.expt_s = flex.linear_interpolation( self.expt.q, self.expt.s, self.q)

    if( weight=='i'):
      self.expt_s = self.expt_I

    self.time_nm =0
    self.time_she=0
    start_name=start_pdb
    self.pdb = PDB(start_name, method=self.method)
    self.she_engine=she.she(start_name,self.q)
    self.natom = self.pdb.natm
    self.scale_factor=backbone_scale
    time1=time.time()
    self.nmode=self.pdb.Hessian(self.cutoff,self.nmode_init,self.scale_factor)
    self.time_nm += (time.time() - time1 )

    self.root=prefix
    self.drmsd = max_rmsd
    self.Rmax2 = self.natom * (self.drmsd) **2.0
    self.step_size= sqrt(self.Rmax2 /self.nmodes)*5.0

    self.new_indx=flex.int(range(self.natom))
    self.stop = False
    self.minscore=1e20
    self.minDev = 0  #minimum deviations of refined structures, compared to refined structure from the previous step
    self.optNum = 10 #number of iterations between geometry optimization
    self.iterate()
    self.log.close()

    print "time used for NM : %d"%self.time_nm
    print "time used for she: %d"%self.time_she

  def interpolation(self, q_array, n_pts=100):
    qmin = q_array[0]
    qmax = q_array[-1]
    q_array = flex.double( range(n_pts) )/float(n_pts)*(qmax-qmin) + qmin
    return q_array

  def iterate(self):
    if(self.Niter > 0): # need to build PDB object from the last PDB file
      iter_name=self.root+str(self.Niter)+".pdb"
      t1 = time.time()
      self.pdb = PDB(iter_name, method=self.method)
      self.nmode=self.pdb.Hessian(self.cutoff,self.nmode_init,self.scale_factor)-1
      self.time_nm += (time.time() - t1)

    self.n1=self.nmodes+1

    score=[]
    candidates=[]
    for kk in range(self.nmode_init-self.nmodes):
      self.modes=flex.int(range(self.nmodes))+7
      self.modes.append(kk+7+self.nmodes)
      self.starting_simplex=[]
      cand=flex.double(self.n1,0)
      for ii in range(self.n1):
        self.starting_simplex.append(flex.double(self.orth(ii,self.n1))*self.step_size+cand)
      self.starting_simplex.append(cand)

      self.optimizer = simplex.simplex_opt( dimension=self.n1,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  monitor_cycle = 4,
                                  tolerance=1e-1)
      self.x = self.optimizer.get_solution()
      candidates.append( self.x.deep_copy() )
      score.append( self.optimizer.get_score() )

    minscore=min(score[0:self.topn])
    print self.Niter,minscore, self.counter

    if((self.Niter % self.optNum) > 1):
      self.stopCheck(minscore)
    self.updateScore(minscore)
    minvec=candidates[score.index(minscore)]
    new_coord=flex.vec3_double(self.pdb.NMPerturb(self.modes,minvec))
    self.Niter=self.Niter+1
    iter_name=self.root+str(self.Niter)+".pdb"
    self.pdb.writePDB(new_coord,iter_name)
    if(self.Niter % self.optNum == 0):
      processed_pdb, pdb_inp = self.pdb_processor.process_pdb_files( pdb_file_names=[iter_name] )
      new_coord = geo_opt(processed_pdb, self.log)
      self.pdb.writePDB(new_coord,iter_name)
    if(not self.stop):
#    if(self.Niter < 50):
      self.iterate()

  def updateScore(self,minscore):
    if((self.Niter % self.optNum) == 0):
      self.minscore=minscore
    elif(minscore < self.minscore):
      self.minDev = self.minscore-minscore
      self.minscore=minscore

  def stopCheck(self, minscore):
    print self.minscore, minscore, self.minDev
    if(self.minscore < minscore):
      self.stop = True
#    if(minscore < (self.minDev*5.0)):
    if(self.minscore-minscore < (self.minDev*0.01)):
      self.stop = True

  def orth(self,indx,n):
    vec=[0]*n
    vec[indx]=1
    return vec

  def target(self, vector):
    self.counter += 1
    result = 0
    length=flex.sum(flex.pow(vector,2))
    if(length > self.Rmax2):
        result = 1e30
    else:
        new_coord = self.pdb.NMPerturb(self.modes,vector)
        t1 = time.time()
        self.she_engine.engine.update_coord(flex.vec3_double(new_coord),self.new_indx)
        new_I = self.she_engine.engine.I()
        self.time_she += (time.time() - t1)
        var = self.expt_s
        s,o = she.linear_fit(new_I, self.expt_I, var)
        result = flex.sum( flex.pow2( (self.expt_I-(s*new_I+o))  /self.expt_s ))
    return result

def geo_opt(processed_pdb ,log):
  sites_cart = geometry_minimization.run(
    processed_pdb_file = processed_pdb,
    log = log )
  return sites_cart



def run(args,log=sys.stdout):
  params = get_input(args, master_params, "pr_ref", banner, print_help )
  if (params is None): exit()

  t1=time.time()
  flex.set_random_seed(0)
  start_pdb=params.pr_ref.start_pdb
  target_i=params.pr_ref.target_I
  method=params.pr_ref.method
  n_modes=params.pr_ref.n_modes
  total_modes=n_modes+4
  max_rmsd=params.pr_ref.max_rmsd
  backbone_scale=params.pr_ref.backbone_force
  weight =params.pr_ref.weight
  prefix=params.pr_ref.prefix
  nmref(start_pdb,target_i,total_modes,n_modes,max_rmsd,backbone_scale,prefix, weight=weight,method=method)
  t2=time.time()
  print "\n start at: ", time.ctime(t1), "\n finished at: ", time.ctime(t2)

def print_help(out):
  print "\nUsage:\n sastbx.refine_i start_pdb=PDB_FILE target_I=target_I n_modes=Num_normal_mode max_rmsd=Max_Rmsd prefix=Prefix\n"
  print
  print "Required:"
  print "  start_pdb: the starting model in PDB format "
  print "  target_I : the experimental SAXS data       "
  print
  print "Optional:"
  print "  n_modes  :  Number of low-frequency normal modes for each perturbation, maximum 10, default=3"
  print "  max_rmsd :  maximum conformational change due to each perturbation "
  print "  prefix   :  prefix of the output "


if __name__ == "__main__":
  run(sys.argv[1:])
