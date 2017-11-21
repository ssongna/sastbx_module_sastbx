import sys,os
from scitbx.array_family import flex
from math import exp,log,sqrt
import random,time
from iotbx import pdb
from mmtbx.command_line import geometry_minimization
from mmtbx.utils import process_pdb_file_srv
from sastbx.intensity import she
from sastbx.md import beads
from sastbx.data_reduction import saxs_read_write
from sastbx.interface import get_input
import iotbx.phil


master_params = iotbx.phil.parse("""\
mc{
  start_pdb = None
  .type=path
  .help="the structure to be refined"

  target_I = None
  .type=path
  .help="the pr derived from scattering profile"

  max_rmsd = 0.5
  .type=float
  .help="maximum rmsd allowed for each iteration"

  backbone_force = 5.0
  .type=float
  .help="force constant of backbone (relative to the rest)"

  weight=i* s
  .type=choice
  .help="weighting scheme for the target function"

  prefix="ref_mc_"
  .type=str
  .help="prefix of the output pdb file"

}

""")

banner = """
============================================================================
       Structure Refinement
============================================================================
"""

class PDB(object):
  def __init__(self,file_name, method='ca', n_res_in_block=2, max_n_block=400, cutoff=12, backbone=10 ):
    self.method=method
    self.n_res_in_block=n_res_in_block
    self.max_n_block=max_n_block
    self.block_start=flex.int()
    self.label=[]
    self.xyz=flex.vec3_double()
    self.natm=0
    self.readPDB(file_name)
#    self.crd=flex.double()
#    for xyz in self.xyz:
#      for value in xyz:
#        self.crd.append(value)
    self.beads = beads( self.xyz, self.block_start, cutoff, backbone)

  def readPDB(self, file_name):
    self.pdbi = pdb.hierarchy.input(file_name=file_name)
    if(len( self.pdbi.hierarchy.models() ) == 0):
      return None

    self.atoms = self.pdbi.hierarchy.models()[0].atoms()

    if(self.method == 'ca'):
      # keep track of the atom types we have encountered
      for atom in self.atoms:
        #if(not atom.hetero):
          self.xyz.append( atom.xyz )
          if(atom.name in [" CA "," P  "," C5 "]):
            self.block_start.append(self.natm)
          self.natm=self.natm+1
      self.n_block = self.block_start.size()
      print "CA-based Elastic Network Model will be used"
      print "number of points used for ENM: ", self.block_start.size()
    else:
      self.residue_groups=self.pdbi.hierarchy.models()[0].residue_groups()
      n_atm_in_res=flex.int()
      res_start=flex.int()
      self.xyz=self.atoms.extract_xyz()
      self.natm = self.xyz.size()
      self.atoms.reset_i_seq()
      for rg in self.residue_groups:
        rg_atoms = rg.atoms()
        n_atm_in_res.append( rg_atoms.size() )
        res_start.append( rg_atoms[0].i_seq )
      n_res=res_start.size()
      self.n_res_in_block=max(n_res/self.max_n_block, self.n_res_in_block)
      ### Merge the residues to block ###
      self.n_atm_in_block=flex.int()
      for ii in range(0,n_res,self.n_res_in_block):
      #  merge_n_atm = flex.sum( n_atm_in_res[ii:ii+self.n_res_in_block] )
      #  self.n_atm_in_block.append(merge_n_atm)
        self.block_start.append( res_start[ii] )
      self.n_block=self.block_start.size()
      ### is there any remaining residues ###
      nres_left = n_res - self.n_block*self.n_res_in_block
      if(nres_left>0):
      #  merge_n_atm = flex.sum(n_atm_in_res[-res_left:])
      #  self.n_atm_in_block.append(merge_n_atm)
        self.block_start.append(res_start[-res_left])
        self.n_block +=1
      print "RTB-based Elastic Network Model will be used"
      print "number of blocks used for ENM: ", self.block_start.size()

  def writePDB(self, crd, file_name):
    for atom,xyz in zip(self.atoms, crd):
      atom.set_xyz(new_xyz=xyz)
    self.pdbi.hierarchy.write_pdb_file( file_name=file_name, open_append=False)


  def perturb(self, dx):
    all_dxyz = self.beads.project(dx)
    new_xyz = self.xyz + all_dxyz
    return new_xyz




class mc_refine(object):
  def __init__(self,start_pdb, target_I, max_rmsd,backbone_scale,prefix, nstep_per_cycle=100, method='ca', weight='i', log='tmp.log'):
    self.counter = 0
    self.topn=3
    self.Niter=0
    self.method=method
    self.cutoff=12
    self.log = open(log, 'w')
    self.nstep_per_cycle = nstep_per_cycle
    self.pdb_obj=PDB(start_pdb, method=self.method)
    crystal_symmetry = self.pdb_obj.pdbi.xray_structure_simple().\
        cubic_unit_cell_around_centered_scatterers(
        buffer_size = 10).crystal_symmetry()
    self.pdb_processor = process_pdb_file_srv( crystal_symmetry=crystal_symmetry )

    self.expt = saxs_read_write.read_standard_ascii_qis(target_I)
    self.q = self.expt.q
    self.expt_I = self.expt.i
    self.expt_s = self.expt.s
    if(self.q.size() > 20):
      self.q = self.interpolation(self.q, n_pts=20)
      self.expt_I = flex.linear_interpolation( self.expt.q, self.expt.i, self.q)
      self.expt_s = flex.linear_interpolation( self.expt.q, self.expt.s, self.q)

    if( weight=='i'):
      self.expt_s = flex.sqrt( self.expt_I )

    self.time_nm =0
    self.time_she=0
    self.she_engine=she.she(start_pdb,self.q)
    self.natom = self.pdb_obj.natm
    self.nbeads= self.pdb_obj.n_block
    self.scale_factor=backbone_scale
    time1=time.time()
    self.time_nm += (time.time() - time1 )

    self.root=prefix
    self.drmsd = max_rmsd
    self.step_size = self.drmsd*3
    self.threshold = self.drmsd**2.0

    self.new_indx=flex.int(range(self.natom))
    self.stop = False
    self.minscore=1e20
    self.minDev = 0  #minimum deviations of refined structures, compared to refined structure from the previous step
    self.optNum = 10 #number of iterations between geometry optimization

    #self.estimate_init_weight()
    #self.restraint_weight *= 8  ## contribute 8x of chi initially
    self.iterate()
    self.log.close()

    print "time used for NM : %d"%self.time_nm
    print "time used for she: %d"%self.time_she

  def interpolation(self, q_array, n_pts=100, qmax=0.25):
    qmin = q_array[0]
    qmax = min(qmax, q_array[-1])
    q_array = flex.double( range(n_pts) )/float(n_pts)*(qmax-qmin) + qmin
    return q_array

  def build_restraint(self, n_restraint):
    r_a = flex.random_size_t(n_restraint)%self.nbeads
    r_b = flex.random_size_t(n_restraint)%self.nbeads
    bead_pairs = flex.tiny_size_t_2()
    for a,b in zip(r_a, r_b):
      bead_pairs.append( (a,b) )
    return bead_pairs

  def estimate_init_weight(self):
    ## initial chi-score ##
    new_I = self.she_engine.engine.I()
    var = self.expt_s
    s,o = she.linear_fit(new_I, self.expt_I, var)
    chi_score = flex.sum( flex.pow2( (self.expt_I-(s*new_I+o))  /self.expt_s ))

    n_restraint = self.nbeads**2/20
    tot_res = 0
    for ii in range(10): ## perturb 10 times to estimate initial restraint ##
      vector = flex.random_double(self.nbeads*3)*self.step_size
      self.restraints = self.build_restraint(n_restraint)
      self.new_xyz = self.pdb_obj.perturb(vector)
      restraint=self.pdb_obj.beads.restraint(self.restraints, self.new_xyz)
      tot_res += restraint
    mean_res = tot_res/10.0
    self.restraint_weight = chi_score/mean_res


  def iterate(self):
    n_restraint = self.nbeads**2/20 #/(self.Niter+1)
    self.restraints = self.build_restraint(n_restraint)
    best_xyz = self.pdb_obj.xyz
    lowest_score = self.target( flex.vec3_double(self.nbeads) )

    current_score=lowest_score
    for kk in range(self.nstep_per_cycle):
      dx = flex.random_double(self.nbeads*3)*self.step_size
      dx = flex.vec3_double(dx)
      s = self.target(dx)
      if( s<lowest_score ):
        lowest_score=s
        best_xyz = self.new_xyz.deep_copy()
        self.pdb_obj.beads.update( best_xyz )
        current_score = s
      else:
        nothing=0 ## to be filled

    self.Niter=self.Niter+1
    iter_name=self.root+str(self.Niter)+".pdb"
    self.pdb_obj.writePDB(best_xyz,iter_name)
    self.pdb_obj.xyz = best_xyz.deep_copy()
    if(self.Niter % self.optNum == 0):
      processed_pdb, pdb_inp = self.pdb_processor.process_pdb_files( pdb_file_names=[iter_name] )
      new_coord = geo_opt(processed_pdb, self.log)
      self.pdb_obj.writePDB(new_coord,iter_name)
#    if(not self.stop):
    if(self.Niter<20):
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


  def target(self, dxyz):
    self.move_on=False
    chi_score=1e30
    if(True):
      for ii in range(50):
        dxyz = self.pdb_obj.beads.relax( self.restraints, dxyz )
        print ii, ' ',
        if ( self.pdb_obj.beads.restraint( self.restraints, dxyz )< self.threshold ):
          self.new_xyz = self.pdb_obj.perturb(dxyz)
          self.move_on=True
          break
    if(self.move_on):
      t1 = time.time()
      self.she_engine.engine.update_coord(self.new_xyz,self.new_indx)
      new_I = self.she_engine.engine.I()
      self.time_she += (time.time() - t1)
      var = self.expt_s
      s,o = she.linear_fit(new_I, self.expt_I, var)
      chi_score = flex.sum( flex.pow2( (self.expt_I-(s*new_I+o)) )  /self.expt_s )
      #restraint=self.pdb_obj.beads.restraint(self.restraints, self.new_xyz)
    #tot = chi_score + restraint*self.restraint_weight
    print chi_score
    self.counter += 1
    return chi_score

def geo_opt(processed_pdb ,log):
  sites_cart = geometry_minimization.run(
    processed_pdb_file = processed_pdb,
    log = log )
  return sites_cart



def run(args,log=sys.stdout):
  params = get_input(args, master_params, "mc", banner, print_help )
  if (params is None): exit()

  t1=time.time()
  flex.set_random_seed(0)
  start_pdb=params.mc.start_pdb
  target_i=params.mc.target_I
  max_rmsd=params.mc.max_rmsd
  backbone_scale=params.mc.backbone_force
  weight =params.mc.weight
  prefix=params.mc.prefix
  mcref=mc_refine(start_pdb,target_i,max_rmsd,backbone_scale,prefix, weight=weight)
  t2=time.time()
  print "\n start at: ", time.ctime(t1), "\n finished at: ", time.ctime(t2)

def print_help(out):
  print "\nUsage:\n sastbx.mc start_pdb=PDB_FILE target_I=target_I max_rmsd=Max_Rmsd prefix=Prefix\n"
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
