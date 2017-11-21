import sys,os
from scitbx import simplex
from scitbx.array_family import flex
from sastbx import refine
from sastbx.refine.PDB import PDB
from math import exp,log,sqrt
import random,time
from libtbx.test_utils import run_command
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
import libtbx.phil
import libtbx.phil.command_line
from cStringIO import StringIO

master_params = iotbx.phil.parse("""\
pr_ref{
  start_pdb = None
  .type=path
  .help="the structure to be refined"

  target_pr = None
  .type=path
  .help="the pr derived from scattering profile"

  n_modes = 3
  .type=int
  .help="number of normal modes used to refine the structure"

  max_rmsd = 0.2
  .type=float
  .help="maximum rmsd allowed for each iteration"

  backbone_force = 10.0
  .type=float
  .help="force constant of backbone (relative to the rest)"

  prefix="pr_ref"
  .type=str
  .help="prefix of the output pdb file"

}

""")



class nmref(object):
  def __init__(self,start_pdb, target_pr, ntotal,nmodes,max_rmsd,backbone_scale,prefix):
    self.counter = 0
    self.nmode_init = ntotal
    self.nmodes = nmodes
    self.topn=10
    self.Niter=0
    self.modes=flex.int(range(self.nmode_init))+7
    self.cutoff=10
    self.weighted = True

##### for histogram ##
    r,self.expt = self.readPr(target_pr)
    self.expt2 = flex.pow(self.expt,2.0)
    #print list(self.expt)
    start_name=start_pdb
    self.pdb = PDB(start_name)
    self.natom = self.pdb.natm
    self.scale_factor=backbone_scale
    self.pdb.Hessian=self.pdb.Hessian(self.cutoff,self.nmode_init,self.scale_factor)

    self.root=prefix
    self.dMax=max(r)
    self.n_slot=r.size()
    self.scale=0
    self.drmsd = max_rmsd
    self.Rmax2 = self.natom * (self.drmsd) **2.0
    self.step_size= sqrt(self.Rmax2 /self.nmodes)*2.0

    self.new_indx=flex.int(range(self.natom))
    self.stop = False
    self.minscore=1e20
    self.minDev = 0  #minimum deviations of refined structures, compared to refined structure from the previous step
    self.optNum = 20 #number of iterations between geometry optimization
    self.iterate()

  def iterate(self):
    self.n1=self.nmodes
    if(self.Niter > 0): # need to build PDB object from the last PDB file
      iter_name=self.root+str(self.Niter)+".pdb"
      self.pdb = PDB(iter_name)
      self.pdb.Hessian=self.pdb.Hessian(self.cutoff,self.nmode_init,self.scale_factor)

# Generate random normal modes
    self.modes=flex.int(range(self.nmodes-1))+7
    self.modes.append(int(random.random()*(self.nmode_init-self.nmodes))+7+self.nmodes)
    self.scale=0

    candidates=[]
    score=[]
    for kk in range(self.topn*10):
      if(kk == 0):
        vec=flex.random_double(self.nmodes)*0
      else:
        vec=(flex.random_double(self.nmodes)-0.5)*2*self.step_size
      result=self.target(vec)
      insert = 0
      for ii in range(len(score)):
        if(score[ii] > result):
          score.insert(ii,result)
          candidates.insert(ii,vec)
          insert=1
          break
      if(insert==0):
        score.append(result)
        candidates.append(vec)

    for kk in range(self.topn):
      self.starting_simplex=[]
      cand=candidates[kk]
      for ii in range(self.n1):
        self.starting_simplex.append(flex.double(self.orth(ii,self.n1))*self.step_size+cand)
      self.starting_simplex.append(cand)

      self.optimizer = simplex.simplex_opt( dimension=self.n1,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-8)
      self.x = self.optimizer.get_solution()
      candidates[kk]=self.x.deep_copy()
      score[kk]=self.optimizer.get_score()

    minscore=min(score[0:self.topn])

    if((self.Niter % self.optNum) > 1):
      self.stopCheck(minscore)
    self.updateScore(minscore)
    minvec=candidates[score.index(minscore)]
    new_coord=flex.vec3_double(self.pdb.NMPerturb(self.modes,minvec))
    self.Niter=self.Niter+1
    iter_name=self.root+str(self.Niter)+".pdb"
    self.pdb.writePDB(new_coord,iter_name)
    if(self.Niter % self.optNum == 0):
       geo_opt(iter_name,iter_name)
    if(not self.stop):
#    if(self.Niter < 100):
      self.iterate()

  def readPr(self,filename):
    r = flex.double()
    pr = flex.double()
    file = open(filename,'r')
    total=0.0
    for line in file:
      keys = line.split()
      r.append(float(keys[0]))
      pr.append(float(keys[1]))
      total+=float(keys[1])
    return r,pr/total

  def updateScore(self,minscore):
    if((self.Niter % self.optNum) == 0):
      self.minscore=minscore
    elif(minscore < self.minscore):
      self.minDev += (self.minscore-minscore)
      self.minscore=minscore

  def stopCheck(self, minscore):
    print self.minscore, minscore, self.minDev, self.counter
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
        result = 100000000000000
    else:
        new_coord = self.pdb.NMPerturb(self.modes,vector)
        self.pdb.model.updateDistArray(flex.vec3_double(new_coord))
        dist_array=self.pdb.model.getDistArray()
        new_pr=self.pdb.model.Histogram(dist_array, self.dMax,self.n_slot)
        if(self.weighted):
          if(self.scale==0):
             self.scale=float(flex.mean(flex.abs(self.expt-new_pr)))/float(flex.mean(self.expt))
             self.scale2 = self.scale*self.scale
          result = flex.pow( (self.expt-new_pr), 2.0)
          result = flex.sum( flex.exp(-self.scale2*self.expt2/(result+10-12)) * result )
        else:
          result = flex.sum( flex.pow( (self.expt-new_pr), 2.0) )
        result=result*100
    return result

def geo_opt(input,output):
  log = "tmp.log"
  cmd = 'phenix.pdbtools %s --geometry-regularization modify.output.file_name=%s > %s' % (input,output, log)
  run_command(command=cmd)

def run(params,log):
  t1=time.time()
  flex.set_random_seed(0)
  start_pdb=params.pr_ref.start_pdb
  target_pr=params.pr_ref.target_pr
  n_modes=params.pr_ref.n_modes
  total_modes=n_modes+4
  max_rmsd=params.pr_ref.max_rmsd
  backbone_scale=params.pr_ref.backbone_force
  prefix=params.pr_ref.prefix
  nmref(start_pdb,target_pr,total_modes,n_modes,max_rmsd,backbone_scale,prefix)
  t2=time.time()
  print "\n start at: ", time.ctime(t1), "\n finished at: ", time.ctime(t2)

def get_input(args):
  if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil=master_params,
      home_scope="pr_ref")

    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log,expert_level=1)
    print >> log, "#phil __END__"
    run( params, log )

def print_help():
  print "\nUsage:\n sastbx.refine_pr start_pdb=PDB_FILE target_pr=target_pr n_modes=Num_normal_mode max_rmsd=Max_Rmsd prefix=Prefix\n"


if __name__ == "__main__":
  get_input(sys.argv[1:])
