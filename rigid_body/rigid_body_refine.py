import os, sys
from scitbx.array_family import flex
import math
import random
from scitbx import simplex
from scitbx import differential_evolution as de
from sastbx.data_reduction import saxs_read_write
from sastbx.rigid_body import rigidbody as rb
from sastbx.rigid_body import rb_engine as rbe

class refine_rb(object):
  def __init__(self,rbe_obj,expt_data, simplex_trial=10, n_step=20):
    self.data=expt_data
    self.pr = expt_data.i
    self.rbe=rbe_obj
    self.nbody=rbe_obj.nbody
    self.n = 3

    self.count=0 
    self.simplex_trial = simplex_trial
    self.n_step = n_step

#    tmp_v = flex.double(6,0)
#    self.target(tmp_v)   
    self.x=None
    self.x = flex.double([0,0,0])
    self.vector=flex.double([5,5,5])*0
    self.angle=flex.double(3,0)
   
#    self.de()
    for ii in range(self.n_step):
      if(ii % 2 == 0):
        self.translate=True
        self.x = self.vector.deep_copy()
      else:
        self.translate = False
	self.x = self.angle.deep_copy()

      self.simplex()


  def de(self):
    self.domain = [ (-0.5,0.5) ]*self.n
    self.optimizer = de.differential_evolution_optimizer( self,
                                                        population_size=self.n*2,
                                                        n_cross=1,
                                                        f=0.85,
                                                        eps=1e-5,
                                                        monitor_cycle=50,
                                                        show_progress=False)
    print '#', list(self.x)

  def simplex(self):
    self.simplex_scores = []
    self.simplex_solutions = []
    for ii in xrange(self.simplex_trial):
      #make a random simplex
      self.starting_simplex = []
      for ii in range(self.n):
        self.starting_simplex.append(random.random()*(flex.random_double(3)*2-1.0)*2+self.x)
      self.starting_simplex.append(self.x)

      self.optimizer = simplex.simplex_opt( dimension=self.n,
                                    matrix  = self.starting_simplex,
                                    evaluator = self,
				    max_iter=50,
                                    tolerance=1e-4)

      self.solution = self.optimizer.GetResult()
      self.score = self.target( self.solution )
      self.simplex_scores.append( self.score )
      self.simplex_solutions.append( self.solution )

    best_index = flex.min_index( flex.double(self.simplex_scores) )
    self.score = self.simplex_scores[ best_index ]
    self.solution = self.simplex_solutions[ best_index ]
    if(self.translate):
      self.vector=self.solution.deep_copy()
      print "translate:", list(self.vector)
    else:
      self.angle=self.solution.deep_copy()
      print "rotate:", list(self.angle)
    self.target(self.solution)
     

  def target(self, vector):
#    self.rbe.rotate_translate( vector[0:3], vector[3:], 1)
    if(self.translate):
      self.rbe.rotate_translate( vector, self.angle, 1)
    else:
      self.rbe.rotate_translate( self.vector, vector, 1)

    calc_pr = self.rbe.get_norm_pr()
    print '&',list(vector)
#    for cp,ep in zip(calc_pr, self.pr):
#      print cp,ep
    score = flex.sum( flex.abs(calc_pr-self.pr)*flex.exp(-self.pr) )
    print '#',score
    return score


#=================================================================================================#
#
#=================================================================================================#

def test(args):
  rbs=[]
  file=args[0]
  expt_data=saxs_read_write.read_standard_ascii_qis( file )
  dmax=expt_data.q[-1]+1
  pdb=None
  for arg in args[1:]:
    pdb= rbe.PDB(arg)
    rbs.append(rb(pdb.xyz, pdb.CA_indx, dmax))

  rb_eng = rbe.rb_engine(rbs,int(dmax) )

  refine = refine_rb( rb_eng, expt_data)
  
  pdb.writePDB(refine.rbe.rbs[1].get_crd(),'refined.pdb')




#==============================================================================================#
#
#==============================================================================================#
if __name__ == "__main__":
  test(sys.argv[1:])

