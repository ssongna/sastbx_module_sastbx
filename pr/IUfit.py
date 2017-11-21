import sys
from scitbx import simplex
import math
from scitbx.array_family import flex
import pr_tools
import random,time
from sastbx.data_reduction import saxs_read_write
import bounds as b

class rcs_fitter(object):
  def __init__(self, n_params, d_max, data, n_int=35,alpha=1.0):
    self.n = n_params
    self.d_max = d_max
    self.data = data
    self.n_int = n_int
    self.alpha = alpha
    self.gaussian = 2.0
    self.bounds=b.bounds(self.n,"cheb_coef.dat")
    if(self.gaussian > 0):
      self.bounds.min = self.bounds.mean - self.gaussian * self.bounds.var
      self.bounds.max = self.bounds.mean + self.gaussian * self.bounds.var

    # make a pofr please
    self.pofr = pr_tools.pofr(self.d_max,self.n,n_int=self.n_int, m_int=self.n_int)
  
    #make a random simplex please
    self.starting_simplex = []
    for ii in range(self.n+1):
      self.starting_simplex.append(flex.random_double(self.n)*2-1)

    self.optimizer = simplex.simplex_opt( dimension=self.n,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-3)     
    self.solution = self.optimizer.GetResult()
    


  def target(self, vector):
    self.pofr.update( vector )
    calc_data = self.pofr.i( self.data.q )
    score = calc_data - self.data.i
    score = flex.abs(score)/self.data.s
    score = flex.sum(score)/calc_data.size()
    t,a,b = self.pofr.entropy_simple()
    total_score = score+t*self.alpha
    #print "#",list(vector), total_score
    return total_score

  def outofbox(self,vector):
    left = vector > self.bounds.min
    right = vector < self.bounds.max
    for ll in left:
      if( not ll):
        return True
    for rr in right:
      if(not rr):
        return True
    return False



  def get_best_pofr(self):
    self.pofr.update( self.solution )
    return self.pofr

  def get_scores(self):
    self.pofr.update( self.solution )
    calc_data = self.pofr.i( self.data.q )
    score = calc_data - self.data.i
    score = flex.abs(score)/self.data.s
    score = flex.sum(score)/calc_data.size()
    t,a,b = self.pofr.entropy_simple()
    return score, t, a, b
    

    
    
    


class random_start_fixed_dmax(object):
  def __init__(self, data, dmax, n_param, n_trial=10, alpha=1):
    self.trials = []
    self.dmax=dmax
    for ii in range( n_trial ):
      tmp_object = rcs_fitter(n_param, dmax, data, alpha=alpha)
      self.trials.append( tmp_object )
    self.collect_scores()
    self.show_best_pr()


  def collect_scores(self):
    self.chi=[]
    self.ent=[]
    for oo in self.trials:
      c,t,a,b= oo.get_scores()
      self.chi.append(c)
      self.ent.append(t)
      print "#", self.dmax, c,t,a,b
    self.chi_index = flex.min_index( flex.double(self.chi) )
    self.ent_index = flex.min_index( flex.double(self.ent) )



  def show_best_pr(self):
    r = self.dmax*flex.double( range(101) )/100.0
    f = self.trials[ self.chi_index ]
    f = f.get_best_pofr().f( r )
    for rr, ff in zip(r,f):
      print rr, ff

 


    

  


 

  








def run(filename, dmax0):
  t1=time.time()
  flex.set_random_seed(0)
  data = saxs_read_write.read_standard_ascii_qis(filename) 
  data.multiply_add(1.0/data.i[0] , 0.0)
  # fitttit = rcs_fitter( 6, dmax, data, alpha=1e1)
  for ii in xrange(1):
    dmax=dmax0+ii*10
    fitters = random_start_fixed_dmax( data, dmax, 6, 15, 0.2 )
    fitters.collect_scores()


 



def test(args):
  filename=args[0]
  d=int(args[1])
  run(filename,d)

if __name__ == "__main__":
  test(sys.argv[1:])
#  print "OK"
