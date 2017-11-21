import sys
from scitbx import simplex
import math
from scitbx.array_family import flex
import pr_tools
import random,time
from sastbx.data_reduction import saxs_read_write
import bounds as b

class rcs_fitter(object):
  def __init__(self, n_params, d_max, data, bounds, n_int=25,alpha=1.0):
    self.n = n_params
    self.d_max = d_max
    self.data = data
    self.n_int = n_int
    self.alpha = alpha
    self.gaussian = 0.0
    self.bounds=bounds
    self.step_size=2
    if(self.gaussian > 0):
      self.bounds.min = self.bounds.mean - self.gaussian * self.bounds.var
      self.bounds.max = self.bounds.mean + self.gaussian * self.bounds.var

    # make a pofr please
    self.pofr = pr_tools.pofr(self.d_max,self.n,n_int=self.n_int, m_int=self.n_int)
  
    #make a random simplex please
    self.starting_simplex = []
    cand = flex.random_double(self.n)
    for ii in range(self.n):
       self.starting_simplex.append(flex.double(self.orth(ii,self.n))*self.step_size + cand)
    self.starting_simplex.append(cand)

#      self.starting_simplex.append((flex.random_double(self.n)*2-1)*self.step_size)

    self.optimizer = simplex.simplex_opt( dimension=self.n,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-8)     
    self.solution = self.optimizer.GetResult()
    
  def orth(self,indx,n):
    vec=[0]*n
    vec[indx]=1
    return vec

  def target(self, vector):
    self.pofr.update( vector )
    if( self.outofbox(vector) ):
#      print "&"
      return 10e12
    calc_data = self.pofr.i( self.data.q )
    score = calc_data - self.data.i
    score = flex.pow(score,2.0)
    #score = flex.pow(score/self.data.s,2.0)
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
    self.n_trial = n_trial
    self.dmax=dmax
    bounds = b.bounds(n_param, "cheb_coef.dat")
    self.r = self.dmax*flex.double( range(101) )/100.0
    self.best_pr = flex.double()
    self.mean_pr = self.r * 0.0
    self.mean2 = self.mean_pr
    self.error = self.mean_pr

    for ii in range( n_trial ):
      tmp_object = rcs_fitter(n_param, dmax, data, bounds, alpha=alpha)
      self.trials.append( tmp_object )
    self.collect_scores()
    self.estimate_error()
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
    print '&'
    for rr, ff, ee in zip(self.r,self.best_pr,self.error):
      print rr, ff,ee
    print '&'
    for rr, ff, ee in zip(self.r,self.mean_pr,self.error):
      print rr, ff,ee

 
  def estimate_error(self):
    for ii in xrange(self.n_trial):     
      pr = self.trials[ii].get_best_pofr()
      pr.pofx.normalize()
      pofr = pr.f (self.r)
      self.mean_pr = self.mean_pr + pofr
      self.mean2 = self.mean2 + flex.pow(pofr,2.0)
      if (ii == self.chi_index):
         self.best_pr = pofr
    self.mean_pr /= self.n_trial
    self.mean2 /= self.n_trial
    self.error = flex.pow(self.mean2 - flex.pow(self.mean_pr,2.0),0.5)

  


 

  








def run(filename, dmax0):
  t1=time.time()
  flex.set_random_seed(0)
  data = saxs_read_write.read_standard_ascii_qis(filename) 
  data.multiply_add(1.0/data.i[0] , 0.0)
  dmax = dmax0
  for ii in xrange(1):
    fitters = random_start_fixed_dmax( data, dmax, 6, 10, 0 )
    fitters.collect_scores()


 



def get_input(args):
  filename=args[0]
  d=int(args[1])
#  t1 = time.time()
  run(filename,d)
#  t2 = time.time()
#  print t2-t1

if __name__ == "__main__":
  get_input(sys.argv[1:])
#  print "OK"
