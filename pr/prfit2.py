import sys
from simplex import *
from math import exp,log,sqrt
from integrating_a_weighted_sinc_function import *
import random,time
from libtbx.test_utils import run_command

def pr2I(ni,kernel,q_array):
  I_q = flex.double()
  for q in q_array:
    iq = kernel.integrate(q,ni)
    I_q.append(iq)
  return I_q

class fit(object):
  def __init__(self,Nb,dMax,I_file,npr_for_gl_int):
    self.Nb = Nb
    self.dMax = dMax
    self.N_pr = npr_for_gl_int
    self.topn=1
    self.step_size=1.0
    self.weighted = False
    self.scale = 1
    self.minscore=100000
    self.minDev=100000
    self.stop = False

##### for intensity ##
    self.q,self.expt_I = self.readExpt(I_file)
    self.a_s=flex.double([1.0,1.0])
    self.b_s=flex.random_double(self.Nb)
    self.fpr = function(self.Nb,self.N_pr,1.0,self.dMax)
    self.r = flex.double()
    r,pr = self.readExpt("mbp.pr")
    rmax = max(r)
    self.r = r/rmax *2.0 - 1.0

    self.optimize()

  def optimize(self):
    candidates=[]
    score=[]
    for kk in range(self.topn):
      ab_s=(flex.random_double(self.Nb))
      result=self.target(self.expt_I,ab_s)
      insert = 0
      for ii in range(len(score)):
        if(score[ii] > result):
          score.insert(ii,result)
          candidates.insert(ii,ab_s)
	  insert=1
          break
      if(insert==0):
	score.append(result)
	candidates.append(ab_s)

    for kk in range(self.topn):
      self.starting_simplex=[]
      cand=candidates[kk]
      for ii in range(self.Nb):
        self.starting_simplex.append(flex.double(self.orth(ii,self.Nb))*self.step_size+cand)
      self.starting_simplex.append(cand)

      self.optimizer = simplex_opt( dimension=self.Nb,
                                  matrix  = self.starting_simplex,
                                  expt_value = self.expt_I,
                                  evaluator = self,
                                  tolerance=1e-8)

      self.x = self.optimizer.GetResult()
      candidates[kk]=self.x.deep_copy()
      score[kk]=self.target(self.expt_I,self.x)

    minscore=min(score[0:self.topn])
    minvec=candidates[score.index(minscore)]
    self.fpr.load_coefs(minvec)
    new_I = self.fpr.get_p_of_r(self.r)
    print '&',list(minvec)
    for r,o,n in zip(self.q,self.expt_I,new_I):
      print r,o,n
    self.q = flex.double(range(2,50))/100.0
    new_I = pr2I(self.N_pr,self.fpr,self.q)
    print '&'
    for q,i in zip(self.q,new_I):
      print q,i*1000000 

#    self.Niter=self.Niter+1
#    if(not self.stop): 
#    if(self.Niter < 100): 
#      self.optimize()

  def readExpt(self,filename):
    q = flex.double()
    Iq = flex.double()
    file = open(filename,'r')
    total=0.0
    for line in file:
      keys = line.split()
      q.append(float(keys[0]))
      Iq.append(float(keys[1]))
    return q,Iq

  def orth(self,indx,n):
    vec=[0]*n
    vec[indx]=1
    return vec 

  def target(self, expt, vec):
    result = 0
    self.fpr.load_coefs(vec)
    new_I = self.fpr.get_p_of_r(self.r)
    #print list(new_I), "I"
    if(self.scale==0):
      self.scale=float(flex.mean(flex.abs(self.expt_I-new_I)))/float(flex.mean(self.expt_I))
    for old,new in zip(self.expt_I,new_I):
      if(old != new):
	if(self.weighted):
          result += exp(-(self.scale*old/(old-new))**2.0)*(old-new)**2.0
	else:
          result += (old-new)**2.0
    result=result*100
  #  print list(vec), result, "R"
    return result

def calc_pr(function_pr,r):
  pr = flex.double()
  for ri in r:
    pri = function_pr.get_p_of_r(ri)
  pr.append(pri)
  return pr


def run(filename):
  t1=time.time()
  flex.set_random_seed(0)
  fit(6,66,filename,10)
  t2=time.time()
#  print "\n start at: ", time.ctime(t1), "\n finished at: ", time.ctime(t2)

if __name__ == "__main__":
  run(sys.argv[1])
#  print "OK"
