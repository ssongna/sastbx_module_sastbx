import sys
from scitbx import simplex
from math import exp,log,sqrt
from integrating_a_weighted_sinc_function import *
import random,time
from libtbx.test_utils import run_command

def pr2I(ni,kernel,q_array,k1):
  I_q = flex.double()
  k1 = 0.424
  for q in q_array:
    iq = kernel.integrate(q,ni)
    iq = iq*exp(-k1**2.0*q**2.0)
    I_q.append(iq)
  return I_q



class fit(object):
  def __init__(self,Nb,dMax,I_file,npr_for_gl_int):
    self.Nb = Nb
    self.dMax = dMax
    self.N_pr = npr_for_gl_int
    self.topn=3
    self.step_size=5
    self.weighted = False
    self.scale = 0
    self.minscore=100000
    self.minDev=100000
    self.stop = False
    self.Niter = 0

##### for intensity ##
    self.q,self.expt_I = self.readExpt(I_file)
    self.b_s=flex.random_double(self.Nb)*0.0
    self.fpr = function(self.Nb,self.N_pr,1.0,self.dMax)
    self.fpr.load_coefs(self.b_s[0:])
    k1 = self.b_s[0]
    new_I = pr2I(self.N_pr,self.fpr,self.q,k1)
    self.candidate = self.b_s.deep_copy()
    self.calc_I = new_I.deep_copy()
    self.optimize()

  def optimize(self):
    self.candidate = flex.random_double(self.Nb)
    candidates=[]
    score=[]
    for kk in range(self.topn*10):
      if(kk == 0):
        ab_s=self.candidate
      else:
        ab_s=(flex.random_double(self.Nb)*self.step_size+self.candidate)
      result=self.target(ab_s)
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

      self.optimizer = simplex.simplex_opt( dimension=self.Nb,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-10)

      self.x = self.optimizer.GetResult()
      candidates[kk]=self.x.deep_copy()
      score[kk]=self.target(self.x)

    minscore=min(score[0:self.topn])
    minvec=candidates[score.index(minscore)]
    self.fpr.load_coefs(minvec[0:])
    k1 = minvec[0]
    new_I = pr2I(self.N_pr,self.fpr,self.q,k1)
    print '&',list(minvec), "coef", minscore
    for r,n in zip(self.q,new_I):#,self.calc_I):
      print r,n

    self.candidate=minvec.deep_copy()
    self.calc_I = new_I.deep_copy()
    self.Niter=self.Niter+1
    print '&'
    for r,n in zip(self.q,self.expt_I):#,self.calc_I):
      print r,n

    r = flex.double(range(-50,50))/50.0+1e-7
    pr=self.fpr.get_p_of_r(r)
    sum_pr = flex.sum(pr)*self.dMax
    pr=pr/sum_pr*100
    r=(r+1)*self.dMax/2.0
    print '&'
    for ri,pri in zip(r,pr):
      print ri,pri
    self.prior=pr.deep_copy()
    self.Nb=self.Nb+self.Nb
    self.fpr=function(self.Nb,self.N_pr,1.0,self.dMax)
    if(self.Niter < 2): 
      self.optimize()


  def readExpt(self,filename):
    q = flex.double()
    Iq = flex.double()
    file = open(filename,'r')
    total=0.0
    for line in file:
      keys = line.split()
      if(keys[0] != '&'):
        if(float(keys[0]) < 0.5):
          q.append(float(keys[0]))
          Iq.append(float(keys[1]))
    I0 = max(Iq)
    return q,Iq/I0

  def orth(self,indx,n):
    vec=[0]*n
    vec[indx]=1
    return vec 

  def target(self, vec):
    k1 = vec[0]
    self.fpr.load_coefs(vec[0:])
    new_I = pr2I(self.N_pr,self.fpr,self.q,k1)
    result = flex.sum(flex.pow((self.expt_I-new_I),2.0))
    if(self.Niter > 0):
      r = flex.double(range(-50,50))/50.0+1e-7
      pr=self.fpr.get_p_of_r(r)
      sum_pr = flex.sum(pr)*self.dMax
      pr=pr/sum_pr*100
      entropy=flex.sum(flex.pow((self.prior-pr),2.0))
      #entropy=flex.sum((self.prior-pr)*flex.log((self.prior/pr)))
      print '&', entropy, result
      result += entropy*0.2
    return result

def run(filename, dmax):
  t1=time.time()
  flex.set_random_seed(0)
  fit(5,dmax,filename,20)
  t2=time.time()
#  print "\n start at: ", time.ctime(t1), "\n finished at: ", time.ctime(t2)

def test(args):
  filename=args[0]
  d=int(args[1])
  run(filename,d)

if __name__ == "__main__":
  test(sys.argv[1:])
#  print "OK"
