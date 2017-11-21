from scitbx.array_family import flex
from stdlib import math as smath

class two_stage_basic(object):
  def __init__(self, q_array, int_array, sig_array, s1=0.02, s2=1.0):
    self.q = q_array.deep_copy()
    self.i = int_array.deep_copy()
    self.s =  sig_array.deep_copy()
    self.s1 = s1
    self.s2 = s2

  def error_model(self,rmax):
    q1 = self.q1(rmax)
    q2 = self.q2(rmax)

    sel1 = flex.bool(self.q <= q1 )
    sel2 = flex.bool(self.q > q1 )
    sel3 = flex.bool(self.q > q2 )

    rc = (self.s2-self.s1)/(q2-q1)
    tmp = self.q*rc-q1*rc+self.s1
    tmp.set_selected(sel1,self.s1)
    tmp.set_selected(sel3,self.s2)    

    tmp = tmp*tmp+self.s*self.s
    return tmp

  def error_model_smear(self,rmax,delta=5):
    vv = self.q*0.0
    icount = 0.0
    for ii in range(-delta,delta):
      tmp_rmax = rmax+ii
      vv += self.error_model(tmp_rmax)
      icount+=1
    vv = vv / icount
    return vv

  def q1(self,rmax):
    a=-3.724
    b=1.269
    c=-0.2639
    x = float(rmax)
    if x<10:
      x = 10.0
    if x > 200:
      x = 200.0
    x = smath.log(x)
    y = a+b*x+c*x*x
    return smath.exp(y)

  def q2(self,rmax):
    a = 3.019
    b = -1.28945
    x = float(rmax)
    if x<10:
      x = 10.0
    if x > 200:
      x = 200.0
    x = smath.log(x)
    y = a+b*x
    return smath.exp(y)

  
if __name__ == "__main__":
  q = flex.double(range(100) )/400.0
  i = q*0.0+1
  s = q*0.0+0.0001
  em = two_stage_basic(q,i,s,s2=0.5)
  ns = em.error_model_smear(15.0,4)
  for qq,ss in zip(q,ns):
    print qq,smath.sqrt(ss) 


