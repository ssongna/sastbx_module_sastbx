from sastbx import rigid_body as rb
from scitbx.array_family import flex
import math
import time

def tst_sphere():
  wigner = rb.wigner3j_fast(30)

  value=math.sqrt(1.0/5.0)
  w= rb.wigner3j(1,1,2,1,1,-2)
  value = wigner.compute(1,1,2,1,1,-2)
  assert abs(w.get_value()-value)<1e-6

  value=math.sqrt(1.0/14.0)
  w=rb.wigner3j(2,2,3,1,2,-3)
  value = wigner.compute(2,2,3,1,2,-3)
  assert abs(w.get_value()-value)<1e-6
  
  value=math.sqrt(15.0/52003.0)*(-5.0)
  w=rb.wigner3j(5, 8, 10, 3, -3, 0)
  value = wigner.compute(5, 8, 10, 3, -3, 0)
  assert abs(w.get_value()-value)<1e-6

  w=rb.wigner3j(5,8,10,2,4,-6)
  value=-(32.0/3.0)*math.sqrt(26.0/2860165.0) 
  value = wigner.compute(5,8,10,2,4,-6)
  assert abs(w.get_value()-value)<1e-6

  value = wigner.compute(8,5,10,4,2,-6)
  print value
  value = wigner.compute(8,5,10,4,-2,-6)
  print value
  value = wigner.compute(5,8,10,2,4,-6)
  print value

  
  
def tst2():
  w= rb.wigner3j_zero(1,1,2)

  w=rb.wigner3j_zero(2,2,3)
  
  w=rb.wigner3j_zero(5, 8, 10)

  w=rb.wigner3j_zero(5,8,10)
  
  print rb.wigner3j(9,14,9,0,0,0).get_value()
  print rb.wigner3j_zero(9,14,9).get_value()
  print rb.wigner3j(8,14,9,0,0,0).get_value()
  print rb.wigner3j(8,14,9,1,1,-2).get_value()

  
  


if __name__=="__main__":
  start=time.time()
  for ii in range(1):
    tst_sphere()
    tst2()
  print time.time()-start
  print "OK"
