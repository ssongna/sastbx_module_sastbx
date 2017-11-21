from scitbx.array_family import flex
from scitbx.math import log_factorial_generator
from stdlib import math as smath

class coefs(object):
  def __init__(self, n_max, only_even=True):
    self.n_max = n_max
    self.lgf = log_factorial_generator(n_max*2+2)
    self.dmap = {}
    self.emap = {}
    self.gen_coefs_D()
    self.gen_coefs_E()    

  def is_even(self,k):
    if k%2 == 0:
      return True
    return False

  def Dnmk(self,n,m,k):
    tmp1 = long(1.0)
    if k%2 != 0:
      tmp1=-1.0
    tmp2   = self.lgf.log_fact(n-k)
    tmp3_1 = self.lgf.log_fact(k)
    tmp3_2 = self.lgf.log_fact( (n+m)/2-k)
    tmp3_3 = self.lgf.log_fact( (n-m)/2-k)
    tmp3 = tmp3_1+tmp3_2+tmp3_3
    result = tmp1*smath.exp(tmp2-tmp3)
    return result

  def gen_coefs_D(self):
    # it is slow, but very easy: make a dictionairy
    for nn in range(self.n_max+1):
      for mm in range(nn+1):
        if (self.is_even(nn) == self.is_even(mm)):
          for kk in range((nn-mm)/2+1):
            self.dmap[ (nn,mm,kk) ] = self.Dnmk(nn,mm,kk)
            print nn,mm,kk, self.dmap[ (nn,mm,kk) ], "D COEF"

  def Ennnm(self,n,np,npp,m):
    result = long(0)
    for k in range((n-m)/2+1):
      for kp in range((np-m)/2+1):
        for kpp in range((npp-m)/2+1):
          result += self.dmap[ (n,m,k) ]*self.dmap[ (np,m,kp) ]*self.dmap[ (npp,m,kpp) ]/(n+np+npp-2.0*(k+kp+kpp)+1)
    return result

  def sort_indices(self,a,b,c):
    ia = flex.int( [a,b,c] )
    so = flex.sort_permutation( ia,True )
    ia = ia.select( so )
    return ia[0],ia[1],ia[2]


  def gen_coefs_E(self):
    n_max = self.n_max
    for n in range(n_max+1):
      for np in range(n+1):
        if self.is_even(n) == self.is_even(np):
          for npp in range(np+1):
            if self.is_even(np) == self.is_even(npp):
              for m in range( min(n,np,npp)+1 ):
                if self.is_even(m) == self.is_even(npp):
                  self.emap[ (n,np,npp,m) ] = self.Ennnm(n,np,npp,m)
                  print n,np,npp,m , self.emap[ (n,np,npp,m) ] 
           


def build_tst_image():
   




def tst():
  nmax=40
  cd = coefs(nmax)


if __name__ == """__main__""":
   tst()









      
   







