from scitbx.array_family import flex
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
import scitbx.math
from cStringIO import StringIO
import math


def parabola_powered(x, k=None):
  x_min = flex.min(x)
  x_max = flex.max(x)
  assert(x_min>=-1)
  assert(x_max<=1)
  base = 1 - x*x
  if k is not None:
    base = flex.pow( base, k )
  return base


def sphere_pr(x):
  x_min = flex.min(x)
  x_max = flex.max(x)
  assert(x_min>=-1)
  assert(x_max<=1)
  r = (x+1)*0.5
  base = 12.0*r*r*(2-3.0*r+r*r*r)
  return base

base_types = {"sphere":sphere_pr, "parabola":parabola_powered }

def base_generator(base_type="sphere"):
  assert base_types.has_key( base_type )
  return base_types[ base_type ]



class function(object):
  def __init__(self, n, m=100, k=2.5, d_max=45.0, base_type="sphere"):
    self.n = n
    self.m = m
    self.k = k
    self.base_type = base_type
    self.d_max=d_max
    self.x = 1.0-2.0*(flex.double(range(m+1))/m)
    self.r = 0.5*(1+self.x)*self.d_max
    self.r[0] = 1e-8

    self.coefs = (flex.random_double(self.n)-0.5)*0.0
    self.load_coefs()
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)

  def show(self):
    result = get_p_of_r(self.x)
    for r,y in zip(self.r, result):
      print r, y

  def load_coefs(self, coefs=None):
    if coefs is None:
      self.coefs = (flex.random_double(self.n)-0.5)*2.0
    else:
      assert len(coefs)==self.n
      self.coefs = coefs
    # no means to refresh the coefficients yet in an elegant manner
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)


  def get_p_of_r(self,x):
    base =  base_generator(self.base_type)
    base = base(x)
    exp_pol = flex.exp( -self.polynome.f( x ) )
    result = exp_pol*base
    return result

  def get_sinc(self, q, x ):
    r = 0.5*(x+1)*self.d_max
    sinc = flex.sin( r*q )/(r*q)
    return sinc

  def integrate(self, q, ni):
    gle = scitbx.math.gauss_legendre_engine(ni)
    x_int = gle.x()
    w_int = gle.w()
    p_of_r = self.get_p_of_r(x_int)
    sinc = self.get_sinc( q, x_int )
    tbi = p_of_r*sinc
    wtbi = tbi*w_int
    result = flex.sum(wtbi)
    return result

def example():
  f = function(5,100)
  q_trials = [0.001, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

  ref_integrals = []
  for q in q_trials:
    ref_integrals.append( f.integrate(q,90) )

  for q, jj in zip(q_trials, range(len(q_trials))):
    print q,jj,
    for ii in range(2,90):
      print 100.0*abs(f.integrate(q,ii)-ref_integrals[jj])/abs(ref_integrals[jj]+1e-13),
    print





if (__name__ == "__main__"):
  example()
