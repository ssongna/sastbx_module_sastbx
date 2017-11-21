from scitbx import fftpack
from scitbx.array_family import flex
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from scitbx.python_utils import random_transform
import scitbx.math
from cStringIO import StringIO
import math
from sastbx import intensity



class parabola_base(object):
  def __init__(self, k=None ):
    p="^"+str(k)
    if k is None:
      p=" "
    self.base_name="parabola%s"%p
    self.k = k
  def base(self, x):
    x_min = flex.min(x)
    x_max = flex.max(x)
    assert(x_min>=-1)
    assert(x_max<=1)
    base_function = (1 - x*x)*3.0/4.0
    if self.k is not None:
      base_function = flex.pow( base_function, k )
    return base_function

class sphere_base(object):
  def __init__(self):
    self.base_name="sphere"
  def base(self,x):
    x_min = flex.min(x)
    x_max = flex.max(x)
    assert(x_min>=-1)
    assert(x_max<=1)
    r = (x+1)*0.5
    base_function = 6.0*r*r*(2-3.0*r+r*r*r)
    return base_function


class pofx(object):
  def __init__(self,prior,n,coefs=None,n_int=25):
    self.prior = prior
    self.coefs=coefs
    self.n=n
    self.base_name="pofx"
    if self.coefs is None:
      self.coefs = flex.double( [0]*self.n )
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)

    gle = scitbx.math.gauss_legendre_engine(n_int)
    self.x_int = gle.x()
    self.w_int = gle.w()
    self.normalize()

  def load_coefs(self, coefs=None,random_scale=1.0, thres=600):
    if coefs is not None:
      assert len(self.coefs)==self.n
      self.coefs = coefs
    else:
      self.coefs = random_scale*( flex.random_double( self.n )*2-1.0 )
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)
    self.normalize()

  def base(self,x, thres=600):
    x_min = flex.min(x)
    x_max = flex.max(x)
    assert(x_min>=-1)
    assert(x_max<=1)
    modi = self.polynome.f( x )
    sel = flex.bool(modi > thres)
    modi = modi.set_selected(sel, thres)
    base_function = flex.exp( modi )
    base_function = base_function*self.prior.base( x )
    return base_function

  def compute_normalisation_constant(self):
    px = self.base(self.x_int)
    result = px*self.w_int
    result = flex.sum(result)
    return result

  def normalize(self):
    result = self.compute_normalisation_constant()
    self.coefs[0]=self.coefs[0]-2.0*math.log(result)
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)
    result = self.compute_normalisation_constant()
    if abs(result-1.0)>1e-6:
      self.normalize()


  def entropy_simple(self):
    # chebyshev moment of posterior
    c_m_post = self.polynome.f( self.x_int )*self.base(self.x_int)
    c_m_post = c_m_post*self.w_int
    c_m_post = flex.sum(c_m_post)

    # chebyshev moment of prior
    c_m_prior = -self.polynome.f( self.x_int )*self.prior.base( self.x_int )
    c_m_prior = c_m_prior*self.w_int
    c_m_prior = flex.sum(c_m_prior)
    return c_m_post+c_m_prior, c_m_post, c_m_prior


  def relative_entropy(self, other_pofx):
    this_one = self.base(self.x_int)
    that_one = other_pofx.base(self.x_int)
    this_one_log = flex.log( this_one + 1e-12 )
    that_one_log = flex.log( that_one + 1e-12 )
    this_that = this_one*(this_one_log-that_one_log)
    this_that = this_that*self.w_int
    this_that = flex.sum( this_that )


    that_this = that_one*(-this_one_log+that_one_log)
    that_this = that_this*self.w_int
    that_this = flex.sum( that_this )
    return that_this + this_that

  def show(self):
    x  = self.x_int
    y  = self.prior.base(x)
    yy = self.base(x)
    ff = self.polynome.f( x )
    for a,b,c,f in zip(x,y,yy,ff):
      print a,b,c,f

class pofr(object):
  def __init__(self, dmax, n_params, prior=None,scale=1.0,n_int=24,m_int=24):
    self.dmax = dmax
    self.n_params=n_params
    self.prior=prior
    self.n_int=n_int
    self.m_int=m_int
    if self.prior is None:
      self.prior =  sphere_base()
    self.pofx = None
    self.pofx = pofx(self.prior,n_params,n_int=n_int)

    gle = scitbx.math.gauss_legendre_engine(self.m_int)
    self.x_int = gle.x()
    self.w_int = gle.w()
    self.r_intg = (self.x_int+1)*0.5*self.dmax

    self.r_block  = flex.double( range(int(dmax)*4 ))/4.0
    self.coefs = flex.double([0]*n_params)
    self.pofx.load_coefs(self.coefs)
    self.scale=scale

  def entropy_simple(self):
    return self.pofx.entropy_simple()

  def get_pofr(self,r):
    return self.pofx.base(r/self.dmax*2.0-1.0)*self.dmax/2.0

  def get_intensity(self, q, eps=1e-3):
    if q < eps:
      q = eps
    p_of_r = self.pofx.base( self.x_int )
    sinc_function = self.get_sinc( q, self.r_intg )
    tbi = p_of_r*sinc_function
    wtbi = tbi*self.w_int
    result = flex.sum(wtbi)
    return result*self.scale

  def get_sinc(self, q, r ):
    sinc = flex.sin( r*q )/(r*q)
    return sinc

  def update(self, params, add_scale=True):
    # the parameters we update are only coefs[0:]; however coefs[0] is strictly for normalisatuion purposes
    # the parameter named scale takes care of scaling against the intensity
    if add_scale:
      self.scale = abs(float(params[0])) # scale has to be positive please
    self.coefs = params.deep_copy()
    #now update pofxi, but reset the first value, it will be determined by normalisation
    self.coefs[0]=0.0
    self.pofx.load_coefs(self.coefs)


  def i(self, q_array):
      result = flex.double()
      for qq in q_array:
        result.append( self.get_intensity(qq) )
      return result

  def fast_i(self, int_object):
    return int_object.get_intensity( self )

  def f(self,r):
    x = r*2.0/self.dmax-1.0
    return self.pofx.base( x )/(self.dmax/2.0)

  def get_rg(self):
    r = self.dmax*flex.double(range(1,101))/100.0
    pr = self.f(r)
    rg2 = flex.sum(flex.pow(r,2.0)*pr)
    norma = flex.sum(pr)
    rg = math.sqrt(rg2/norma)/1.414
    return rg


class fast_integrator(object):
  def __init__(self, dmax,q,div=1,q_step=0.01):
    self.r = flex.double( range(int(dmax)*div) )/float(div) + 1.0/(div*2.0)
    self.dmax = dmax
    self.q = q
    self.div = div
    self.q_max=flex.max(q)+0.1
    self.q_step = q_step
    self.bli = intensity.block_integrator(self.dmax,1.0/div,self.q_max,self.q_step )
    self.bli.setup_arrays(self.q)

  def get_intensity(self, pr):
    ints = self.bli.get_intensity( pr.f(self.r) )#*2.0/self.dmax
    return ints*pr.scale

  def get_intensity_from_pr_array(self, pr_array):
    ints = self.bli.get_intensity( pr_array )#*2.0/self.dmax
    return ints



def test_fast_integrator(dmax=50):
  # make a pofx of a sphere
  import time
  pr = pofr(dmax,10)
  q = flex.double( range(350) )/1000.0
  fi = fast_integrator(dmax, q, div=1,q_step=0.01)
  a = time.time()
  i_slow = None
  i_fast = None
  for ii in range(100):
    i_fast = fi.get_intensity( pr )
  aa = time.time()
  for ii in range(100):
    i_slow = pr.i(q)
  aaa = time.time()
  print "#", aa-a, aaa-aa, ( aaa-aa)/ (aa-a)
  #for q, i,ii in zip(q,i_slow,i_fast):
  #  print q,i,ii






class pofr_variances(object):
  def __init__(self, pr, q):
    self.pr = pr
    self.q = q

  def perturb(self,variance_array, ntrials=100):
    assert self.pr.n_params == len(variance_array)
    intensity = None
    var = None
    calc_i = None
    for trial in xrange(ntrials):
      gauss = random_transform.normal_variate(N=self.pr.n_params)
      new_coefs = self.pr.coefs + gauss*variance_array
      new_pr = pofr(self.pr.dmax, self.pr.n_params, self.pr.prior, self.pr.scale,self.pr.n_int,self.pr.m_int)
      new_pr.update( new_coefs, add_scale=False)
      if intensity is None:
        tmp = new_pr.i( self.q )
        intensity = tmp.deep_copy()
        var = (tmp*tmp).deep_copy()
        calc_i = tmp.deep_copy()
      else:
        tmp = new_pr.i( self.q )
        intensity += tmp
        var += tmp*tmp

    mean_i = intensity/ntrials
    var_i = var/ntrials-mean_i*mean_i

    return calc_i, mean_i, var_i,  mean_i/calc_i, var_i/mean_i


class coefficient_sigma_model(object):
  def __init__(self, n=9, start=0.45, stop=0.04):
    self.n = n
    self.start = start
    self.stop=stop
    self.indices = flex.double( range(self.n) )
    self.sigma = self.compute_initial_sigma_array()

  def compute_initial_sigma_array(self):
    dy = self.start-self.stop
    dx = self.n+1
    rc = dy/dx
    result = self.indices*rc + self.start
    return result

  def damp_sigma(self, a, b):
    g = 1.0/( 1 + flex.exp( -(self.indices-a)*b ) )
    return self.sigma*g


class estimate_sigma_model(object):
  def __init__(self,obs_data, calc_data, pr):
    self.q = obs_data.q
    self.o = obs_data.i
    self.c = calc_data
    self.v = obs_data.s
    self.v = v*v
    self.pr = pr
    self.coef_model = coefficient_sigma_model(n=self.pr.n, start=0.45, stop=0.04)
    self.pr_vars = pofr_variances(self.pr, self.q)

    self.x = None
    self.n = 3
    self.domain = [(0,0.2), (-1,10), (1,50) ]
    self.optimizer = de.differential_evolution_optimizer(self, population_size=6, n_cross=1, f=0.80, eps=1e-5, monitor_cycle=50, show_progress=True)

  def target(self,vector):
    scale = abs(vector[0])
    a = vector[1]
    b = vector[2]
    # make a set of variances please
    sig_scales = self.coef_model.damp_sigma(a,b)
    sigs = self.pr_vars.damp(sig_scales,50)
    score = flex.abs(self.o-self.c*scale)/flex.sqrt(self.s*self.s + sigs*sigs)
    score = flex.sum(score)/score.size()
    return score

  def print_status(self,ms,mms,vector,txt):
    print ms, mms, list(vector), txt

  def get_solution(self):
    return abs(self.x[0]), self.x[1], self.x[2]





def test_pofx(n=8):
  p = parabola_base()
  s = sphere_base()

  pofx_p = pofx(p,n=12,n_int=45)
  pofx_s = pofx(s,n=12,n_int=45)


  for ii in xrange(10):
    s = ii/5.0
    for jj in range(10):
      pofx_p.load_coefs(None,s)
      pofx_s.load_coefs(None,s)
      t,a,b =  pofx_p.entropy_simple()
      assert(a>-1e-6)
      assert(b>-1e-6)
      t,a,b =  pofx_s.entropy_simple()
      assert(a>-1e-6)
      assert(b>-1e-6)

  this_pofr = pofr(45.0, 3)
  q_array = flex.double( range(150) )/150.0
  q_array[0]=0.0001
  i_array = []
  for q in q_array:
    tmp_i = this_pofr.get_intensity(q)
    i_array.append( tmp_i )
  ii_array = this_pofr.i( q_array )



if (__name__ == "__main__"):
  test_pofx()
  test_fast_integrator()
