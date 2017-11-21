import math

from sastbx.fXS import fourier_legendre_series, set_negative_to_zero
from scitbx import lbfgsb
from scitbx.array_family import flex

class fit_parameters(object):
  def __init__(self):
    self.fls_data = None
    self.minimize = True
    self.standardize = False
    self.wavelength = 1.0
    self.fit_order = 40
    self.q = None
    self.x = None
    self.ac = None
    self.v = None

def read_ac(filename=None):

  f = open(filename, 'r')
  l = f.readlines()
  f.close()

  q = list()
  x = list()
  c2 = list()
  v = list()

  for i in xrange(len(l)):
    if (l[i][0] == '#'):
      q.append(l[i].rstrip())
      tmp_x = flex.double()
      tmp_c2 = flex.double()
      tmp_v = flex.double()
    elif (l[i][0] == '&'):
      x.append(tmp_x[6:-5].deep_copy())
      c2.append(tmp_c2[6:-5].deep_copy())
      v.append(tmp_v[6:-5].deep_copy())
    else:
      split_line = l[i].split()
      tmp_x.append(float(split_line[0]))
      tmp_c2.append(float(split_line[1]))
      if (len(split_line) > 2):
        s = float(split_line[2])
      else:
        s = 1.0e-6
      tmp_v.append(s*s)  # variance / n

  return q, x, c2, v

def read_sc(filename=None):

  split_filename = filename.split('.')
  new_filename = split_filename[0] + '.sc'

  f = open(new_filename, 'r')
  l = f.readlines()
  f.close()

  q = flex.double(len(l))
  sc = flex.double(len(l))

  for i in xrange(len(l)):
    split_line = l[i].split()
    q[i] = float(split_line[0])
    sc[i] = float(split_line[1])

  return q, sc

def set_odd_coefficients_to_zero(c):
  for i in xrange(1,len(c),2):
    c[i] = 0.0
  return c

class lbfgs_optimizer(object):
  def __init__(self,f,s,c,x,fls):
    self.f = f
    self.s = s
    self.xx = x
    self.x = set_negative_to_zero(c.deep_copy())
    self.n = len(c)
    self.fls = fls
    self.max_d = None
    self.dx = 1.0e-6
    self.tdx = 2.0*self.dx

    assert (len(self.f) == len(self.s))

    l = flex.double(self.n, 0)
    u = flex.double(self.n, 0)
    nbd = flex.int(self.n, 1)
    self.minimizer = lbfgsb.minimizer(n=self.n, m=20, l=l, u=u, nbd=nbd)
    self.fc, self.g = self.compute_functional_and_gradients()
    while True:
      if (self.minimizer.process(self.x, self.fc, self.g)):
        self.fc, self.g = self.compute_functional_and_gradients()
      elif (self.minimizer.is_terminated()):
        break

  def objective_function(self,x):
    c_new = flex.double(2*len(self.x)-1,0.0)
    count = 0
    for i in xrange(0,len(c_new),2):
      c_new[i] = x[count]
      count += 1
    f_new = self.fls.compute_function(c_new,self.xx)
    d = (self.f - f_new)/self.s
    d = flex.sum(d*d)
    return d

  def gradient(self,x):
    g = flex.double(self.n, 1.0)
    for i in xrange(self.n):
      new_x = x.deep_copy()
      new_x[i] = new_x[i] + self.dx
      f1 = self.objective_function(new_x)
      new_x = x.deep_copy()
      new_x[i] = new_x[i] - self.dx
      f0 = self.objective_function(new_x)
      g[i] = (f1 - f0)/self.tdx
    return g

  def compute_functional_and_gradients(self):
    new_x = self.x.deep_copy()
    obj_fun = self.objective_function(new_x)
    g = self.gradient(new_x)

    self.max_d = obj_fun

    return obj_fun, g

  def estimate_asymptotic_variance(self,x):
    v = flex.double(self.n,1.0)
    for i in xrange(self.n):
      new_x = x.deep_copy()
      new_x[i] = new_x[i] + self.dx
      g1 = self.gradient(new_x)
      new_x = x.deep_copy()
      new_x[i] = new_x[i] - self.dx
      g0 = self.gradient(new_x)
      d = g1[i] - g0[i]
      if (d > 0.0):
        v[i] = self.tdx/d
    return v

def fit_data(p):

  fls = fourier_legendre_series()
  fls.read_polynomials(p.fls_data)

  cos_sq = p.q*p.wavelength/(4.0*math.pi)
  cos_sq = cos_sq*cos_sq
  sin_sq = 1.0 - cos_sq
  fit_x = cos_sq + sin_sq*flex.cos(p.x)
  fit_ac = p.ac.deep_copy()
  fit_x, fit_ac = zip(*sorted(zip(fit_x,fit_ac)))
  fit_x = flex.double(fit_x)
  fit_ac = flex.double(fit_ac)

  fit_c = fls.compute_coefficients(p.fit_order,fit_ac,fit_x)
  fit_c = set_odd_coefficients_to_zero(fit_c)

  fit_v = fls.compute_coefficient_variances(p.fit_order,p.v,fit_x)
  fit_v = set_odd_coefficients_to_zero(fit_v)

  if (p.minimize):
    nz_c = flex.double()
    for k in xrange(0,len(fit_c),2):
      if (fit_c[k] < 0.0):
        fit_c[k] = -fit_c[k]
      nz_c.append(fit_c[k])
    m = lbfgs_optimizer(fit_ac,flex.sqrt(p.v),nz_c,fit_x,fls)
    nz_v = m.estimate_asymptotic_variance(m.x)
    assert(nz_c.size() == nz_v.size())
    count = 0
    for k in xrange(0,len(fit_c),2):
      fit_c[k] = m.x[count]
      fit_v[k] = nz_v[count]
      count += 1

  f_min = 1.0
  f_max = 1.0

  if (p.standardize):
    # standardize fitted curve to have a min of 1.0 and max of 2.0
    old_f = fls.compute_function(fit_c,fit_x)
    # assume f is positive
    f_min = flex.min(old_f)
    f_max = flex.max(flex.fabs(old_f / f_min - 1.0))
  scales = (f_min,f_max)

  return fit_c, fit_v, scales
