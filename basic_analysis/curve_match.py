from scitbx.array_family import flex
import sys, math

from scitbx import simplex
from scitbx.math import chebyshev_polynome
from sastbx.basic_analysis import pofr_data

class nonlinear_fit(object):
  def __init__(self, x, y, order):
    self.x = x
    self.y = y
    self.n = order
    coefs = flex.random_double( self.n )
    self.polynome = chebyshev_polynome( self.n, -1.0, 1.0, coefs )
    self.optimize()

  def optimize(self):
    start_matrix = []
    for ii in range( self.n + 1 ):
      start_matrix.append( flex.random_double(self.n) )

    optimizer = simplex.simplex_opt(dimension = self.n,
				    matrix = start_matrix,
				    evaluator = self,
				    tolerance=1e-4
				   )
    self.solution = optimizer.get_solution()
    self.score = self.target( self.solution )

  def target(self, vector):
    polynome = chebyshev_polynome( self.n, -1.0, 1.0, vector )
    self.fit_y = polynome.f( self.x )
    score = flex.sum( flex.pow2( self.fit_y - self.y ) )
    return score

def calibration( pr_1, pr_2, n_params):
  x1 = pr_1.r/pr_1.r[-1]
  x2 = pr_2.r/pr_2.r[-1]
  cdf_1 = pr_1.pr2cdf()
  cdf_2 = pr_2.pr2cdf()
  new_cdf_2 = flex.linear_interpolation( x2, cdf_2, x1).deep_copy()
  fitting = nonlinear_fit(  new_cdf_2, cdf_1, n_params )
  solution = fitting.solution
  fn = chebyshev_polynome( n_params, -1.0, 1.0, solution)
  return fn

def test(args):
  pr_1 = pofr_data.pofr(filename=args[0])
  #pr_2 = pofr_data.pofr(filename=args[1])
  pr_2 = pofr_data.pofr(filename=args[1], from_gnom=False)
  num_params = int(args[2])
  #pr_1.print_data(pr_1.pr)
  #pr_2.print_data(pr_2.pr)
  new_r = flex.double( range(101) ) /100.0

  pr_1.pr = pr_1.linear_interpolation( new_r * pr_1.r[-1])
  pr_2.pr = pr_2.linear_interpolation( new_r * pr_2.r[-1])
  pr_1.r = new_r * pr_1.r[-1]
  pr_2.r = new_r * pr_2.r[-1]

  fn = calibration( pr_1, pr_2, num_params )
  cdf_2 = pr_2.pr2cdf()
  calibrated_cdf_2 = fn.f( cdf_2 )

  pr_2.cdf = calibrated_cdf_2

  fitted_pr = pr_2.cdf2pr()
  fitted_pr = fitted_pr * pr_2.r[-1] / pr_1.r[-1]
  #new_pr_1 = pr_1.linear_interpolation( pr_2.r )
  new_pr_1 = pr_1.pr
  print '&'
  for ri, p1 in zip( pr_1.r, fitted_pr):
    print ri, p1
  pr_1.print_data( pr_1.pr)
  score = flex.sum(flex.pow2( fitted_pr - new_pr_1 ))
  score = math.sqrt( score / pr_2.r.size() )
  
  print '#',
  for s in fn.coefs():
    print s,
  print pr_1.r[-1], score


  #print '#',list(solution)
  #for xi,yi,fi in zip( fitting.x, fitting.y, fitting.fit_y):
  #  print xi,yi,fi
 

if __name__ == "__main__":
  test(sys.argv[1:])
