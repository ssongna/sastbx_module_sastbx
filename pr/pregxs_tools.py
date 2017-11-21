import sys, os, math
from scitbx import simplex
from scitbx import differential_evolution as de
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
from scitbx.array_family import flex
from sastbx.pr import pr_tools
import random,time
from sastbx.data_reduction import saxs_read_write
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from libtbx import easy_pickle
import time
from sastbx.data_reduction import curves
from sastbx.basic_analysis import guinier_analyses


def get_q_array_uniform_body(data, level=1.0e-2, n_last=5, q_min=0.0, q_max=0.25):
  # first we need to define a 'background level'
  n_points = data.q.size()
  #bkg = flex.mean(data.i[n_points-n_last:])
  bkg = 0
  i_tmp = data.i-bkg
  location = n_points
  maxi = flex.max( data.i )-bkg
  for ii in range(n_points):
    this_i = data.i[n_points-1 - ii] - bkg
#    print data.q[n_points-1 - ii],this_i,this_i/maxi,level
    if this_i >= maxi*level:
      location = ii
      break

  # now cut the data please
  location = n_points-1-location
  if q_max > data.q[location]:
    q_max = data.q[location]
  new_data = data.cut_data(q_min,q_max)
  return new_data








class multi_fit(object):
  def __init__(self, data, d_max, n_params=5):
    self.data     = data
    self.d_max    = d_max
    self.n_params = n_params

    self.pofr = pr_tools.pofr(self.d_max, self.n_params)
    self.integrator_object = pr_tools.fast_integrator(q=data.q, dmax=d_max)

    self.matrix = [ flex.random_double(n_params)*2.0-0.5 ]
    for ii in range(n_params):
      self.matrix.append( flex.random_double(n_params)*2.0-0.5 )


    self.opti = dssa(n_params, self.matrix, self)

  def target(self, vector):
    self.pofr.update( vector )
    calc_data = self.pofr.fast_i( self.integrator_object )
    score = calc_data - self.data.i
    score = flex.pow2(score/self.data.i)
    #score = flex.pow2(score)
    score = flex.sum(score)/calc_data.size()
    return score










class fixed_dmax_fitter_with_np(object):
  def __init__(self, prior, n_params, d_max, data, integrator, n_fst_pass=4,simplex_trial=5):
    self.prior = prior
    self.n_params = n_params
    self.n_fst_pass = n_fst_pass
    self.d_max = d_max
    self.data = data
    self.integrator_object = integrator
    self.integrator_object.bli.setup_arrays(data.q)
    self.simplex_trial = simplex_trial

    self.x = None
    self.f = None
    t1 = time.time()
#    self.counter=0
    if(self.n_params < 6):
      self.de()
    self.simplex()
    self.calc_entropy()
#   update r, pr for the final solution
    self.r = flex.double( range(int(self.d_max)) )
    self.best_pr = self.pofr.f( self.r )
    self.f = self.best_pr
    self.calc_i = self.pofr.fast_i( self.integrator_object )
#    print self.counter, time.time() - t1

  def de(self):
    self.n = self.n_fst_pass
    self.pofr = pr_tools.pofr(self.d_max, self.n, self.prior)
    self.domain = [ (-0.1,0.1) ]*self.n
    self.optimizer = de.differential_evolution_optimizer( self,
                                                        population_size=self.n,
                                                        n_cross=1,
                                                        f=0.85,
                                                        eps=1e-5,
                                                        monitor_cycle=50,
                                                        show_progress=False)
    # self.x would have been updated by this optimizer


  def simplex(self):
    self.n = self.n_params

    if(self.x is None):
      self.x = flex.double([1]+[0]*(self.n_params-1))
    else:
      self.x = self.x.concatenate( flex.double([0]*(self.n_params - self.n_fst_pass)) )
    self.pofr = pr_tools.pofr(self.d_max,self.n,self.prior)

    self.simplex_scores = []
    self.simplex_solutions = []
    for ii in xrange(self.simplex_trial):
      #make a random simplex
      self.starting_simplex = []
      for ii in range(self.n):
        self.starting_simplex.append((flex.random_double(self.n)*2-1.0)+self.x)
      self.starting_simplex.append(self.x)

      self.optimizer = simplex.simplex_opt( dimension=self.n,
                                    matrix  = self.starting_simplex,
                                    evaluator = self,
                                    tolerance=1e-4)

      self.solution = self.optimizer.get_solution()
      self.score = self.target( self.solution )
      self.simplex_scores.append( self.score )
      self.simplex_solutions.append( self.solution )

    best_index = flex.min_index( flex.double(self.simplex_scores) )
    self.solution = self.simplex_solutions[ best_index ]
    #self.cma(m=self.solution)
    self.score = self.simplex_scores[ best_index ]
    self.pofr.update( self.solution )

  def cma(self, m=None):
    if(m is None):
      m = flex.double([1]*self.n_params)
    s = flex.double([1]*self.n_params)
    cma_obj = cma_es_driver( self.n_params, m, s, self.target )
    self.pofr.update( cma_obj.x_final )
    self.solution = cma_obj.x_final
    self.score = self.target( cma_obj.x_final )


  def target(self, vector):
    self.pofr.update( vector )
    calc_data = self.pofr.fast_i( self.integrator_object )
    score = calc_data - self.data.i
    score = flex.pow2(score /self.data.i)
    score = flex.sum(score)/calc_data.size()
#    self.counter +=1
    return score

  def show_pr(self,outfile):
    for rr, ff in zip(self.r,self.best_pr):
      if(ff > 10e-5):
        print >> outfile, rr, ff

  def show_obs_vs_calc(self,outfile):
    calc = self.pofr.fast_i( self.integrator_object )
    q = self.data.q
    obs = self.data.i
    for qq,cc,oo in zip(q,calc,obs):
       print >> outfile, qq, oo, cc

  def get_cdf(self):
    if(self.f is None):
      r = flex.double( range(int(self.d_max)) )
      self.f = self.pofr.f( r )
    cdf = flex.double()
    sum = 0.0
    cdf.append(sum)
    for ff in self.f:
      sum += ff
      cdf.append(sum)
    return cdf

  def calc_entropy(self):
    self.entropy = self.pofr.pofx.relative_entropy(pr_tools.sphere_base())


#========================================================================#
class fixed_dmax_fitter(object):
  def __init__(self, prior, data, d_max, n_params, n_fst_pass, n_trial=1, n_simplex=10):
    self.prior = prior
    self.data = data
    self.d_max = d_max
    self.n_params = n_params
    self.n_fst_pass = n_fst_pass
    self.n_trial = n_trial
    self.n_simplex = n_simplex
    self.integrator = pr_tools.fast_integrator(self.d_max, self.data.q, div=1)
    self.trials = []

    self.sub_data()
    self.fit()
    self.get_best_fit()

  def sub_data(self):
#    self.increment=math.pi / self.d_max
    self.q_min = self.data.q[0]
    self.q_max0 = self.data.q[self.data.q.size() - 1]
    self.increment=(self.q_max0 - self.q_min) / self.n_params
    self.increment=max(self.increment, math.pi / self.d_max)
    self.data_array = []
    for nn in range(self.n_fst_pass, self.n_params):
      q_max = min(self.q_min + self.increment * nn, self.q_max0)
      stop_at = q_stop(self.data.q, q_max)
      self.data_array.append(curves.simple_saxs_data(self.data.q[:stop_at], self.data.i[:stop_at], self.data.s[:stop_at]))
    self.data_array.append(self.data)


  def fit(self):
    for ii in range(self.n_trial):
      prior = self.prior
      sphere = pr_tools.sphere_base()
      for nn in range(self.n_fst_pass, self.n_params+1):
        tmp_object = fixed_dmax_fitter_with_np(prior, nn, self.d_max, self.data_array[nn-self.n_fst_pass], self.integrator,
                                               n_fst_pass=self.n_fst_pass,
                                               simplex_trial=self.n_simplex)
        prior = update_prior(tmp_object.pofr.pofx, sphere)
      self.trials.append( tmp_object )

  def get_best_fit(self):
    scores = []
    for ii in range(self.n_trial):
      scores.append( self.trials[ii].score )
    best_index = flex.min_index( flex.double(scores) )
    self.best_fit = self.trials[ best_index ]
    self.best_score = scores[ best_index ]
    ### update best fitting info for print out ###
    self.calc_i = self.trials[best_index].calc_i
    self.r = self.trials[best_index].r
    self.pr = self.trials[best_index].best_pr
    ### end of update ####



#=============================================================================#
# this class implements the dmax_scanning
#=============================================================================#
class dmax_scan(object):
  def __init__(self, prior, data, d_max, delta, step, rg, n_params, n_fst_pass, n_trial=1, n_simplex=10, entropy_thresh=1.24,outfile=''):
    self.prior = prior
    self.data = data
    self.d_max = d_max
    self.rg = rg
    self.n_params = n_params
    self.n_fst_pass = n_fst_pass
    self.n_trial = n_trial
    self.n_simplex = n_simplex
    self.entropy_threshold=entropy_thresh
    self.d_array=[]
    self.outfile = outfile 
    d_start = d_max - delta
    d_end = d_max + delta
    self.n_step = int(2.0 * delta / step + 0.5)
    for ii in range(self.n_step):
      d = d_start + ii * step
      self.d_array.append(d)

    self.scan()
    self.get_ks_dist()
    self.get_pr_statistics()

  def scan(self):
    self.fitters = []
    self.cdfs = []
    self.prs = []
    self.entropy = flex.double()
    for d in self.d_array:
      fitter_d = fixed_dmax_fitter(self.prior,
                                   self.data,
                                   d,
                                   self.n_params-1,
                                  # self.n_fst_pass+3,
                                   self.n_fst_pass,
                                   self.n_trial,
                                   self.n_simplex)
      self.entropy.append( fitter_d.best_fit.entropy )
      if(fitter_d.best_fit.entropy > self.entropy_threshold):
        print  "Attention: "
        print  "the shape of the molecule might be non-globular;"
        print  "There might be (partial) unfolding, if it is a protein"
        print  "Given dmax="+str(d)

        with open(self.outfile,"a") as log:
          print >>log, "Attention: "
          print >>log, "the shape of the molecule might be non-globular;"
          print >>log, "There might be (partial) unfolding, if it is a protein"
          print >>log, "Given dmax="+str(d)



      self.fitters.append( fitter_d )
      tmp_cdf = fitter_d.best_fit.get_cdf()
      self.cdfs.append( tmp_cdf.concatenate( flex.double([1]*int(self.d_array[-1]-d)) ) )
      tmp_pr = fitter_d.best_fit.f
      self.prs.append( tmp_pr.concatenate( flex.double([0]*int(self.d_array[-1]-d)) ) )


  def get_ks_dist(self):
    self.dist_mat=[]
    n_trials = len(self.fitters)
    self.saved_trials=n_trials
    for ii in range(n_trials):
      self.dist_mat.append( flex.double( [0] * n_trials) )

    for ii in range(n_trials):
      for jj in range(ii):
        d_cdf = self.cdfs[ii] - self.cdfs[jj]
        max_d_cdf = flex.max( flex.abs( d_cdf ) )
        self.dist_mat[ii][jj] = max_d_cdf
        self.dist_mat[jj][ii] = max_d_cdf

    average_mcdf=flex.double()
    for ii in range(n_trials):
      average_mcdf.append( flex.mean( self.dist_mat[ii] ) )

    self.best_index = flex.min_index( average_mcdf )
    self.dmax_best=self.d_array[ self.best_index ]
    self.mcdf_mean = average_mcdf[ self.best_index ]
    self.mcdf_var = flex.mean( flex.pow2( self.dist_mat[self.best_index] - self.mcdf_mean ) )
    self.mcdf_sigma = math.sqrt( self.mcdf_var )


  def get_pr_statistics(self):
    sel = flex.bool( self.dist_mat[ self.best_index ] < self.mcdf_mean )

    self.average_pr = flex.double(int(self.d_array[-1]),0)
    tmp_pr2 = self.average_pr*0

    refined_d_array = flex.double()
    for ii in range(self.saved_trials):
      if(sel[ii]):
        self.average_pr = self.average_pr + self.prs[ii]
        tmp_pr2 = tmp_pr2 + flex.pow2( self.prs[ii] )
        refined_d_array.append( self.d_array[ii] )

    #print refined_d_array.size()
    self.average_pr = self.average_pr / refined_d_array.size()
    self.sigma_pr = flex.sqrt( (tmp_pr2 / refined_d_array.size() - flex.pow2(self.average_pr) ) )
    self.d_mean = flex.mean(refined_d_array)
    self.d_sigma = math.sqrt( flex.mean( flex.pow2( refined_d_array - self.d_mean ) ) )

    print  "DMAX(based on k-s distance)=",
    print "BEST DMAX = ", self.d_array[self.best_index]
    print "MEAN DMAX(sigma) = ", self.d_mean, "("+str(self.d_sigma)+")"
    
    with open (self.outfile,"a") as log:
      print >>log, "DMAX(based on k-s distance)=",
      print >>log,"BEST DMAX = ", self.d_array[self.best_index]
      print >>log,"MEAN DMAX(sigma) = ", self.d_mean, "("+str(self.d_sigma)+")"
    self.r = flex.double(range(self.average_pr.size() ))

  def print_pr(self, outfile, threshold=1e-6):
    for rr, pp,ss in zip(self.r,self.average_pr,self.sigma_pr):
      if(pp > threshold):
        print >> outfile, rr,pp,ss


  def get_scores(self):
    scores = []
    rg_s = []
    for fitter in self.fitters:
      rg = fitter.best_fit.pofr.get_rg()
      s = abs( self.rg - rg )
      rg_s.append(rg)
      scores.append(s)
    self.best_index = flex.min_index( flex.double(scores) )
    self.dmax_best=self.d_array[ self.best_index ]
    self.calc_i = self.fitters[self.best_index].calc_i
    print "#DMAX=", self.dmax_best, self.d_max, self.rg, rg_s[self.best_index]

  def get_best_dmax(self):
    prior = self.fitters[self.best_index].best_fit.pofr.pofx
    fitter = fixed_dmax_fitter(prior,
                                   self.data,
                                   self.dmax_best,
                                   self.n_params,
                                   self.n_params,
                                   self.n_trial,
                                   self.n_simplex)
    #fitter.best_fit.show_pr( open("best_scan.pr", 'w') )
    #fitter.best_fit.show_obs_vs_calc( open("best_scan.qii", 'w') )
    self.fitter = fitter
    self.calc_i = fitter.best_fit.calc_i
    threshold = 1e-4
    dmax_best0 = self.dmax_best
    with open(self.outfile,"a") as log:
      print >>log,"DMAX_CORRECTED:" 
      print >>log,"dmax(threshold)"
    for rr in range(int(dmax_best0/2.0),int(dmax_best0) ):
      if(fitter.best_fit.f[rr] < threshold):
        self.dmax_best = rr
        with open(self.outfile,"a") as log:
          print >>log, str(self.dmax_best)+"("+str(threshold)+")"
        threshold = threshold / 2.0
        if(threshold < 1e-5):
          break


#===============================================================================#
# Global utility functions
#===============================================================================#
def q_stop(q_array, stop_value):
   for ii in range(q_array.size()):
     if(q_array[ii] > stop_value):
       return ii

def update_prior(prior,sphere):
  if(prior.prior.base_name == 'pofx'):
    prior_coefs = prior.prior.coefs
    self_coefs = prior.coefs
    diff_len = len(self_coefs) - len(prior_coefs)
    if(diff_len > 0):
      coefs = self_coefs + prior_coefs.concatenate(flex.double( [0]*diff_len ))
    else: # diff_len = 0
      coefs = self_coefs + prior_coefs
    prior = pr_tools.pofx(sphere, len(self_coefs), coefs)
  return prior

#===============================================================================#


def get_inputs(args):
  data = args[0]
  d_max = float(args[1])
  if(len(args)>2):
    rg = float(args[2])
  else:
    rg = 0
  prior = None
  data = saxs_read_write.read_standard_ascii_qis(data)
  if(rg == 0):
    msga = guinier_analyses.multi_step_rg_engine( data )
    rg = msga.median_rg

  m = 1.0/data.i[0]
  data.multiply_add(m,0.0)
  n_params=10
  n_fst_pass=4

#  fitter = fixed_dmax_fitter(prior, data, d_max, n_params, n_fst_pass, n_trial=4, n_simplex=10)
#  fitter.best_fit.show_pr( open("best.pr",'w') )
#  fitter.best_fit.show_obs_vs_calc( open("best.qii",'w') )

  delta=rg
  step=2
  d_max_scan = dmax_scan(prior, data, d_max, delta, step, rg, n_params, n_fst_pass, n_trial=1, n_simplex=10)
  d_max_scan.get_best_dmax()


if __name__ == "__main__":
  t1 = time.time()
  get_inputs( sys.argv[1:] )
  print time.time() - t1
