import sys, os, math
from scitbx import simplex
from scitbx import differential_evolution as de
#from scitbx import simulated_annealing
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
from scitbx.array_family import flex
import pr_tools
import random,time
from sastbx.data_reduction import saxs_read_write
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
import time

import entropic_restraints

class rcs_fitter(object):
  def __init__(self, n_params, n_fst_pass, d_max, data,  n_int=35, simplex_trial=5):
    self.n = n_fst_pass
    self.n_coeff = n_params
    self.delta = self.n_coeff - self.n
    self.data = data
    self.d_max = max(self.data.q)
    self.n_int = n_int
    self.x = None

    self.ent = entropic_restraints.entropy_restraint()

    # make a pofr please
    self.pofr = pr_tools.pofr(self.d_max,self.n,n_int=self.n_int, m_int=self.n_int)

    self.q_weight = flex.bool(self.data.q < 0.1 ).as_double()

    self.weight = 1.0
    # first we do a global optimisation using a diferentail evolution search
    self.domain = []
    for ii in range(self.n):
      self.domain.append( (-0.1, 0.1) )
    self.optimizer = de.differential_evolution_optimizer(self, population_size=self.n, n_cross=1, f=0.85, eps=1e-5, monitor_cycle=50, show_progress=False)    
    self.q_weight = self.q_weight*0.0+1.0 
    self.x = self.x.concatenate( flex.double([0]*self.delta) )
    self.n = self.n+self.delta
    self.pofr = pr_tools.pofr(self.d_max,self.n,n_int=self.n_int, m_int=self.n_int)
    self.simplex_trial=simplex_trial
    self.simplex_scores = []
    self.simplex_solutions = []
    for ii in xrange(self.simplex_trial):
      #make a random simplex please
      self.weight = 1.0
      self.starting_simplex = []
      for ii in range(self.n+1):
        self.starting_simplex.append(0.10*(flex.random_double(self.n)*2-1.0)+self.x)

      self.optimizer = simplex.simplex_opt( dimension=self.n,
                                    matrix  = self.starting_simplex,
                                    evaluator = self,
                                    tolerance=1e-3)     
      self.solution = self.optimizer.get_solution()
      self.score = self.target( self.solution )
      self.simplex_scores.append( self.score )
      self.simplex_solutions.append( self.solution )

    best_simplex_score = self.simplex_scores[0]
    this_simplex = 0    
    for ii in xrange(self.simplex_trial):
      if self.simplex_scores[ii] < best_simplex_score:
        best_simplex_score = self.simplex_scores[ii] 
        this_simplex = ii

    self.solution = self.simplex_solutions[ this_simplex ]

    #self.optimizer = simulated_annealing.sa_optimizer( self, self.solution, flex.double( self.n*[0.0051] ), start_t=2.1, end_t=0.001, burn_in=100, burn_out=50000, steps=5000 , show_progress=True)
    #self.solution, self.score = self.optimizer.get_solution()
    self.pofr.update( self.solution )

    self.calc_data = self.pofr.f( self.data.q )
        


  def print_status(self, mins, means, vec, txt):
    print mins, means, list(vec), txt
  

  def target(self, vector):
    self.pofr.update( vector )
    calc_data = self.pofr.f( self.data.q )
    score = flex.pow2((calc_data - self.data.i)/(calc_data+self.data.i))
    score = flex.sum(score)/calc_data.size()
    total_score = score 
    #print "#",total_score, t, score, self.ent.rstr(t), list(vector)
    return total_score

  def get_best_pofr(self):
    self.pofr.update( self.solution )
    return self.pofr



class random_start_fixed_dmax(object):
  def __init__(self, data, dmax, n_param, n_fst_pass, n_trial=10, n_simplex=10, data_q=None):
    self.trials = []
    self.data=data
    self.dmax=max(data.q)
    if(data_q is None):
      self.data_q = flex.double( range(51) )/ 100.0
    else:
      self.data_q = data_q
    self.integrator_obj = pr_tools.fast_integrator(self.dmax, self.data_q)
    self.new_r = flex.double(range(1,int(self.dmax+1)) )
    intepolated_pr = flex.linear_interpolation( data.q, data.i, self.new_r)
    #print list(intepolated_pr)
    self.i_expt_pr = self.integrator_obj.get_intensity_from_pr_array( intepolated_pr )
    
    for ii in range( n_trial ):
      tmp_object = rcs_fitter(n_param, n_fst_pass, dmax, data, simplex_trial=n_simplex)
      self.trials.append( tmp_object )
    self.collect_scores()


  def collect_scores(self):
    self.chi=[]
    self.ent=[]
    for oo in self.trials:
      c= oo.target(oo.solution)
      self.chi.append(c)
    self.chi_index = flex.min_index( flex.double(self.chi) )

  def print_entropy(self): 
    f = self.trials[ self.chi_index ]
    t,a,b = f.get_best_pofr().entropy_simple()
    print t,a,b,"entropy",

  def show_calc_i(self,outfile):
#    q = flex.double( range(1,51) )/100.0
    pr_fit = self.trials[self.chi_index].get_best_pofr().f( self.new_r )
    i_fit = self.integrator_obj.get_intensity_from_pr_array( pr_fit )
    for qq,ii,i2 in zip(self.data_q,self.i_expt_pr, i_fit):
      print >> outfile, qq,ii,i2

  def show_best_pr(self,outfile):
    r = self.dmax*flex.double( range(101) )/100.0
    f = self.trials[ self.chi_index ]
    f = f.get_best_pofr().f( r )
    for rr, ff in zip(r,f):
      print >> outfile, rr, ff

  def show_obs_vs_calc(self,outfile):
    calc = self.trials[ self.chi_index ].calc_data
    rg = self.trials[ self.chi_index ].get_best_pofr().get_rg()
    q = self.data.q
    obs = self.data.i
    print >>outfile, "& RG = ", rg, self.dmax
    for qq,cc,oo in zip(q,calc,obs):
       print >> outfile, qq, oo, cc 




class d_max_scan(object):
  def __init__( self, data, n_param, n_fst_pass, n_trial, d_max_start, d_max_stop, n_step=10, n_simplex=5):
    self.data = data
    self.n_param = n_param
    self.n_fst_pass = n_fst_pass
    self.n_trial = n_trial
    self.n_simplex=n_simplex
    delta = (d_max_stop-d_max_start)/n_step
    self.d_max_array = []
    for ii in range(n_step+1):
      self.d_max_array.append( d_max_start+ii*delta )

    self.chi_scores = []
    self.ent_scores = []
    self.objects = []

    for dmax in self.d_max_array:
      o,c,e = self.this_dmax(dmax)
      self.chi_scores.append( c )
      self.ent_scores.append( e )
      self.objects.append( o )
      print dmax, c, e 


  def this_dmax(self, dmax):
    oo = random_start_fixed_dmax( self.data, dmax, self.n_param, self.n_fst_pass, self.n_trial, self.n_simplex )
    chi = oo.chi[ oo.chi_index ]
    ent = oo.ent[ oo.chi_index ]
    return  oo, chi, ent
 










master_params = iotbx.phil.parse("""\
pregxs{
  data = None
  .type=path
  .help = "q Intensity Sigma"
  .multiple=True

  d_max = None
  .type=float
  .help="Maximum distance in particle"

  scan = False
  .help="When True, a dmax scan will be performed"
  .type=bool

   

 
  fitting{
    n_coeff = 8
    .type=int
    .help="Number of coefficients describing p(r)"
    n_fst_pass = 4
    .type = int
    .help="Number of coefficients fitted in the first optimisation pass."
    n_trials=2
    .type=int
    .help="n_trials fittings will be performed, only the best will be taken."
    n_trials_simplex=10
    .type=int
    .help="n_trials_simplex random restarts of simplex will be performed. make sure this aint 1"

    delta = 15
    .type=float
    .help="scan range:  dmax-delta - > dmax+delta"
    
    n_step = 6
    .type=int
    .help="Number of steps in dmax scan"

  }
  output = "pregxs"
  .type=path
  .help = "Output base. expect a .pr and a .qii file."

}
""")



def go(params, out=None):
  if out is None:
    out = sys.stdout

  data_array = []
  multies = []
  dmax=params.pregxs.d_max 
  nparam=params.pregxs.fitting.n_coeff
  nfst=params.pregxs.fitting.n_fst_pass
  ntrials=params.pregxs.fitting.n_trials
  strials=params.pregxs.fitting.n_trials_simplex
  for item in params.pregxs.data:
    data = saxs_read_write.read_standard_ascii_qis(item) 
  #  m = 1.0/data.i[0]
  #  data.multiply_add(m,0.0)
    data_array.append( data )
  #  multies.append( m )


    if params.pregxs.scan:
      d_max_start = dmax-params.pregxs.fitting.delta
      d_max_stop = dmax+params.pregxs.fitting.delta
      n_step = params.pregxs.fitting.n_step
      scanner =  d_max_scan(data,nparam,nfst,ntrials,d_max_start,d_max_stop,n_step,strials)
 


    else:
      fitters = random_start_fixed_dmax( data, dmax, nparam, nfst, ntrials,n_simplex=strials)
      coefs = fitters.trials[fitters.chi_index].solution
      for cc in coefs:
        print cc,
      print item, "COEF"
      pr_fit = fitters.trials[fitters.chi_index].get_best_pofr().f( data.q )
      print flex.mean(flex.pow2( (data.i-pr_fit)/(data.i+pr_fit) ) )*4.0, "CHI2"
      #fitters.show_calc_i(outfile=open(params.pregxs.output+'.qi_pr','w'))
      #fitters.show_obs_vs_calc(outfile=open(params.pregxs.output+'.pr_fit','w'))
      #print fitters.chi[fitters.chi_index],
      #fitters.print_entropy()
      #print params.pregxs.output


def get_input(args):
  if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil=master_params,
      home_scope="sas_I")

    print >> log, "#phil __OFF__"
    print >> log, "======================================================================"
    print >> log, "                              pregxs                                  "
    print >> log, "                P(r) Estimation Given Xray Scattering data            " 
    print >> log, "======================================================================"
    print >> log


    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log,expert_level=1)
    print >> log, "#phil __END__"

    go( params, log) 
    return True


def run(args):
  params = get_input(args)
  
def print_help():
  print "\nUsage:\nsastbx.Prfit data=target_pr"

  #run(filename,d, n)

if __name__ == "__main__":
  get_input( sys.argv[1:] )
#  print "OK"
