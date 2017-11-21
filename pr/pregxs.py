import pregxs_tools as pretls
from sastbx.intensity.sas_I import reduce_raw_data
from sastbx.intensity.sas_I import write_json
import sys, os, math
import iotbx.phil
import random,time
from sastbx.data_reduction import saxs_read_write
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from sastbx.basic_analysis import guinier_analyses
from libtbx.utils import Sorry, date_and_time, multi_out
from sastbx.interface import get_input
from scitbx.array_family import flex


base_path = os.path.split(sys.path[0])[0]

global targetfile
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
    n_trials=4
    .type=int
    .help="n_trials fittings will be performed, only the best will be taken."
    n_trials_simplex=5
    .type=int
    .help="n_trials_simplex random restarts of simplex will be performed. make sure this aint 1"

    rg = 0
    .type=float
    .help="Radius of Gyration"

    delta = 15
    .type=float
    .help="scan range:  dmax-delta - > dmax+delta"

    step = 3
    .type=int
    .help="Number of steps in dmax scan"

  }
  output = "pregxs"
  .type=path
  .help = "Output base. expect a .pr and a .qii file."

}
""")

banner = "-------------------P(r) Estimation using parametric function-------------------"

def write_pr_json(filename, x,y1):

  out=open(filename, 'w')

  head = '{"elements":['
  ele1 = '{"type": "line", "dot-style":{"type": "dot", "dot-size": 3, "colour":"#DFC329"}, "width":3, "colour": "DFC329", "text":"model", "font-size":10,'
  #string1 = out+head+ele1+'"values":\n'
  # with open(targetfile,"w") as f:
  #   f.write(str(out)+str(head)+str(ele1)+'"values":')

  print>>out, head, ele1, '"values":',
 
  y_min=0
  y_max=flex.max(y1)

  # with open(targetfile,"w") as f:
  #   f.write(str(out)+'['+ ', '.join('%5.6f' % v for v in y1) + ']')
  #   f.write(str(out)+'}'+']\n')
  print>>out, '[' + ', '.join('%5.6f' % v for v in y1) + ']',
  print>>out, '}',  # end of y1

  print>>out, ']',  # end of elements
  ### now setup axis ###
  # with open(targetfile,"w") as f:
  #   f.write(str(out)+',"y_axis":{"min":%f, "max":%f}\n'%(y_min, y_max))
  print>>out,',"y_axis":{"min":%f, "max":%f}'%(y_min, y_max),

  steps = 10
  x_labels = '["'
  for xx in x:
    x_labels = x_labels + str(xx) + '","'
  x_labels=x_labels[0:-2]+']'
  # with open(targetfile,"w") as f:
  #   f.write(str(out)+',"x_axis":{"min":%d, "max":%d, "steps":%d, "labels":{"labels":%s,"steps":%d}}'%(0,x.size(),steps,x_labels, steps))
  #   f.write(str(out)+'}\n')
  print>>out,',"x_axis":{"min":%d, "max":%d, "steps":%d, "labels":{"labels":%s,"steps":%d}}'%(0,x.size(),steps,x_labels, steps),
  print>>out,'}' ##end of the file
  out.close()


def go(params, log=None):
  if log is None:
    log = sys.stdout

  log2 = sys.stdout
  data_array = []
  multies = []
  d_max=int(params.pregxs.d_max+0.5)
  rg=params.pregxs.fitting.rg
  increment = math.pi / d_max
  n_params=params.pregxs.fitting.n_coeff
  n_fst_pass=params.pregxs.fitting.n_fst_pass
  n_trial=params.pregxs.fitting.n_trials
  n_simplex=params.pregxs.fitting.n_trials_simplex
  prior = None

  targetfile = os.path.join(os.path.split(sys.path[0])[0],"pregxs.txt") 
  for item in params.pregxs.data:
    data = saxs_read_write.read_standard_ascii_qis(item)
    qmax = data.q[-1]
    bandwidth = 0.01  # can be changed is needed
    data = reduce_raw_data( data, qmax, bandwidth, level=0.0001, outfile=targetfile)
    m = 1.0/data.i[0]
    data.multiply_add(m,0.0)
    data_array.append( data )
    multies.append( m )

    q_min = data.q[0]
    q_max0 = data.q[ data.q.size() - 1 ]
    nparams_max = int ((q_max0 - q_min) / increment)
    # print "n_params: ",n_params
    # print "nparams_max: ", nparams_max
    if(nparams_max < n_params):
      with open (targetfile,"a") as f:
        f.write("WARNING: number of parameters is larger than maximum number of shannon channels covered by expt data\n")
      print  "WARNING: number of parameters is larger than maximum number of shannon channels covered by expt data"

    if params.pregxs.scan:
      delta = int(params.pregxs.fitting.delta)
      step = params.pregxs.fitting.step
      scanner  = pretls.dmax_scan(prior, data, d_max, delta, step, rg, n_params, n_fst_pass, n_trial=n_trial, n_simplex=n_simplex,entropy_thresh=1.24, outfile=targetfile)
      scanner.print_pr( open(params.pregxs.output+"average.pr", 'w') )
      scanner.get_best_dmax()
      fitter = scanner.fitter
      fitter.best_fit.show_pr( open(params.pregxs.output+"best.pr",'w') )
      fitter.best_fit.show_obs_vs_calc( open(params.pregxs.output+"best.qii",'w') )
      write_pr_json(params.pregxs.output+"data.json", scanner.r, scanner.average_pr)
      write_json( params.pregxs.output+"qii.json", data.q, scanner.calc_i, data.i )

    else:
      fitter = pretls.fixed_dmax_fitter(prior, data, d_max, n_params, n_fst_pass, n_trial=n_trial, n_simplex=n_simplex)
      fitter.best_fit.show_pr( open(params.pregxs.output+"best.pr",'w') )
      fitter.best_fit.show_obs_vs_calc( open(params.pregxs.output+"best.qii",'w') )
      write_pr_json(params.pregxs.output+"data.json", fitter.r, fitter.pr)
      write_json( params.pregxs.output+"qii.json", data.q, fitter.calc_i, data.i)

    

def run(args):
  t1 = time.time()
  targetfile = os.path.join(os.path.split(sys.path[0])[0],"pregxs.txt")
  with open(targetfile,"w") as f:
    f.truncate()
    
  tempf = open(targetfile,'a')
  params = get_input(args, master_params, "pregxs", banner, print_help, tempf)
  tempf.close()
  if params is None: return
  go( params)
  with open(targetfile,"a") as f:
    f.write(str(time.time() - t1)+" time use\n")
    f.write("__END__")

  return

def print_help(out=sys.stdout):
  print "\nUsage:\n  sastbx.pregxs data=data_file d_max=dmax scan=True/False output=output_prefix\n"

if __name__ == "__main__":
  run( sys.argv[1:] )
