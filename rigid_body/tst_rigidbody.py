import os,sys,math
from scitbx.array_family import flex
from sastbx.rigid_body import rigidbody as rb
from sastbx.rigid_body import rb_engine as rbe
import time
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
import libtbx.phil
import libtbx.phil.command_line
from cStringIO import StringIO
from sastbx.basic_analysis import pofr_data
from scitbx import simplex
from sastbx.pr import pr_tools

master_params = iotbx.phil.parse("""\
pr{
  pdb_file= None
  .type=path
  .help="PDB filename to be calculated"

  pr_file= None
  .type=path
  .help=""

  Max = 200
  .type=float
  .help="upper limit of the histgram building"

  bin_size = 1.0
  .type=float
  .help="bin size for the histgram"

  output="modified.pdb"
  .type=path
  .help="Output file name for calculated P(r)"

}

""")

class fit(object):
  def __init__(self, body, expt_pr):
    self.body=body
    self.expt_pr = expt_pr
    self.dmax = int(self.expt_pr.r[-1]) + 1
    self.n = 2
    self.optimize()

  def optimize(self):
    start_matrix = []
    for ii in range( self.n + 1):
      start_matrix.append( (flex.random_double(self.n)-1.0)*0.2)
    optimizer = simplex.simplex_opt( dimension= self.n,
				     matrix = start_matrix,
				     evaluator = self,
				     tolerance=1e-4)
    self.solution = optimizer.get_solution()
    self.score = self.target( self.solution )

  def target(self,v):
    v = flex.abs(v + 1.0)
    self.body.set_weights(v[0],v[1])
    self.pr = self.body.get_hist()
    self.pr = self.pr[:self.dmax] / sum(self.pr)
    self.score = flex.sum(flex.pow2( self.pr - self.expt_pr.pr ))
    return self.score
    


def run(params,log):
  file=params.pr.pdb_file
  pr_file=params.pr.pr_file
  expt_pr = pofr_data.pofr(filename=pr_file, from_gnom=True)
  dMax = int(expt_pr.r[-1]) + 1
  r = flex.double(range(dMax) )
  new_pr = expt_pr.linear_interpolation( r )
  expt_pr = pofr_data.pofr( r = r, pr=new_pr)
  #dMax=params.pr.Max
  scale = 1.0
  pdb=rbe.PDB(file,scale)
  cutoff = 14
  Nmodes = 10
  maxi = 16
  grid = False
  q_array = flex.double( range( 101 )) / 200.0
  integrator=pr_tools.fast_integrator(dMax, q_array, div=1)

  indx = flex.int( range( pdb.xyz.size() ) )
  body=rb(pdb.xyz, indx, dMax, maxi, grid)
  refine = fit( body, expt_pr )
  solution = refine.solution + 1.0
  i_expt = integrator.get_intensity_from_pr_array( expt_pr.pr )
  i_fit  = integrator.get_intensity_from_pr_array( refine.pr )
  qmin, ratio = compare_i( i_expt, i_fit, q_array )
  print "#", dMax, qmin, 1.0, 
  for s in solution:
    print s,
  print math.sqrt(refine.score/dMax)

  for qq, ii, jj in zip( q_array, i_expt, i_fit):
    print qq, ii, jj



def compare_i( i1, i2, q, threshold=0.05):
  i1 = i1/i1[0]
  i2 = i2/i2[0]
  ratio = flex.abs(i1-i2)/(i1+i2)
  threshold = threshold /2.0

  for qq, rr in zip(q, ratio):
    if rr> threshold:
      return qq, rr



#  t=[0]*3
#  r=[0.2,0.5,0.5]
#  r=[0]*3
#  surface = body.get_surface_atoms()
#  for s in (surface):
#    print s
#  h = body.get_hist()
#  r = range(h.size())
#  sum = flex.sum( h )
#  h = h/sum
#  for ri, hi in zip(r,h):
#    if(hi > 0):
#      print ri, hi
  #crd=body.rotate_translate(t,r)
  #pdb.writePDB(crd, params.pr.output)


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
      home_scope="pr")

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
#    print >> log, "#phil __ON__"
#    new_params =  master_params.format(python_object=params)
#    new_params.show(out=log,expert_level=1)
#    print >> log, "#phil __END__"
    run( params, log )

def print_help():
  print "Usage:\n sastbx.pr pdb_file=PDB_FILE Max=MAX type=type(\"all\" or \"CA\") output=outputfile\n"

    
if __name__ == "__main__":
  get_input(sys.argv[1:])
