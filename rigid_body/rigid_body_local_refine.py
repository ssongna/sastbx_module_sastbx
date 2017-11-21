import os, sys
from scitbx.array_family import flex
from scitbx import direct_search_simulated_annealing
from math import pi, exp
import random, time
from scitbx import lbfgs, simplex
from sastbx.data_reduction import saxs_read_write
from sastbx.rigid_body import rigidbody as rb
from sastbx.rigid_body import rb_engine as rbe
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
import libtbx.phil
import libtbx.phil.command_line
from cStringIO import StringIO
from iotbx import pdb

master_params = iotbx.phil.parse("""\
refine{
  data
  {
    type = Pr *Iq
    .type = choice
    .help = "the target data type: I(q) or P(r), default: I(q) "

    file = None
    .type = path
    .help = "File with data"
  }


  model
  .multiple=True
  {
    file=None
    .type=path
    .help="File name with model"
    rigid_body
    .multiple=True
    {
     selection=None
     .type=str
     .help="rigid_body definition as atom selection"

     fixed_orientation=False
     .type=bool
     .help="Fixed orientation?"
  
     fixed_position=False
     .type=bool
     .help="Fixed position?"

     max_rotation = 20.0
     .type=float
     .help="maximum orientation shift"

     max_translation=10.0
     .type=float
     .help="Maximum translation shift"
    }
  }

  settings{
    bin_size = 1
    .type=float
    .help="bin size for the histgram"

    sub_sample_ratio=1
    .type = int
    .help = "Determines sub sample ratio of histogram for p(r) calculation (group size)" 

  }
  
  output="refined"
  .type=path
  .help="prefix of the output file names"

}

""")




class refine_rb(object):
  def __init__(self,rbe_obj,expt_data, data_type='pr', shift=None, both=True, simplex_trial=10, n_step=20):
    self.data = expt_data
    self.rbe=rbe_obj
    self.nbody=rbe_obj.nbody
    self.n_refine=self.nbody - 1

    self.count=0
    self.simplex_trial = simplex_trial
    self.data_type = data_type

    self.translate = True
    self.both = both

    if(self.both):
      self.n = 6
    else:
      self.n = 3

    if(shift is None):
      self.x = flex.double(self.n*self.n_refine, 0)
    else:
      self.x =flex.double()
      for ii in range(1, self.nbody):
        self.x = self.x.concatenate( shift[ii] )
        if(self.both):
          self.x = self.x.concatenate( flex.double(3,0) )

    self.angle = flex.double(3,0)
    self.vector = flex.double(3,0)
    self.scores = flex.double()

    self.simplex()

  def run_dssa(self):
    self.starting_simplex = []
    for ii in range(self.n):
      self.starting_simplex.append((flex.random_double(self.n*self.n_refine)*2-1.0)*1.0+self.x)
    self.starting_simplex.append(self.x)

    dssa_optimizer = direct_search_simulated_annealing.dssa(     dimension=self.n,
                                    matrix  = self.starting_simplex,
                                    evaluator = self,
                                    max_iter=5000,
                                    further_opt=True,
                                    tolerance=1e-4)
    candidates = dssa_optimizer.get_candi()

    self.solution = dssa_optimizer.get_solution()
    self.score = self.target(self.solution)


  def simplex(self):
    self.simplex_scores = []
    self.simplex_solutions = []
    for ii in xrange(self.simplex_trial):
      #make a random simplex
      self.starting_simplex = []
      for ii in range(self.n*self.n_refine):
        self.starting_simplex.append((flex.random_double(self.n*self.n_refine)*2-1.0)*1.0+self.x)
      self.starting_simplex.append(self.x)

      self.optimizer = simplex.simplex_opt( dimension=self.n*self.n_refine,
                                    matrix  = self.starting_simplex,
                                    evaluator = self,
                                    max_iter=500,
                                    tolerance=1e-4)

      self.solution = self.optimizer.get_solution()
      self.score = self.target( self.solution )
      self.simplex_scores.append( self.score )
      self.simplex_solutions.append( self.solution )

    best_index = flex.min_index( flex.double(self.simplex_scores) )
    self.solution = self.simplex_solutions[ best_index ]
    self.score = self.target(self.solution)



  def target(self, combined_vector):
    vector = []
    for ii in range(self.n_refine):
      vector.append( combined_vector[ii*self.n:(ii+1)*self.n] )

    if(self.both):
      for ii in range(self.n_refine):
        self.rbe.rotate_translate( vector[ii][0:3], vector[ii][3:], ii+1 )
    elif(self.translate):
      for ii in range(self.n_refine):
        self.rbe.rotate_translate( vector[ii], self.angle, ii+1)
    else:
      for ii in range(self.n_refine):
        self.rbe.rotate_translate( self.vector, vector[ii], ii+1 )


    if(self.data_type == 'pr'):
      calc_pr = self.rbe.get_norm_pr()
      score = 1000*flex.sum( flex.abs(calc_pr-self.data)*flex.exp(-self.data*self.data) )
      clash = self.rbe.get_clash()
      score += clash/score
    else:
      calc_i = self.rbe.get_intensity()
      scale = flex.sum( calc_i*self.data.i ) / flex.sum( flex.pow2(calc_i) )
      score = flex.sum( flex.abs( scale*calc_i - self.data.i ) )
    return score


#=================================================================================================#
#
#=================================================================================================#

def run_multiple_pdb(params, log):
  dmax=params.refine.dmax
  group_size = params.refine.group_size
  max_num_fibonacci = 12

  target_data = saxs_read_write.read_standard_ascii_qis(params.refine.target)
  if(params.refine.data_type == 'pr'):
    dmax = int(target_data.q[-1]+0.5)
    new_r = flex.double( range(dmax) )
    target_data = flex.linear_interpolation( target_data.q, target_data.i, new_r )

  rbs=[]
  center = []
  pdb_objects = []

  main_body = True

  max_size = 0
  for item in params.refine.model:
    pdb_obj = pdb.hierarchy.input(item)
    size = pdb_obj.hierarchy.atoms_size()
    atom_indx = flex.int( range(0, size, group_size ) )
    xyz = flex.vec3_double()
    atoms = pdb_obj.hierarchy.atoms()
    for a in atoms:
      xyz.append( a.xyz )
    if(size > max_size):
      max_size = size
      rbs.insert(0, rb(xyz, atom_indx, max_num_fibonacci, shift, rotate))
      pdb_objects.insert(0, pdb_obj )
    else:
      rbs.append(rb(xyz, atom_indx, max_num_fibonacci, shift, rotate))
      pdb_objects.append( pdb_obj )

  num_body = len( pdb_objects )


  shift = [(0,0,0)]
  if(params.refine.data_type == 'pr'):
    rb_eng = rbe.rb_engine(rbs,int(dmax) )
  else:
    rb_eng = rbe.rb_engine(rbs,int(dmax), q_array=target_data.q )

  for ii in range(1, num_body):
    shift.append( flex.double( rbs[ii].center() ) - flex.double( rbs[0].center() ) )
  #  rb_eng.rbs[ii].rotate_only(list( flex.random_double(3,) ), 10.0/180.0*pi)
    rb_eng.rbs[ii].translate_after_rotation(list(shift[ii]) )

  filename="initial_model.pdb"
  write_pdb(filename, num_body, rb_eng, pdb_objects)

  refine = refine_rb( rb_eng, target_data, data_type=params.refine.data_type, shift=shift, both=True)

  solution = refine.solution
  refine.target( solution )

  filename=params.refine.output+".pdb"
  write_pdb(filename, num_body, rb_eng, pdb_objects)


def run_single_pdb(params, log):
  group_size = params.refine.group_size
  max_num_fibonacci = 12

  target_data = saxs_read_write.read_standard_ascii_qis(params.refine.target)
  if(params.refine.data_type == 'pr'):
    dmax = int(target_data.q[-1]+0.5)
    new_r = flex.double( range(dmax) )
    target_data = flex.linear_interpolation( target_data.q, target_data.i, new_r )

  rbs=[]
  center = []
  pdb_objects = []

  main_body = True

  pdb_inp = pdb.hierarchy.input(params.refine.model[0])
  cache = pdb_inp.hierarchy.atom_selection_cache()
  max_size = 0
  for item in params.refine.rigid_body:
    fix_location = item.fixed_position
    fix_orientation = item.fixed_orientation
    cache_selected = cache.selection(string=item.selection)
    pdb_obj = pdb_inp.hierarchy.atoms().select(cache_selected)
    size = pdb_obj.size()
    atom_indx = flex.int( range(0, size, group_size ) )
    xyz = flex.vec3_double()
    atoms = pdb_obj
    for a in atoms:
      xyz.append( a.xyz )
    if(size > max_size):
      max_size = size
      rbs.insert(0, rb(xyz, atom_indx, dmax, max_num_fibonacci, fix_location, fix_orientation))
      pdb_objects.insert(0, pdb_obj )
    else:
      rbs.append(rb(xyz, atom_indx, dmax, max_num_fibonacci, fix_location, fix_orientation))
      pdb_objects.append( pdb_obj )

  num_body = len( pdb_objects )


  shift = [(0,0,0)]
  if(params.refine.data_type == 'pr'):
    rb_eng = rbe.rb_engine(rbs,int(dmax) )
  else:
    rb_eng = rbe.rb_engine(rbs,int(dmax), q_array=target_data.q )

  for ii in range(1, num_body):
    shift.append( flex.double( rbs[ii].center() ) - flex.double( rbs[0].center() ) )
  #  rb_eng.rbs[ii].rotate_only(list( flex.random_double(3,) ), 10.0/180.0*pi)
    rb_eng.rbs[ii].translate_after_rotation(list(shift[ii]) )

  filename="initial_model.pdb"
  write_pdb_single(filename, num_body, rb_eng, pdb_inp, pdb_objects)

  refine = refine_rb( rb_eng, target_data, data_type=params.refine.data_type, shift=shift, both=True)

  solution = refine.solution
  refine.target( solution )

  filename=params.refine.output+".pdb"
  write_pdb_single(filename, num_body, rb_eng, pdb_inp, pdb_objects)

#===============================================================================

def write_pdb(filename, num_body, rb_eng, pdb_objects):
  for ii in range(num_body):
    new_xyz = rb_eng.rbs[ii].get_crd()
    for a, xyz in zip( pdb_objects[ii].hierarchy.atoms(), new_xyz ):
      a.set_xyz( new_xyz=xyz )
    if(ii==0):
      pdb_objects[ii].hierarchy.write_pdb_file( file_name=filename, open_append=False)
    else:
      pdb_objects[ii].hierarchy.write_pdb_file( file_name=filename, open_append=True)


def write_pdb_single(filename, num_body, rb_eng, pdb_inp, pdb_objects):
  for ii in range(num_body):
    new_xyz = rb_eng.rbs[ii].get_crd()
    for a, xyz in zip( pdb_objects[ii], new_xyz ):
      a.set_xyz( new_xyz=xyz )
  pdb_inp.hierarchy.write_pdb_file( file_name=filename, open_append=False)


#===============================================================================
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
      home_scope="refine")

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
    if( len(params.refine.model) > 1):








#===============================================================================
      run_multiple_pdb( params, log )
    elif( len(params.refine.model) == 1):
      run_single_pdb( params, log )

def print_help():
  print "Usage:\n "
  print "sastbx.refine_rb model=single_pdb_file rigid_body={rb1, rb2} target=target.dat data_type=*pr or i"
  print "or:"
  print "sastbx.refine_rb model={pdb_file1,pdb_file2} target=target.dat data_type=*pr or i"
  print



#==================================================================



#==============================================================================================#
#
#==============================================================================================#
if __name__ == "__main__":
  t1 = time.time()
  get_input(sys.argv[1:])
  print "Total Time Used: ", time.time() - t1
