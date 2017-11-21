from iotbx import ranges
import iotbx.phil
import iotbx.pdb
from cctbx.array_family import flex
from scitbx.math import matrix
from libtbx.utils import Sorry, multi_out
import libtbx.phil.command_line
from cStringIO import StringIO
import sys, os
import experimental_info, integrator, average, saxs_read_write
from iotbx import ranges
import pickle


def print_help():
  print "No help yet"

class specimen(object):
  def __init__(self, specimen_params):
    self.params = experimental_info.specimen_specs(specimen_params)
    self.integrated_data = []

  def image_names(self):
    return self.params.file_names

  def add_data(self, integrated_data):
    self.integrated_data.append( integrated_data )

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "----   Specimen information   ----"
    self.params.show()

def process_specimen(params):
  specimens = []
  for specimen_spec in params.specimen:
    tmp = specimen( specimen_spec )
    specimens.append( tmp )
  return specimens

class data_set(object):
  def __init__(self,params):
    self.specimen_array = process_specimen( params )
    self.hutch = experimental_info.hutch_specs( params.hutch_setup )
    sample_image = self.specimen_array[0].image_names()[0]
    self.hutch.set_values( sample_image )
    self.q_range = None
    self.simple_saxs_data_sets = []
    self.transmission_factors = []
    self.this_buffer, self.this_window = self.buffer_and_window_id()
    self.corrected_simple_data_sets = []

  def buffer_and_window_id(self):
    buffer_count = 0
    window_count = 0
    buffer_id = None
    window_id = None
    for spec, ii in zip(self.specimen_array, xrange(len(self.specimen_array)) ):
      if spec.params.type=="buffer":
        buffer_count += 1
        buffer_id = ii
      if spec.params.type=="window":
        window_count += 1 
        window_id = ii

    #assert buffer_count == 1
    #assert window_count == 1
    return buffer_id, window_id

  def simple_background_subtract(self, return_sets=True):
    n = len(self.simple_saxs_data_sets)
    print  self.this_buffer, self.this_window
    for ssds, ii in zip(self.simple_saxs_data_sets,xrange(n)):
       
      if not ( (ii == self.this_buffer) or (ii == self.this_window) ) :
        new_data = ssds.deep_copy()
        print flex.mean( new_data.i ), 
        new_buffer = self.simple_saxs_data_sets[self.this_buffer].deep_copy()
        print flex.mean( new_buffer.i )
        for q,ic,im in zip( new_data.q, new_data.i, new_buffer.i):
          print q,ic,im,ic-im

        new_data.add_curve( new_buffer )
        new_data.data_id = new_data.data_id+"_bckgrndsbtrctd"
        self.corrected_simple_data_sets.append( new_data )
        new_data.show_summary()
    if return_sets:
      return self.corrected_simple_data_sets

  def integrate_data(self, processing_options, out=None):
    if out is None:
      out = sys.stdout

    integration_object = integrator.integrate( self.hutch, processing_options.integration )

    q_low, q_mean, q_step = integration_object.get_q_arrays()
    self.q_range = integrator.q_range( q_low, q_mean, q_step )
    for this_specimen in self.specimen_array:
      for this_image_name in this_specimen.image_names():
        print >> out, "Integrating %s "%(this_image_name)
        mi,vi = integration_object.integrate_image( this_image_name  )
        new_radially_integrated_dataset = integrator.radially_integrated_data( self.q_range )
        new_radially_integrated_dataset.load_data(mi,vi)
        transmission, i_tot,i_after = integration_object.get_aux_data(this_image_name)
        new_radially_integrated_dataset.set_aux_info(transmission, i_tot, i_after)
        tmp_image_name = this_image_name.split("/")
        tmp_image_name = tmp_image_name[ len(tmp_image_name)-1 ]
        #if  
        new_radially_integrated_dataset.write_data( file_name = tmp_image_name+".raw_data.qis"  ) 
        this_specimen.add_data( new_radially_integrated_dataset )
        print >> out, "    ---->  I1/I0: %8.5f    I0: %5.3e   I1: %5.3e     "%(transmission, i_tot, i_after)
        print >> out

  def return_specimen(self):
    """This should be all that is needed for averaging curves"""
    return self.specimen_array

  def analyze_and_merge_specimens(self, scaling_options):
    #first we need a reference scale
    # lets make that equal to I0 for the first image of the first specimen encountered
    ## ref_scale = self.specimen_array[0].
    ref_scale = 1.0 #self.specimen_array[0].integrated_data[0].i_tot

    for spec in self.specimen_array:
      avg = average.averaging_engine(spec, scaling_options)
      # this is done in scale to reference
      avg.scale_to_reference(0)
      #show me the data
      avg.show()
      #write average intensities
      avg.write_average_intensity()
      avg_data = avg.get_average_simple_saxs_data()
      self.simple_saxs_data_sets.append( avg_data )
      self.transmission_factors.append( (avg.mean_trans,avg.std_trans) )
    return self.simple_saxs_data_sets

def process_data_sets(params):
  data_sets = []
  for ds in params.experiment.data_set:
    data_sets.append( data_set( ds ) )
  return data_sets

def show_data_sets(data_sets,out=None):
    if out is None:
      out = sys.stdout
    for ds in data_sets:
      print "########### DATA SET ###########"
      for spec in ds.specimen_array:
        spec.show(out)
      ds.hutch.show(out)

    print "################################"


def run(params,out):
  data_sets = process_data_sets( params )
  show_data_sets( data_sets )
  integrated_data_sets = []
  simple_saxs_curves = []
  corrected_simple_saxs_curves = []
  print dir( params.processing_options )
  for data_set in data_sets:
    data_set.integrate_data( params.processing_options )
    integrated_data_sets.append( data_set.return_specimen() )
    ssds = data_set.analyze_and_merge_specimens( params.processing_options.scaling.averaging )
    cssds = data_set.simple_background_subtract()
    saxs_read_write.write_array_of_saxs_data( cssds )
    simple_saxs_curves.append( ssds )
    corrected_simple_saxs_curves.append( cssds )

  for ssds in simple_saxs_curves:
    for ssd in ssds:
      ssd.show_summary()



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
      master_phil=experimental_info.master_params,
      home_scope="experiment")

    print >> log, "#phil __OFF__"
    print >> log, "=========================="
    print >> log, "           SDP            "
    print >> log, "   SAXS Data Processing   "
    print >> log, "=========================="
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

    effective_params = experimental_info.master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    print >> log, "#phil __ON__"
    new_params =  experimental_info.master_params.format(python_object=params)
    new_params.show(out=log,expert_level=1)
    print >> log, "#phil __END__"
    run( params, log )


if __name__ == "__main__":
  get_input( sys.argv[1:] )
