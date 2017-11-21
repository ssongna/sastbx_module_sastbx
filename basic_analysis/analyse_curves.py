from sastbx.data_reduction import saxs_read_write
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import pdb, ranges, detectors
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from scitbx.math import scale_curves
import sys, os, math
from sastbx.data_reduction import saxs_read_write
import guinier_analyses, unified_fit, porod_analyses, kratky_tools


master_params =  iotbx.phil.parse("""\
analyze_curves{
  specimen
  .multiple=True
  {
    id = None
    .type = str
    .help = "Id of specimen"
    data
    .multiple=True
    {
      file_name = None
      .type=path
      .help="A data set"

      concentration = None
      .type=float
      .help="Concentration in mg/ml"
    }
  }
}

""")



def uniform_limit(q,i,i0=None,log_rat=-2.0):
  if i0 is None:
    i0 = i[0]
  for qq,ii in zip(q,i):
    if ii < i0*pow(10.0, log_rat ):
      return qq
  return qq




class simple_specimen(object):
  def __init__(self, data_def, out=None):
    self.out=out
    if self.out is None:
      self.out = None

    self.data_def = data_def
    self.data_id = self.data_def.id
    self.file_names = []
    self.concentration = []
    self.ssd = []
    for data in self.data_def.data:
      file_name = data.file_name
      conc = data.concentration
      tmp_sd = saxs_read_write.read_standard_ascii_qis( file_name )
      tmp_sd.show_summary()
      self.ssd.append( tmp_sd )
      self.file_names.append( file_name )
      self.concentration.append( conc )
      msga = guinier_analyses.multi_step_rg_engine( tmp_sd , out=self.out)
      kratky_tools.kratky_analyses(tmp_sd)



  def scale_data(self, ref_id=0):
    m = []
    v = []
    for ddd in self.ssd:
      m.append( ddd.s )
      v.append( ddd.s*ddd.s  )
    scales, offsets = scale_curves.scale_it_pairwise( m,v,ref_id,factor=10,show_progress=False,add=False)
    for scale, ssd in zip(scales,self.ssd):
      self.ssd.multiply( scale )







def run(params, log):
  for spec in params.analyze_curves.specimen:
    spec_spec = simple_specimen( spec, out=log )




def print_help():
  print "No help yet"


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
      home_scope="analyze_curves")

    print >> log, "#phil __OFF__"
    print >> log, "=========================="
    print >> log, "            AC            "
    print >> log, "      Analyze Curves      "
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

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log,expert_level=1)
    print >> log, "#phil __END__"
    run( params, log )

if __name__ == "__main__":
  get_input( sys.argv[1:] )
