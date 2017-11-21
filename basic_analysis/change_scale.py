from sastbx.data_reduction import saxs_read_write
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import pdb, ranges, detectors
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from libtbx import easy_pickle
from scitbx.math import scale_curves
import sys, os, math
from sastbx.data_reduction import saxs_read_write
import guinier_analyses, unified_fit, porod_analyses, kratky_tools


master_params =  iotbx.phil.parse("""\
change_scale{
  input{
    data=None
    .type=path
    .help="Filename of data"
    scale = 4pistol_A 4pistol_nm 2stol_A 2stol_nm
    .type=choice
    .help="Angular scale of input"
  }

  output{
    postfix="rescaled.qis"
    .type=str
    .help="Append this to the filename given earlier"
    scale=*4pistol_A 4pistol_nm 2stol_A 2stol_nm
    .type=choice
    .help="Output angular scale"
  }
}

""")


angular_scales = { "4pistol_A"  : 1.0,
                   "4pistol_nm" : 10.0,
                   "2stol_A"    : 1.0/6.28318531,
                   "2stol_nm"   : 10.0/6.28318531    }

angular_names = { "4pistol_A"   : "4*Pi*Sin(theta)/Lambda (A^-1)",
                   "4pistol_nm" : "4*Pi*Sin(theta)/Lambda (nm^-1)",
                   "2stol_A"    : "2*Sin(theta)/Lambda (A^-1)",
                   "2stol_nm"   : "2*Sin(theta)/Lambda (nm^-1)"    }


def run(params, log):
  print >>log, "Changing from input scale :      %s"%(angular_names[params.change_scale.input.scale])
  print >>log, "to output scale           :      %s"%(angular_names[params.change_scale.output.scale])
  # read in the SAXS data
  mydata = saxs_read_write.read_standard_ascii_qis( params.change_scale.input.data )
  #first change the scale to 2stol_A
  mydata.q = mydata.q/angular_scales[params.change_scale.input.scale]
  #now go to where we have to go
  mydata.q = mydata.q*angular_scales[params.change_scale.output.scale]
  new_file_name = params.change_scale.input.data+"."+params.change_scale.output.postfix
  print >>log, "Writing output file: ",new_file_name
  saxs_read_write.write_standard_ascii_qis(mydata, new_file_name)


def print_help():
  print "Rescale SAXS curves"


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
    print >> log, "      sastbx.rescale      "
    print >> log, "   Rescale angular values "
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
