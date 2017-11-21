from mmtbx import xmanip_tasks
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import ranges
import iotbx.phil
import iotbx.pdb
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx.math import matrix
from libtbx.utils import Sorry, multi_out
import libtbx.phil.command_line
from cStringIO import StringIO
import sys, os
import experimental_info


def print_help():
  print "No help yet"


#class essential_info(object):
  

#class single_image_integrator(object):
#  def 




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
    process_it = process_data( params )      



if __name__ == "__main__":
  get_input( sys.argv[1:] )

