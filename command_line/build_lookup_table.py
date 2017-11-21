# LIBTBX_SET_DISPATCHER_NAME sastbx.build_lookup_table

from sastbx.data_reduction import build_lookup_file
import sys

if (__name__ == "__main__"):
    build_lookup_file.run(args=sys.argv[1:])
