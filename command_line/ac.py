# LIBTBX_SET_DISPATCHER_NAME sastbx.ac

from sastbx.basic_analysis import analyse_curves
import sys 

if (__name__ == "__main__"):
    analyse_curves.get_input( args=sys.argv[1:] )

