# LIBTBX_SET_DISPATCHER_NAME sastbx.sdp

from sastbx.data_reduction import sdp 
import sys

if (__name__ == "__main__"):
    sdp.get_input( args=sys.argv[1:] )

