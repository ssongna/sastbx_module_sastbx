# LIBTBX_SET_DISPATCHER_NAME sastbx.refine_rb

from sastbx.rigid_body import rigid_body_local_refine
import sys

if (__name__ == "__main__"):
  rigid_body_local_refine.get_input( sys.argv[1:] )

