# LIBTBX_SET_DISPATCHER_NAME sastbx.superpose

from sastbx.zernike_model import zalign
import sys

if (__name__ == "__main__"):
  zalign.run( sys.argv[1:] )

