# LIBTBX_SET_DISPATCHER_NAME sastbx.buildmap

from sastbx.zernike_model import model_interface
import sys

if (__name__ == "__main__"):
  model_interface.run( sys.argv[1:] )

