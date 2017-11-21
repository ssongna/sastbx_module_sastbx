# LIBTBX_SET_DISPATCHER_NAME sastbx.image_simulator

from sastbx.fXS import image_simulator
import sys

if (__name__ == "__main__"):
  image_simulator.run( sys.argv[1:] )
