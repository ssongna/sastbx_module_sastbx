# LIBTBX_SET_DISPATCHER_NAME sastbx.build_db

from sastbx.zernike_model import build_db
import sys

if (__name__ == "__main__"):
  build_db.run( sys.argv[1:] )

