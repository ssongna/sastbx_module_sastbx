# LIBTBX_SET_DISPATCHER_NAME sastbx.shapeup

from sastbx.zernike_model import search_pdb
import sys

if (__name__ == "__main__"):
  search_pdb.run( sys.argv[1:] )

