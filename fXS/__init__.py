import boost.python
ext = boost.python.import_ext("sastbx_fXS_ext")
from sastbx_fXS_ext import *

import libtbx.load_env
if (libtbx.env.build_options.enable_cuda):
  cuda_ext = boost.python.import_ext("sastbx_fXS_cuda_ext")
  from sastbx_fXS_cuda_ext import *
