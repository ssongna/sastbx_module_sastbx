import libtbx.load_env
import os
Import("env_etc")

env_etc.sastbx_dist = libtbx.env.dist_path("sastbx")
env_etc.sastbx_include = os.path.dirname(env_etc.sastbx_dist)

env_etc.sastbx_common_includes = [
  env_etc.libtbx_include,
  env_etc.sastbx_include,
  env_etc.cctbx_include,
  env_etc.scitbx_include,
  env_etc.boost_include,
]

if (not env_etc.no_boost_python):
  Import("env_scitbx_boost_python_ext")
  env_bpl = env_scitbx_boost_python_ext.Copy()
  env_etc.include_registry.append(
    env=env_bpl,
    paths=[env_etc.sastbx_include] + env_etc.scitbx_common_includes)

SConscript("data_reduction/SConscript")
SConscript("intensity/SConscript")
SConscript("refine/SConscript")
SConscript("rigid_body/SConscript")
SConscript("zernike_model/SConscript")
SConscript("pr/SConscript")
SConscript("fXS/SConscript")
SConscript("md/SConscript")
