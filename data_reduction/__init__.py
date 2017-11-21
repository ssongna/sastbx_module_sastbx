import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("sastbx_data_reduction_ext")
from sastbx_data_reduction_ext import *

# make objects picklable
# perpendicular detector
# getinitargs
def perpendicular_detector_getinitargs(self):
  return ( self.nx(), self.ny(), self.scale_x(), self.scale_y(), self.origin_x(), self.origin_y(), self.distance(), self.wavelength() )
perpendicular_detector.__getinitargs__ = perpendicular_detector_getinitargs
