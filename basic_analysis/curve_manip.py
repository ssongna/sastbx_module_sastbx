from scitbx.array_family import flex
from scitbx.math import scale_curves
from sastbx.data_reduction import curves



class curve_interpolator(object):
  def __init__(self, data_sets):
    self.data_sets = data_sets
    self.min_q = None 
    self.max_q = None
    for ds self.data_sets:
       




