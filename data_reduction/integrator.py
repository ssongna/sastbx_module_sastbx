from iotbx import ranges
import iotbx.phil
import iotbx.pdb
from cctbx.array_family import flex
from scitbx.math import matrix
from libtbx.utils import Sorry, multi_out
import libtbx.phil.command_line
from cStringIO import StringIO
import sys, os
import experimental_info
import sastbx.data_reduction
from sastbx.data_reduction import curves, saxs_read_write
from iotbx import detectors
import pickle

non_pedestal_image_types = ["EDF"]

def image_as_flex_int(image_name,image_type=None,subtract_pedestal=True):
  img_factory = detectors.TrySingleImageType(image_name,image_type)
  img_factory.readHeader(external_keys=[('PEDESTAL','PEDESTAL',int)])
  img_factory.read()
  result = img_factory.linearintdata
  result = flex.int( list(result) )
  if subtract_pedestal:
    if image_type not in non_pedestal_image_types:
      result = result - img_factory.parameters['PEDESTAL']

  return result

def get_header(image_name, image_type, external_keys):
  img_factory = detectors.TrySingleImageType(image_name,image_type)
  img_factory.readHeader(external_keys=external_keys)
  return img_factory.parameters

def get_values_from_header(image_name, header_keywords, image_type, external_keys=None ):
  img_factory = detectors.TrySingleImageType(image_name,image_type)
  img_factory.readHeader(external_keys=external_keys)
  result = []
  for key in header_keywords:
    result.append( img_factory.parameters[key] )
  return result






def parse_ion_chamber_keys( txt ):
  result = []
  keys = txt.split(",")
  for key in keys:
    tmp = ""
    for letter in key:
      if letter is not " ":
        tmp+=letter
    result.append( tmp  )
  return result

def get_total_ion_chamber_value( image_name, image_type, keys):
  # first we check if the header has these values
  ion_chamber_keys = [ ('I0_1', 'I0_1', float),
                       ('I0_2', 'I0_2', float),
                       ('I1_1', 'I1_1', float),
                       ('I1_2', 'I1_2', float) ]


  header = get_header( image_name, image_type, ion_chamber_keys )
  all_there = True
  for key in keys:
    if not header.has_key( key):
      all_there = False
  if all_there:
    result = 0
    for key in keys:
      result += header[ key ]
    return result
  else:
    raise Sorry( "Not all keys from %s can be found"%str(keys) )





class integrate(object):
  def __init__(self, hutch, integration_options):
    # input
    self.hutch = hutch
    self.integration_options = integration_options

    self.aux_info_object = None

    # make a detector
    self.detector = sastbx.data_reduction.perpendicular_detector(self.hutch.nx,self.hutch.ny,
                             self.hutch.scale_x, self.hutch.scale_y,
                             self.hutch.beam_origin_x, self.hutch.beam_origin_y,
                             self.hutch.distance, self.hutch.wavelength )
    # build a mask
    self.mask_object = sastbx.data_reduction.image_mask(self.hutch.nx,self.hutch.ny)
    self.build_mask()
    self.mask_data = self.mask_object.mask()

    # build the shadow area
    self.shadow_object = sastbx.data_reduction.image_mask(self.hutch.nx,self.hutch.ny)
    self.build_shadow()


    # get some bin info
    self.min_q = self.integration_options.q_range.min_q
    self.max_q = self.integration_options.q_range.max_q
    self.q_step = self.integration_options.q_range.q_step

    self.image_q_range = self.detector.q_range()

    if self.min_q is None:
      self.min_q = int(self.image_q_range[0]/self.q_step)*self.q_step
    if self.max_q is None:
      self.max_q = self.image_q_range[1]
    if self.q_step is None:
      self.q_step = 0.001

    # precompute the bins
    self.detector.compute_bin_indices( self.q_step,
                                       self.min_q,
                                       self.max_q )
    #get q bin values
    self.q_low  = self.detector.q_low_bin_values()
    self.q_mean = self.detector.q_mean_bin_values()
    self.q_step = self.integration_options.q_range.q_step

    # make an integrator object
    self.integrator_engine = sastbx.data_reduction.radial_integrator( self.detector, self.mask_data, None )

    


  def build_mask(self):
    # XXX WE NEED TO DO SOMETHING ABOUT 'INSIDE/OUTSIDE' DEFINITIONS XXX #
    # first add polygons please
    for polyg in self.hutch.exclude_polygons:
      tmp = sastbx.data_reduction.image_mask(self.hutch.nx, self.hutch.ny)
      tmp.add_polygon( flex.int(polyg.corners_x), flex.int(polyg.corners_y), 0,0 ) # we need a solution for this please
      self.mask_object.add( tmp.mask() )

    # now add circle
    tmp = sastbx.data_reduction.image_mask(self.hutch.nx, self.hutch.ny)
    tmp_circle = self.hutch.exclude_circle
    if tmp_circle is not None:
      tmp.add_circle(tmp_circle.x, tmp_circle.y, tmp_circle.radius)
      self.mask_object.add( tmp.mask() )

    # now get the mask definition from the image
    tmp = sastbx.data_reduction.image_mask(self.hutch.nx, self.hutch.ny)
    img_name = self.hutch.exclude_from_image
    if img_name is not None:
      tmp.detect_dead_pixels( image_as_flex_int(img_name,self.hutch.detector_type,False) )
      self.mask_object.add( tmp.mask() )

  def build_shadow(self):
    # XXX WE NEED TO DO SOMETHING ABOUT 'INSIDE/OUTSIDE' DEFINITIONS XXX #
    # first add polygons please
    for polyg in self.hutch.shadow_polygons:
      tmp = sastbx.data_reduction.image_mask(self.hutch.nx, self.hutch.ny)
      tmp.add_polygon( flex.int(polyg.corners_x), flex.int(polyg.corners_y), 0,0 ) # we need a solution for this please
      self.shadow_object.add( tmp.mask() )

    # now add circle
    tmp = sastbx.data_reduction.image_mask(self.hutch.nx, self.hutch.ny)
    tmp_circle = self.hutch.shadow_circle
    if tmp_circle is not None:
      tmp.add_circle(tmp_circle.x, tmp_circle.y, tmp_circle.radius)
      self.shadow_object.add( tmp.mask() )


  def integrate_image(self,image_name,get_shadow_data=False):
    data = image_as_flex_int( image_name,self.hutch.detector_type  )
    self.integrator_engine.integrate(data,get_shadow_data)
    mi = self.integrator_engine.mean_intensity()
    vi = self.integrator_engine.variance_intensity()
    return mi, vi

  def get_aux_data(self, image_name):
    if self.hutch.aux_info_object is None:
      before = parse_ion_chamber_keys( self.hutch.params.detector.header_keys.I0 )
      after  = parse_ion_chamber_keys( self.hutch.params.detector.header_keys.I1 )
      assert len(before) == len(after)
      before = get_total_ion_chamber_value( image_name, self.hutch.detector_type, before )
      after  = get_total_ion_chamber_value( image_name, self.hutch.detector_type, after  )
    else:
      before, after = self.hutch.aux_info_object.get_i0_i1( image_name )

    transmission_correction = after/max(before, 1e-13)
    i_total = before
    i_after = after
    return transmission_correction, i_total, i_after

  def get_q_arrays(self):
    return self.q_low, self.q_mean, self.q_step




class q_range(object):
  def __init__(self, q_low, q_mean, q_step):
    self.q_low  = q_low
    self.q_mean = q_mean
    self.q_step = q_step

class radially_integrated_data(object):
  def __init__(self, q_range, name=None):
    self.q_range = q_range
    self.name = name
    self.n_points = None
    self.mean = None
    self.var  = None
    self.thrd = None
    self.frth = None

    self.transmission_factor = None
    self.i_tot = None

  def load_data(self,fst,snd,thrd=None,frth=None):
    self.mean = fst
    self.var = snd
    self.thrd = thrd
    self.frth = frth

  def set_aux_info( self, transmission, i_tot,i_after):
    self.transmission_factor = transmission
    self.i_tot = i_tot
    self.i_after = i_after


  def multiply(self, factor):
    self.mean = self.mean*factor
    self.var = self.var*factor*factor
    if self.thrd is not None:
      self.thrd = self.thrd*factor*factor*factor
    if self.frth is not None:
      self.frth = self.frth*factor*factor*factor*factor

  def add(self, value):
    self.mean = self.mean+value

  def mean_var(self):
    return self.mean, self.var

  def write_data(self, file_name=None, out=None ):
    if out is None:
      out = sys.stdout
    if file_name is None:
      file_name = self.name
    print >> out, "   Writing:", file_name
    print >> out, "   Scaling by I1: ", 1.0/self.i_after
    #scale by I1
    scale = 1.0/self.i_after
    data = curves.simple_saxs_data(self.q_range.q_mean, self.mean*scale, flex.sqrt(self.var)*scale)
    saxs_read_write.write_standard_ascii_qis(data, file_name)     


