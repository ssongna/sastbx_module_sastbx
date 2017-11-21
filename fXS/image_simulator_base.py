try:
  from PIL import Image, ImageOps
  pil_is_available = True
except Exception:
  pass

from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex

from sastbx.fXS import direct_sum_structure_factors
gpu_available = False

import libtbx.load_env
if (libtbx.env.build_options.enable_cuda):
  from sastbx.fXS import cuda_direct_summation
  gpu_available = True
else:
  def cuda_direct_summation():
    pass

# =============================================================================
class beam_properties(object):

  def __init__(self):
    self.wavelength = None
    self.flux = None
    self.t = None

class detector_properties(object):

  def __init__(self):
    self.detector_size = None
    self.corner_position = None
    self.pixel_size = None

class model_properties(object):

  def __init__(self):
    self.beam_properties = None
    self.detector_properties = None
    self.random_state = None
    self.structure_generator = None
    self.n_images = None

# =============================================================================
class image_simulator(object):
  """
  =============================================================================
  Image Simulator Class
  Simulates a scattering pattern given some structures

  Arguments:
    None

  Useful accessible methods:
    sum_structure_factors - calculates the total structure factors for all
                            structures in all orientations
    make_image - returns an array of scaled pixel intensities

  Useful accessible data members:
    structures - the structure_generator object containing all the structures
                 and orientations
    structure_factors - the total structure factors are stored after calling
                        sum_structure_factors
    image_base - the image_base object containing the geometry parameters

  Notes:
    Both structures and image_base need to be set before any methods are
      called.  No checks are done.
    The algorithm argument for sum_structure_factors is either 'fft', 'direct',
      or 'None'.  When set to 'None', the structure factor calculation methods
      (cctbx.xray.structure_factors.from_scatterers) will automatically try to
      determine the best method, 'fft' or 'direct'.  Since image_base will
      generate a list Miller indices that are used for the image, 'direct'
      seems to be the preferred method because the list is short.

  -----------------------------------------------------------------------------
  """
  def __init__(self):
    self.structures = None
    self.structure_factors = None
    self.intensities = None
    self.solvent_model = None
    self.bulk_solvent_image = None
    self.image_base = None
    self.n_cpu = 1
    self.cached_h = None
    self.cached_q = None

  def sum_structure_factors_cpu(self):

    self.structure_factors = flex.complex_double(len(self.cached_h),
                                                 complex(0.0,0.0))

    # sum structure factors
    for i_structure in xrange(len(self.structures.species)):
      for j_orientation in xrange(self.structures.species[i_structure].n_copies):
        tmp_xyz = self.structures.species[i_structure].xyz.deep_copy()
        translation = flex.double(self.structures.translations[i_structure]
                                  [j_orientation])
        rotation = matrix.sqr(self.structures.rotations[i_structure]
                              [j_orientation])
        for i in xrange(len(tmp_xyz)):
          xyz = flex.double(rotation * tmp_xyz[i])
          tmp_xyz[i] = tuple(xyz + translation)
        self.structure_factors += direct_sum_structure_factors\
                                  (self.structures.species[i_structure].
                                   scattering_types,
                                   tmp_xyz,
                                   self.structures.species[i_structure].
                                   boundary_layer_scaling_factors,
                                   self.cached_h,
                                   self.structures.species[i_structure].
                                   scattering_type_registry)

  def incoherent_sum_structure_factors_cpu(self):

    self.intensities = flex.complex_double(len(self.cached_h),complex(0.0,0.0))

    # sum structure factors incoherently
    for i_structure in xrange(len(self.structures.species)):
      for j_orientation in xrange(self.structures.species[i_structure].n_copies):
        tmp_xyz = self.structures.species[i_structure].xyz.deep_copy()
        translation = flex.double(self.structures.translations[i_structure]
                                  [j_orientation])
        rotation = matrix.sqr(self.structures.rotations[i_structure]
                              [j_orientation])
        for i in xrange(len(tmp_xyz)):
          xyz = flex.double(rotation * tmp_xyz[i])
          tmp_xyz[i] = tuple(xyz + translation)
        structure_factors = direct_sum_structure_factors\
                            (self.structures.species[i_structure].
                             scattering_types,
                             tmp_xyz,
                             self.structures.species[i_structure].
                             boundary_layer_scaling_factors,
                             self.cached_h,
                             self.structures.species[i_structure].
                              scattering_type_registry)
        self.intensities += flex.norm(structure_factors)

  def sum_structure_factors_gpu(self):

    # sum structure factors
    sf = cuda_direct_summation()
    for i in xrange(len(self.structures)):
      rot = flex.double()
      trans = flex.vec3_double()
      for rt in xrange(len(self.structures.rotations[i])):
        for j in xrange(len(self.structures.rotations[i][rt])):
          rot.append(self.structures.rotations[i][rt][j])
        trans.append(self.structures.translations[i][rt])
      sf.add(self.structures.species[i].scattering_types,
             self.structures.species[i].xyz,
             self.structures.species[i].\
             boundary_layer_scaling_factors,
             self.cached_h,rot,trans,
             self.structures.species[i].\
             scattering_type_registry)

    self.structure_factors = sf.get_sum()

  def incoherent_sum_structure_factors_gpu(self):

    self.intensities = flex.double(len(self.cached_h),0.0)

    # sum structure factors incoherently
    for i in xrange(len(self.structures)):
      for rt in xrange(len(self.structures.rotations[i])):
        sf = cuda_direct_summation()
        rot = flex.double()
        trans = flex.vec3_double()
        for j in xrange(len(self.structures.rotations[i][rt])):
          rot.append(self.structures.rotations[i][rt][j])
        trans.append(self.structures.translations[i][rt])
        sf.add(self.structures.species[i].scattering_types,
               self.structures.species[i].xyz,
               self.structures.species[i].\
               boundary_layer_scaling_factors,
               self.cached_h,rot,trans,
               self.structures.species[i].\
               scattering_type_registry)
        self.intensities += flex.norm(sf.get_sum())

  def build_image(self,n_photons=1e20,coherent=True):
    """
    From the International Tables of Crystallography (2006). Vol. C, Chapter 2.6,
    the scattered intensity from one electron is given the by Thomson formula

                    r_e^2   / 1 + cos^2(2 theta) \
         I_e = I_0 ------- | -------------------- |          (Eq 2.6.1.1)
                     r^2    \         2          /

    where r_e is the classical radius of the electron (2.818 x 10^(-15) m),
    r is the distance between the sample and the detector, and 2 theta is
    the scattering angle.  The total scattering intensity is thus I_e
    multiplied by the total number of scattering electrons.
    """
    # cache Miller indices
    if (self.cached_h is None):
      self.cached_h = self.image_base.get_corner_h()
      self.cached_q = self.image_base.get_center_q()

    # calculate intensities
    if (gpu_available):
      if (coherent):
        self.sum_structure_factors_gpu()
        self.intensities = flex.norm(self.structure_factors)
      else:
        self.incoherent_sum_structure_factors_gpu()
    else:
      if (coherent):
        self.sum_structure_factors_cpu()
        self.intensities = flex.norm(self.structure_factors)
      else:
        self.incoherent_sum_structure_factors_cpu()

    # scale intensities
    r_e = 2.818e-5  # radius of electron in Angstroms
    d = self.image_base.get_ewald_sphere().get_distance()  # detector distance
    i000 = 0.0
    for i in xrange(len(self.structures)):
      i000 += self.structures.species[i].scattering_type_registry.\
              sum_of_scattering_factors_at_diffraction_angle_0() *\
              self.structures.species[i].n_copies
    if (approx_equal(i000,0.0,out=None)):
      i000 = 1.0
    scale = self.structures.total_electrons/i000 * n_photons*(r_e*r_e)/(d*d)
    self.intensities = scale * self.intensities
    return self.intensities

# =============================================================================
def write_image(file_name='image.png',detector_size=None,image_data=None,
                max_value=2.0**8-1):
  if (pil_is_available):
    assert((detector_size[0]*detector_size[1]) == len(image_data))
    working_image_data = image_data.deep_copy()
    min_image_value = flex.min(working_image_data)
    if (min_image_value > 0.0):
      working_image_data = flex.fabs(working_image_data / min_image_value - 1.0)
    max_image_value = flex.max(working_image_data)
    working_image_data = max_value/max_image_value * working_image_data
    image = Image.new('L',detector_size)
    image.putdata(working_image_data)
    image = ImageOps.invert(image)
    image.save(file_name)
