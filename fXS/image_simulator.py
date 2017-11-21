import copy, math, os, pickle, random, time
from multiprocessing import Pool

from PIL import Image, ImageOps

from cctbx import miller
from cctbx.xray import structure_factors
from iotbx import pdb
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex

from sastbx.fXS import ewald_sphere, image_composer, mean_and_variance_by_q,\
     direct_sum_structure_factors, I_q
from sastbx.fXS.structure_generator import structure_generator

import libtbx.load_env
if (libtbx.env.build_options.enable_cuda):
  from sastbx.fXS import cuda_direct_summation

# =============================================================================
class beam_properties(object):

  def __init__(self):
    self.wavelength = None
    self.flux = None
    self.t = None

class detector_properties(object):

  def __init__(self):
    self.detector_size = None
    self.beam_center = None
    self.pixel_size = None
    self.distance = None

class model_properties(object):

  def __init__(self):
    self.beam_properties = None
    self.detector_properties = None
    self.random_state = None
    self.structure_generator = None

    self.n_images = None
    self.n_bins = None
    self.q_min = None
    self.q_max = None
    self.dq = None

    self.gpu = None

  def show(self):
    print 'Detector Properties'
    print '==================='
    print 'Size =', self.detector_properties.detector_size
    print 'Center =', self.detector_properties.beam_center
    print 'Pixel size =', self.detector_properties.pixel_size
    print 'Distance =', self.detector_properties.distance
    print
    print 'Beam Properties'
    print '==============='
    print 'Wavelength =', self.beam_properties.wavelength
    print 'Flux =', self.beam_properties.flux
    print 'Time =', self.beam_properties.t
    print
    print 'Structures'
    print '=========='
    for i in xrange(len(self.structure_generator)):
      print 'Species',i,'has',self.structure_generator.species[i].n_copies,\
            'orientation(s)'
    print
    print 'Model Parameters'
    print '================'
    print 'N images =', self.n_images
    print 'N bins =', self.n_bins
    print 'q_min =', self.q_min
    print 'q_max =', self.q_max
    print 'dq =', self.dq
    print

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
    image_composer - the image_composer object containing the geometry
                     parameters

  Notes:
    Both structures and image_composer need to be set before any methods are
      called.  No checks are done.
    The algorithm argument for sum_structure_factors is either 'fft', 'direct',
      or 'None'.  When set to 'None', the structure factor calculation methods
      (cctbx.xray.structure_factors.from_scatterers) will automatically try to
      determine the best method, 'fft' or 'direct'.  Since image_composer will
      generate a list Miller indices that are used for the image, 'direct'
      seems to be the preferred method because the list is short.

  -----------------------------------------------------------------------------
  """
  def __init__(self):
    self.structures = None
    self.structure_factors = None
    self.solvent_model = None
    self.bulk_solvent_image = None
    self.image_composer = None
    self.n_cpu = 1
    self.cached_h = None
    self.cached_q = None

  def sum_structure_factors(self):

    # cache Miller indices
    if (self.cached_h is None):
      self.cached_h = self.image_composer.cache_h()
      self.cached_q = self.image_composer.get_q()

    # sum structure factors
    p = self.map_jobs()
    if (self.n_cpu > 1):
      pool = Pool(processes=self.n_cpu)
      sf = pool.map(map_structure_factor,p)
      pool.close()
      pool.join()
    else:
      sf = [ None for i in xrange(len(p)) ]
      for i in xrange(len(p)):
        sf_i = single_structure_factor(p[i][0])
        for j in xrange(1,len(p[i])):
          sf_i += single_structure_factor(p[i][j])
        sf[i] = sf_i

    self.structure_factors = sf[0]
    for i in xrange(1,len(sf)):
      if (sf[i] is not None):
        self.structure_factors += sf[i]

  def sum_structure_factors_gpu(self,gpu=0):

    # cache Miller indices
    if (self.cached_h is None):
      self.cached_h = self.image_composer.cache_h()
      self.cached_q = self.image_composer.get_q()

    # sum structure factors
    sf = cuda_direct_summation()
    for i in xrange(len(self.structures)):
      rot = flex.double()
      trans = flex.vec3_double()
      for rt in xrange(len(self.structures.rotations[i])):
        for j in xrange(len(self.structures.rotations[i][rt])):
          rot.append(self.structures.rotations[i][rt][j])
        trans.append(self.structures.translations[i][rt])
        #trans.append((0.0,0.0,0.0))
      sf.add(self.structures.species[i].scattering_types,
             self.structures.species[i].xyz,
             self.structures.species[i].\
             boundary_layer_scaling_factors,
             self.cached_h,rot,trans,
             self.structures.species[0].\
             scattering_type_registry)

    self.structure_factors = sf.get_sum()

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
    r_e = 2.818e-5  # radius of electron in Angstroms
    d = self.image_composer.get_distance()  # detector distance
    if (self.structure_factors is None):
      self.sum_structure_factors_gpu()
    if (coherent):
      intensities = flex.norm(self.structure_factors)
    else:
      intensities = self.structure_factors
    bc = self.image_composer.get_beam_center()
    ds = self.image_composer.get_detector_size()
    i000 = intensities[bc[1]*ds[0] + bc[0]]
    if (approx_equal(i000,0.0,out=None)):
      i000 = 1.0
    scale = self.structures.total_electrons/i000 * n_photons*(r_e*r_e)/(d*d)
    intensities = scale * intensities
    image = self.image_composer.build_image(intensities)
    return image

  def map_jobs(self):
    # calculate number of items per processor
    n_total = 0
    for i_structure in xrange(len(self.structures.species)):
      n_total += self.structures.species[i_structure].n_copies
    n_per_cpu = int(math.ceil(n_total/self.n_cpu))

    # returns a list of parameters required for calculating structure factor
    parameters = [ list() for i in xrange(self.n_cpu) ]
    current_cpu = 0
    for i_structure in xrange(len(self.structures.species)):
      for j_orientation in xrange(self.structures.species[i_structure].n_copies):
        parameters[current_cpu].append(
          map_structure_factor_parameters\
          (self.structures.species[i_structure].scattering_types,
           self.structures.species[i_structure].xyz,
           self.structures.species[i_structure].boundary_layer_scaling_factors,
           self.cached_h,
           self.structures.species[i_structure].scattering_type_registry,
           self.structures.rotations[i_structure][j_orientation],
           self.structures.translations[i_structure][j_orientation]))
        if (len(parameters[current_cpu]) >= n_per_cpu):
          current_cpu += 1
    return parameters

# =============================================================================
class map_structure_factor_parameters(object):
  """
  =============================================================================
  Structure Factor Parmeters Class
  Data structure to simplify the passing of information required for
  calculating the structure factor of a model

  Useful accessible data members:
    structure - the original, unmodified cctbx.xray.structure
    rotation - the rotation (tuple) to be applied to the structure
    translation - the translation (list) to be applied to the structure
    resolution - the resolution for the structure factor calculation
    miller_set - the Miller indices (cctbx.miller.set) for the structure factor
                 calculation
  -----------------------------------------------------------------------------
  """

  def __init__(self,scattering_types,xyz,boundary_layer_scaling_factors,h,
               scattering_type_registry,rotation,translation):
    self.scattering_types = scattering_types
    self.xyz = xyz
    self.boundary_layer_scaling_factors = boundary_layer_scaling_factors
    self.h = h
    self.scattering_type_registry = scattering_type_registry
    self.rotation = rotation
    self.translation = translation

# =============================================================================
def map_structure_factor(g):

  if (len(g) > 0):
    structure_factor = single_structure_factor(g[0])
    for i in xrange(1,len(g)):
      structure_factor += single_structure_factor(g[i]).data()
    return structure_factor
  return None

# =============================================================================
def single_structure_factor(p):

  # orient model
  new_xyz = p.xyz.deep_copy()
  translation = flex.double(p.translation)
  rotation = matrix.sqr(p.rotation)
  for i in xrange(len(new_xyz)):
    xyz = new_xyz[i]
    xyz = flex.double(rotation * xyz)
    new_xyz[i] = tuple(xyz + translation)

  # calculate structure factor
  structure_factor = direct_sum_structure_factors\
                     (p.scattering_types,new_xyz,
                      p.boundary_layer_scaling_factors,p.h,
                      p.scattering_type_registry)
  return structure_factor

# =============================================================================
def map_model_I_q(mp):

  # set up parameters
  mp.structure_generator.random.setstate(mp.random_state)

  es = ewald_sphere()
  es.set_wavelength(mp.beam_properties.wavelength)
  es.set_distance(mp.detector_properties.distance)

  ic = image_composer()
  ic.set_detector_size(mp.detector_properties.detector_size)
  ic.set_beam_center(mp.detector_properties.beam_center)
  ic.set_pixel_size(mp.detector_properties.pixel_size)
  ic.set_ewald_sphere(es)

  # initialize image simulator
  im = image_simulator()
  im.structures = mp.structure_generator
  im.image_composer = ic

  # make images
  sums = flex.vec2_double(mp.n_bins)
  for i in xrange(mp.n_images):
    im.structures.randomize()
    im.sum_structure_factors_gpu(gpu=mp.gpu)
    image_data = im.build_image(n_photons=(mp.beam_properties.flux*
                                           mp.beam_properties.t))
    if (False):
      for q_i in [0.01, 0.1, 0.2]:
        print q_i, I_q(im.cached_q,image_data,q_i,mp.dq)
      write_image(detector_size=mp.detector_properties.detector_size,
                  image_data=image_data)
    if (True):
      file_name = '/scratch/home/bkpoon/images/' + str(os.getpid()) + '_' + str(i)
      f = open(file_name,'wb')
      pickle.dump(im.structure_factors,f,2)
      f.close()

    # convert image data to I(q)
    mv = mean_and_variance_by_q(mp.q_min,mp.q_max,mp.dq,mp.n_bins,
                                im.cached_q,image_data)
    sums += mv

  return sums

# =============================================================================
def write_image(file_name='image.png',detector_size=None,image_data=None,
                max_value=2.0**8-1):
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

# =============================================================================
def generate_parameters(mp=None,n_images=None,n_cpu=None):
  n_structures = 0
  for i in xrange(len(mp.structure_generator)):
    n_structures += mp.structure_generator.species[i].n_copies

  r = random.Random()
  r.seed()
  if (n_images > n_cpu):
    n_images_per_cpu = int(math.floor(n_images/n_cpu))
    n_cpu = int(math.ceil(n_images/n_images_per_cpu))
    parameters = [copy.deepcopy(mp) for i in xrange(n_cpu)]
    remaining_images = n_images
    gpu = 0
    max_gpu = 4
    for i in xrange(len(parameters)):
      parameters[i].gpu = gpu
      gpu += 1
      if (gpu == max_gpu):
        gpu = 0
      if (remaining_images > n_images_per_cpu):
        parameters[i].n_images = n_images_per_cpu
        remaining_images -= n_images_per_cpu

        # move random state
        # one rotation requires 3 random numbers and one translation requires 3
        # random numbers, so 6 calls per structure
        # use 7 to be safe
        parameters[i].random_state = r.getstate()
        r.jumpahead(7*n_structures*n_images_per_cpu)
      else:
        parameters[i].n_images = remaining_images
        parameters[i].random_state = r.getstate()
        r.jumpahead(7*n_structures*n_images_per_cpu)
        break
  else:
    n_cpu = n_images
    parameters = [copy.deepcopy(mp) for i in xrange(n_cpu)]
    for i in xrange(n_cpu):
      parameters[i].gpu = i
      parameters[i].n_images = 1
      parameters[i].random_state = r.getstate()
      r.jumpahead(7*n_structures)
  return n_cpu,parameters

# =============================================================================
def run(args):
  t0 = time.time()

  # define parameters
  dp = detector_properties()
  ## dp.detector_size = (3072,3072)
  ## dp.beam_center = (1536,1536)
  ## dp.detector_size = (1024,1024)
  ## dp.beam_center = (512,512)
  dp.detector_size = (512,512)
  dp.beam_center = (256,256)
  dp.pixel_size = 6*(102 * 1e4)  # 102 um
  dp.distance = 2.0 * 1e10  # 2.0 m

  bp = beam_properties()
  bp.wavelength = 1.0  # 1 Angstrom
  bp.flux = 1.0e10
  bp.t = 1.0

  mp = model_properties()
  mp.detector_properties = dp
  mp.beam_properties = bp

  mp.n_images = 0
  mp.n_bins = 100
  mp.q_min = 0.005
  mp.q_max = 0.45
  mp.dq = (mp.q_max - mp.q_min)/mp.n_bins

  s = structure_generator()
  s.use_solvent = True
  s.box_size = 1*(10000.0)
  s.min_separation = 500.0 #1*(10000.0)
  s.add_species(pdb_input=pdb.input(file_name='b.pdb'),n_copies=1)
  #s.add_species(pdb_input=pdb.input(file_name='6LYZ.pdb'),n_copies=1)
  #s.add_species(pdb_input=pdb.input(file_name='full_16.pdb'),n_copies=1)
  mp.structure_generator = s

  # construct list of parameters
  n_images = 20000
  n_cpu = 8
  n_cpu,parameters = generate_parameters(mp=mp,n_images=n_images,
                                         n_cpu=n_cpu)

  # run jobs
  #mv = map_model_I_q(parameters[0]); exit()
  t1 = time.time()
  pool = Pool(processes=n_cpu)
  try:
    mv = pool.map(map_model_I_q,parameters)
    pool.close()
  except Exception:
    pool.terminate()
  finally:
    pool.join()

  # calculate means over all images
  t2 = time.time()
  mv_sum = flex.vec2_double(mp.n_bins)
  for i in xrange(len(mv)):
    mv_sum += mv[i]

  f = open('./data/' + str(s.species[0].n_copies) + '_' + str(n_images) +\
           '.' + args[0],'w')
  for i in xrange(mp.n_bins):
    q_i = mp.q_min + 0.5*mp.dq + i*mp.dq
    f.write('%f %f %f\n'%(q_i,mv_sum[i][0]/n_images,mv_sum[i][1]/n_images))

  f.close()

  t3 = time.time()
  print (t1-t0)/60, 'min (setup)'
  print (t2-t1)/60, 'min (structure factors)'
  print (t3-t2)/60, 'min (binning)'
