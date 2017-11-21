import copy, math, os, pickle, random, sys, time
from multiprocessing import Pool
from libtbx.introspection import number_of_processors

from iotbx import pdb
from scitbx.array_family import flex
from scitbx.random import variate, poisson_distribution, set_random_seed
from sastbx.fXS import ewald_sphere, image_composer, mean_and_variance_by_q,\
     multiple_poisson, apply_translation
from sastbx.fXS.image_simulator import image_simulator, write_image,\
     detector_properties, beam_properties
from sastbx.fXS.structure_generator import structure_generator

# =============================================================================
class map_random_image_parameters(object):
  def __init__(self):
    self.base_image_directory = None
    self.image_list = None
    self.structure_generator = None
    self.beam_properties = None
    self.detector_properties = None
    self.random_state = None
    self.n_images = None
    self.n_particles = None

# =============================================================================
def map_random_image(p=None):

  result = [None for i in xrange(p.n_images)]

  # set random state for current thread
  p.structure_generator.random.setstate(p.random_state)

  es = ewald_sphere()
  es.set_wavelength(p.beam_properties.wavelength)
  es.set_distance(p.detector_properties.distance)

  ic = image_composer()
  ic.set_detector_size(p.detector_properties.detector_size)
  ic.set_beam_center(p.detector_properties.beam_center)
  ic.set_pixel_size(p.detector_properties.pixel_size)
  ic.set_ewald_sphere(es)

  im = image_simulator()
  im.structures = p.structure_generator
  im.image_composer = ic

  h = ic.cache_h()
  sf_length = (p.detector_properties.detector_size[0] + 1) *\
              (p.detector_properties.detector_size[1] + 1)
  n_photons = p.beam_properties.flux * p.beam_properties.t

  for i in xrange(p.n_images):
    sf = flex.complex_double(sf_length,0.0)
    im.structures.species[0].n_copies = p.n_particles[i]
    im.structures.randomize()
    for j in xrange(p.n_particles[i]):
      image_name = p.base_image_directory + im.structures.random.choice(image_list)
      t = im.structures.translations[0][j]
      f = open(image_name,'rb')
      sf_tmp = pickle.load(f)
      f.close()
      sf_tmp = apply_translation(sf_tmp,h,t)
      sf += sf_tmp
    im.structure_factors = sf
    image_data = im.build_image(n_photons=n_photons)
    result[i] = image_data.deep_copy()

  return result

# =============================================================================
def generate_parameters(p=None,n_images=None,n_cpu=None,mean=None):
  pv = variate(poisson_distribution(mean))

  r = random.Random()
  r.seed()
  if (n_images > n_cpu):
    n_images_per_cpu = int(math.floor(n_images/n_cpu))
    n_cpu = int(math.ceil(n_images/n_images_per_cpu))
    parameters = [copy.deepcopy(p) for i in xrange(n_cpu)]
    remaining_images = n_images
    for i in xrange(len(parameters)):
      if (remaining_images > n_images_per_cpu):
        parameters[i].n_images = n_images_per_cpu
        remaining_images -= n_images_per_cpu
      else:
        parameters[i].n_images = remaining_images
      n_particles = [0 for k in xrange(parameters[i].n_images)]
      n_jump = 0
      for j in xrange(parameters[i].n_images):
        n_particles[j] = pv()
        n_jump += 6*n_particles[j]
      parameters[i].n_particles = copy.deepcopy(n_particles)
      parameters[i].random_state = r.getstate()
      r.jumpahead(n_jump)
  else:
    n_cpu = n_images
    parameters = [copy.deepcopy(p) for i in xrange(n_cpu)]
    for i in xrange(n_cpu):
      parameters[i].n_images = 1
      n_particles = [0 for k in xrange(parameters[i].n_images)]
      n_jump = 0
      for j in xrange(parameters[i].n_images):
        n_particles[j] = pv()
        n_jump += 6*n_particles[j]
      parameters[i].n_particles = copy.deepcopy(n_particles)
      parameters[i].random_state = r.getstate()
      r.jumpahead(n_jump)

  return n_cpu,parameters

# =============================================================================
if (__name__ == '__main__'):

  t0 = time.time()
  set_random_seed(int(t0))

  # parameters
  dp = detector_properties()
  dp.detector_size = (512,512)
  dp.beam_center = (256,256)
  dp.pixel_size = 6*(102 * 1e4)  # 102 um
  dp.distance = 2.0 * 1e10  # 1.5 m

  bp = beam_properties()
  bp.wavelength = 1.0  # 1 Angstrom
  bp.flux = 1.0e10
  bp.t = 1.0

  es = ewald_sphere()
  es.set_wavelength(bp.wavelength)
  es.set_distance(dp.distance)

  ic = image_composer()
  ic.set_detector_size(dp.detector_size)
  ic.set_beam_center(dp.beam_center)
  ic.set_pixel_size(dp.pixel_size)
  ic.set_ewald_sphere(es)

  s = structure_generator()
  s.add_species(pdb_input=pdb.input('full_16.pdb'),n_copies=1)

  im = image_simulator()
  im.structures = s
  im.image_composer = ic
  h = ic.cache_h()

  q = ic.get_q()
  n_bins = 100
  q_min = 0.005
  q_max = 0.45
  dq = (q_max - q_min)/n_bins

  # can buy 10 mg/mL lysozyme solution -> 42000 molecules/ um^3
  # 10 mg/mL * (1 g / 1000 mg) * (Avogadro's Number particles / 14331.2 g) *
  # (1 cm / 10000 um)^3
  # 10 ug/mL should be 42 molecules / um^3
  mean = 10
  n_images = 10000

  # read file names
  base = '3kfb_images'
  image_list_file = base + '.pickle'
  base_image_directory = '/scratch/home/bkpoon/' + base + '/'
  composite_image_directory = './composite_images/'
  f = open(image_list_file,'rb')
  image_list = pickle.load(f)
  f.close()

  # generate parameters for multiprocessing
  n_cpu = number_of_processors()
  p = map_random_image_parameters()
  p.base_image_directory = base_image_directory
  p.image_list = image_list
  p.structure_generator = s
  p.beam_properties = bp
  p.detector_properties = dp
  p.n_images = n_images
  n_cpu,parameters = generate_parameters(p=p,n_images=n_images,n_cpu=n_cpu,
                                         mean=mean)

  # run jobs
  #map_random_image(parameters[0]); exit()
  t1 = time.time()
  pool = Pool(processes=n_cpu)
  try:
    composite_images = pool.map(map_random_image,parameters)
    pool.close()
  except Exception:
    pool.terminate()
  finally:
    pool.join()

  # apply Poisson noise to counts
  t2 = time.time()
  sums = flex.vec2_double(n_bins,(0.0,0.0))
  image_count = 0;
  for i in xrange(len(composite_images)):
    for j in xrange(len(composite_images[i])):
      image_data = multiple_poisson(composite_images[i][j]).as_double()

      if (True):
        image_name = composite_image_directory + str(os.getpid()) + '_' +\
                     str(image_count).zfill(len(str(n_images))) + '.dat'
        f = open(image_name,'wb')
        pickle.dump(image_data,f,2)
        f.close()
        image_count += 1

      # convert image data to I(q)
      mv = mean_and_variance_by_q(q_min,q_max,dq,n_bins,q,image_data)
      sums += mv

  # output scattering curve
  f = open('./composite_data/' + str(mean) + '_' + str(n_images) + '.dat','w')
  for i in xrange(n_bins):
    q_i = q_min + 0.5*dq + i*dq
    f.write('%f %f %f\n'%(q_i,sums[i][0]/n_images,sums[i][1]/n_images))
  f.close()

  t3 = time.time()
  print (t1-t0)/60, 'min (setup)'
  print (t2-t1)/60, 'min (generation)'
  print (t3-t2)/60, 'min (processing)'
