import copy, math, random

from scitbx.array_family import flex
from scitbx.random import variate, poisson_distribution

from sastbx.fXS import detector_geometry, image_base, ewald_sphere,\
     multiple_poisson, solvent_image,\
     set_negative_to_zero, c2_tiles
from sastbx.fXS.image_simulator_base import image_simulator

# =============================================================================
class map_correlate_results():
  def __init__(self):
    self.ring_q = None
    self.ac_average = None
    self.ac_variance = None
    self.sc_average = None
    self.sc_variance = None

# =============================================================================
class map_correlate_parameters():
  def __init__(self):
    self.model_properties = None
    self.n_particles = None
    self.mean = None
    self.q_min = None
    self.q_max = None
    self.q_pixel_depth = None
    self.phi_pixel_radius = None
    self.coherent = True
    self.photon_poisson_noise = False
    self.particle_count_noise = False
    self.ring_indices = None

# =============================================================================
def generate_parameters(p=None,n_images=None,n_cpu=None):

  # determine number of images per thread
  if (n_images > n_cpu):
    n_images_per_cpu = int(math.floor(n_images/n_cpu))
    n_cpu = int(math.ceil(n_images/n_images_per_cpu))
    parameters = [copy.deepcopy(p) for i in xrange(n_cpu)]
    remaining_images = n_images
    for i in xrange(len(parameters)):
      if (remaining_images > n_images_per_cpu):
        parameters[i].model_properties.n_images = n_images_per_cpu
        remaining_images -= n_images_per_cpu
      else:
        parameters[i].model_properties.n_images = remaining_images
  else:
    n_cpu = n_images
    parameters = [copy.deepcopy(p) for i in xrange(n_cpu)]
    for i in xrange(n_cpu):
      parameters[i].model_properties.n_images = 1

  # jumble random state for each thread
  r = random.Random()
  r.setstate(p.model_properties.random_state)
  pv = list()
  for m in p.mean:
    pv.append(variate(poisson_distribution(m)))
  for i in xrange(len(parameters)):
    n_jump = 0
    parameters[i].n_particles = []
    for j in xrange(len(p.mean)):
      if (p.particle_count_noise):
        parameters[i].n_particles.append\
          (pv[j](parameters[i].model_properties.n_images))
      else:
        parameters[i].n_particles.append(
          flex.int(parameters[i].model_properties.n_images,p.mean[j]))
      for k in xrange(parameters[i].model_properties.n_images):
        n_jump += 7*parameters[i].n_particles[j][k]
      n_jump += parameters[i].model_properties.n_images
    parameters[i].model_properties.random_state = r.getstate()
    r.jumpahead(n_jump)
  p.model_properties.random_state = r.getstate()

  return n_cpu,parameters

# =============================================================================
def map_correlate_gpu(p=None):

  maxint = 2147483647

  # construct image objects
  es = ewald_sphere()
  es.set_wavelength(p.model_properties.beam_properties.wavelength)
  es.set_distance(p.model_properties.detector_properties.distance)

  dg = detector_geometry()
  dg.set_corner_position(p.model_properties.detector_properties.corner_position)
  dg.set_detector_size(p.model_properties.detector_properties.detector_size)
  dg.set_pixel_size(p.model_properties.detector_properties.pixel_size)

  ib = image_base()
  ib.set_detector_geometry(dg)
  ib.set_ewald_sphere(es)

  h = ib.get_corner_h()
  q = ib.get_center_q()
  phi = ib.get_center_phi()

  im = image_simulator()
  im.cached_h = h
  im.cached_q = q
  im.structures = p.model_properties.structure_generator
  im.structures.random.setstate(p.model_properties.random_state)
  im.image_base = ib
  n_photons = p.model_properties.beam_properties.flux *\
              p.model_properties.beam_properties.t
  if (im.structures.use_solvent):
    im.bulk_solvent_image = solvent_image(q)
    r_e = 2.818e-5  # radius of electron in Angstroms
    d = p.model_properties.detector_properties.distance
    i000 = flex.max(im.bulk_solvent_image)
    # 1 g/cm^3 (1 mol/18 g) (N_A HOH/1 mol) (10 e/1 HOH) (1 cm^3/1e24 A^3)
    # 0.33456 e / A^3
    particle_volume = 0.0
    for i in xrange(len(im.structures.species)):
      particle_volume += im.structures.species[i].n_copies*(4.0/3.0)*math.pi*\
                         math.pow(im.structures.species[i].radius,3)
    bulk_solvent_electrons = 0.33456 *(im.structures.box_size *\
                                       im.structures.box_size *\
                                       im.structures.box_size - particle_volume)
    pa = p.model_properties.detector_properties.pixel_size[0] *\
         p.model_properties.detector_properties.pixel_size[1]
    scale = bulk_solvent_electrons/i000 * n_photons*(r_e*r_e)/(d*d) * pa
    im.bulk_solvent_image = scale * im.bulk_solvent_image

  # results array
  ac = c2_tiles()
  ac.set_ewald_sphere(es)
  ac.add_geometry(dg)
  ac.set_q_limits(p.q_min,p.q_max)
  ac.set_ring_pixel_sizes(p.q_pixel_depth,p.phi_pixel_radius)
  ac.initialize()
  n_q = ac.get_n_rings()
  ac_average = [ None for i in xrange(n_q) ]
  ac_variance = [0.0 for i in xrange(n_q) ]
  for i in xrange(n_q):
    ac_average[i] = flex.double(len(ac.get_c2(i)),0.0)
  sc_average = flex.double(n_q,0.0)
  sc_variance = flex.double(n_q,0.0)

  # select rings if available
  if (p.ring_indices is not None):
    tile = 0
    ring_pixel_indices = flex.int()
    for ring_index in p.ring_indices:
      ring_pixel_indices.extend(ac.get_pixel_indices(ring_index,tile))
    ring_h = flex.vec3_double(ring_pixel_indices.size())
    for i in xrange(ring_pixel_indices.size()):
      ring_h[i] = h[ring_pixel_indices[i]]
    im.cached_h = ring_h

  # construct composite images and process
  for i in xrange(p.model_properties.n_images):
    for j in xrange(len(p.mean)):
      im.structures.species[j].n_copies = p.n_particles[j][i]
    im.structures.randomize()
    image_data = im.build_image(n_photons=n_photons,coherent=p.coherent)

    # copy rings into image
    if (p.ring_indices is not None):
      ring_image_data = image_data.deep_copy()
      image_data = flex.double(h.size(),flex.min(ring_image_data))
      for j in xrange(ring_pixel_indices.size()):
        image_data[ring_pixel_indices[j]] = ring_image_data[j]
    image_data = ib.integrate(image_data)

    if (im.structures.use_solvent):
      image_data = image_data + im.bulk_solvent_image

    # apply Poisson noise to photon count
    if (p.photon_poisson_noise):
      image_data = multiple_poisson\
                   (image_data,
                    im.structures.random.randint(0,maxint)).as_double()

    # subtract background solvent
    if (im.structures.use_solvent):
      image_data = image_data - im.bulk_solvent_image
      image_data = set_negative_to_zero(image_data)

    # process image
    ac.reset_intensities()
    ac.add_intensities(image_data)
    if (p.ring_indices is None):
      ac.process_intensities()
      for j in xrange(n_q):
        current_c2 = ac.get_c2(j)
        ac_average[j] += current_c2
        ac_variance[j] += current_c2*current_c2
    else:
      for ring_index in p.ring_indices:
        ac.process_ring(ring_index)
        current_c2 = ac.get_c2(ring_index)
        ac_average[ring_index] += current_c2
        ac_variance[ring_index] += current_c2*current_c2
    current_sc = ac.get_mean_ring_intensities()
    sc_average += current_sc
    sc_variance += current_sc * current_sc

    ## write_image(file_name='test.png',
    ##             detector_size=p.model_properties.detector_properties.detector_size,
    ##             image_data=image_data)
    ## write_image(file_name='rings.png',
    ##             detector_size=p.model_properties.detector_properties.detector_size,
    ##             image_data=ac.bin_mask().as_double())

  ring_q = flex.double(n_q)
  for i in xrange(n_q):
    ring_q[i] = ac.get_ring_q(i)

  result = map_correlate_results()
  result.ring_q = ring_q
  result.ac_average = ac_average
  result.ac_variance = ac_variance
  result.sc_average = sc_average
  result.sc_variance = sc_variance
  return result
