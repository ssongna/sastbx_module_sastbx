import os, pickle, random, time
from multiprocessing import Pool

from iotbx import pdb
from libtbx.introspection import number_of_processors
from scitbx.array_family import flex
from sastbx.fXS import ewald_sphere, image_composer, I_q
from sastbx.fXS.image_simulator import image_simulator, write_image,\
     detector_properties, beam_properties, model_properties,\
     generate_parameters
from sastbx.fXS.structure_generator import structure_generator

# =============================================================================
def map_calculate_rings(mp):

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

  im = image_simulator()
  im.structures = mp.structure_generator
  im.image_composer = ic
  im.cached_h = mp.h

  # calculate structure factors
  for i in xrange(mp.n_images):
    im.structures.randomize()
    im.sum_structure_factors()

    # copy intensities to full image
    sf = im.structure_factors
    all_sf = flex.complex_double((mp.detector_properties.detector_size[0]+1)*\
                                 (mp.detector_properties.detector_size[1]+1),\
                                 complex(0.0,0.0))
    k = 0
    for j in xrange(len(mp.use_index)):
      if (mp.use_index[j]):
        all_sf[j] = sf[k]
        k += 1
    im.structure_factors = all_sf

    image_data = im.build_image(n_photons=bp.flux*bp.t)

    if (True):
      for q_i in [0.01, 0.1, 0.2]:
        print q_i, I_q(ic.get_q(),image_data,q_i,mp.dq)
      write_image(file_name='ring.png',
                  detector_size=mp.detector_properties.detector_size,
                  image_data=image_data)

    if (False):
      file_name = './test_rings/' + str(os.getpid()) + '_' + str(i)
      f = open(file_name,'wb')
      pickle.dump(image_data,f,2)
      f.close()

  return mp.n_images

# =============================================================================
if (__name__ == '__main__'):

  t0 = time.time()

  # parameters
  dp = detector_properties()
  dp.detector_size = (512,512)
  dp.beam_center = (256,256)
  dp.pixel_size = 6*(102 * 1e4)  # 102 um
  dp.distance = 1.5 * 1e10  # 1.5 m

  bp = beam_properties()
  bp.wavelength = 1.0  # 1 Angstrom
  bp.flux = 1.0e5
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
  s.add_species(pdb_input=pdb.input('6LYZ.pdb'),n_copies=2)
  mp.structure_generator = s

  rings = [0.01, 0.1, 0.2]

  es = ewald_sphere()
  es.set_wavelength(bp.wavelength)
  es.set_distance(dp.distance)

  ic = image_composer()
  ic.set_detector_size(dp.detector_size)
  ic.set_beam_center(dp.beam_center)
  ic.set_pixel_size(dp.pixel_size)
  ic.set_ewald_sphere(es)

  q = ic.get_q()

  # flag pixel corners for each ring
  pixel_xy = list()
  all_h = ic.cache_h()
  use_index = flex.bool(len(all_h),False)
  use_index[dp.beam_center[1]*dp.detector_size[0] + dp.beam_center[0]] = True
  n_x = dp.detector_size[0] + 1
  for x in xrange(dp.detector_size[0]):
    for y in xrange(dp.detector_size[1]):
      q_i = y*dp.detector_size[0] + x
      h_i = y*n_x + x
      for q_ring in rings:
        ring_min = q_ring - 0.5*mp.dq
        ring_max = q_ring + 0.5*mp.dq
        if ((q[q_i] >= ring_min) and (q[q_i] < ring_max)):
          pixel_xy.append((x,y))
          use_index[h_i] = True
          use_index[y*n_x + (x+1)] = True
          use_index[(y+1)*n_x + x] = True
          use_index[(y+1)*n_x + (x+1)] = True

  # select reciprocal space vectors
  ring_h = flex.vec3_double()
  ring_h.reserve(use_index.count(True))
  for i in xrange(len(all_h)):
    if (use_index[i]):
      ring_h.append(all_h[i])

  # calculate intensities
  n_images = 50
  n_cpu = number_of_processors()
  n_cpu,parameters = generate_parameters(mp=mp,n_images=n_images,
                                         n_cpu=n_cpu)
  for p in parameters:
    p.h = ring_h
    p.use_index = use_index

  if (False):
    output = map_calculate_rings(parameters[0])
    exit()

  # run jobs
  t1 = time.time()
  pool = Pool(processes=n_cpu)
  output = pool.map(map_calculate_rings,parameters)
  pool.close()
  pool.join()

  t2 = time.time()
  print (t1-t0)/60, 'min (setup)'
  print (t2-t1)/60, 'min (structure factors)'
