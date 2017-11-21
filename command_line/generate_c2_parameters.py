# LIBTBX_SET_DISPATCHER_NAME sastbx.generate_c2_parameters

import cPickle, os, random, time

from iotbx import pdb
from scitbx.random import set_random_seed

from sastbx.fXS import ewald_sphere, detector_geometry
from sastbx.fXS.image_simulator_base import detector_properties, beam_properties,\
     model_properties
from sastbx.fXS.gpu_c2_calculation import map_correlate_parameters,\
     generate_parameters
from sastbx.fXS.structure_generator import structure_generator

# =============================================================================
if (__name__ == '__main__'):

  t0 = time.time()
  set_random_seed(int(t0))

  # user-defined parameters
  # ---------------------------------------------------------------------------
  detector_size = (1 * 1024, 1 * 1024)
  pixel_size = (3 * 102 * 1e4, 3 * 102 * 1e4)  # 306 um
  distance = 2.0 * 1e10                        # 2.0 m
  wavelength = 1.0                             # 1 Angstrom
  #flux = 250e6  # 1e11 photons/s in 20 um x 20 um, bl 5.0.2
  flux = 2e12                                  # 2e12 photons/pulse
  t =  1.0                                     # 1 pulse
  box_size = 1*(10000.0)                       # 10000 A = 1 um
  min_separation = 500.0                       # 500 A

  pdb_files = ['a.pdb']
  particle_counts = [100]

  coherent = True
  use_solvent = True
  photon_poisson_noise = True
  particle_count_noise = True

  ring_indices = [] #range(0,120,10)

  c2_q_min = 0.05                              # binning parameters
  c2_q_max = 0.48
  q_pixel_depth = 2
  phi_pixel_radius = 2

  n_chunks = 2                                 # job parameters
  n_images_bins = 2
  n_images_step_size = 10
  # ---------------------------------------------------------------------------

  # parameters
  dp = detector_properties()
  dp.detector_size = detector_size
  dp.pixel_size = pixel_size
  dp.corner_position = (dp.detector_size[0]/2 * dp.pixel_size[0],
                        dp.detector_size[1]/2 * dp.pixel_size[1], distance)
  dp.distance = dp.corner_position[-1]

  bp = beam_properties()
  bp.wavelength = wavelength
  bp.flux = flux
  bp.t = t

  es = ewald_sphere()
  es.set_wavelength(bp.wavelength)

  dg = detector_geometry()
  dg.set_corner_position(dp.corner_position)
  dg.set_detector_size(dp.detector_size)
  dg.set_pixel_size(dp.pixel_size)

  s = structure_generator()
  s.use_solvent = use_solvent
  s.box_size = box_size
  s.min_separation = min_separation
  for pf in pdb_files:
    s.add_species(pdb_input=pdb.input(pf),n_copies=1)

  mp = model_properties()
  mp.detector_properties = dp
  mp.beam_properties = bp
  mp.structure_generator = s
  mp.random_state = random.getstate()

  # generate parameters for multiprocessing
  p = map_correlate_parameters()
  p.q_min = c2_q_min
  p.q_max = c2_q_max
  p.q_pixel_depth = q_pixel_depth
  p.phi_pixel_radius = phi_pixel_radius
  p.model_properties = mp
  p.coherent = coherent
  p.photon_poisson_noise = photon_poisson_noise
  p.particle_count_noise = particle_count_noise
  p.mean = particle_counts
  if ( (ring_indices is not None) and (len(ring_indices) == 0) ):
    ring_indices = None
  p.ring_indices = ring_indices

  output_directory = './raw_data/'
  if (os.path.exists(output_directory) is False):
    os.mkdir(output_directory)

  # generate parameters for all jobs
  current_image_total = 0
  current_image_bin = 0
  while (current_image_bin < n_images_bins):
    n_chunks,parameters = generate_parameters\
                          (p=p,n_images=n_images_step_size,n_cpu=n_chunks)
    current_image_total += n_images_step_size
    current_image_bin += 1
    t1 = time.time()

    for i in xrange(n_chunks):
      f = open(output_directory + '/' + str(current_image_total).zfill(8) +
               '_' + str(i).zfill(3) + '.par','wb')
      cPickle.dump(parameters[i],f,cPickle.HIGHEST_PROTOCOL)
      f.close()

  t1 = time.time()
  t_setup = t1 - t0

  print t_setup/60, 'min (setup)'
