# LIBTBX_SET_DISPATCHER_NAME sastbx.gpu_c2

import math, os, random, socket, time

from iotbx import pdb
from libtbx.smp_utils import Pool, multiprocessing_params
from scitbx.array_family import flex
from scitbx.random import set_random_seed

from sastbx.fXS import ewald_sphere, detector_geometry, c2_tiles
from sastbx.fXS.image_simulator_base import detector_properties, beam_properties,\
     model_properties
from sastbx.fXS.gpu_c2_calculation import map_correlate_parameters,\
     generate_parameters, map_correlate_gpu
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
  min_separation = 200.0                       # 500 A

  pdb_files = ['a.pdb']
  particle_counts = [100]

  coherent = True
  use_solvent = True
  photon_poisson_noise = True
  particle_count_noise = True

  ring_indices = [40]

  c2_q_min = 0.05                              # binning parameters
  c2_q_max = 0.48
  q_pixel_depth = 2
  phi_pixel_radius = 2

  n_cpu = 10                                   # job parameters
  n_images_bins = 1
  n_images_step_size = 10000
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
  #es.set_distance(dp.distance)

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

  # get binning
  ac = c2_tiles()
  ac.set_ewald_sphere(es)
  ac.add_geometry(dg)
  ac.set_q_limits(c2_q_min,c2_q_max)
  ac.set_ring_pixel_sizes(q_pixel_depth,phi_pixel_radius)
  ac.initialize()

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
  p.ring_indices = ring_indices

  n_q = ac.get_n_rings()

  ## write_image(file_name='rings.png',detector_size=dp.detector_size,
  ##             image_data=ac.bin_mask(0).as_double())

  ac_average = [ None for i in xrange(n_q) ]
  ac_variance = [0.0 for i in xrange(n_q) ]
  for i in xrange(n_q):
    ac_average[i] = flex.double(len(ac.get_c2(i)),0.0)
  sc_average = flex.double(n_q,0.0)
  sc_variance = flex.double(n_q,0.0)

  t1 = time.time()

  t_setup = t1 - t0
  t_processing = 0.0
  t_collecting = 0.0

  output_directory = './correlations/'
  if (os.path.exists(output_directory) is False):
    os.mkdir(output_directory)

  # run jobs
  current_image_total = 0
  current_image_bin = 0
  while (current_image_bin < n_images_bins):
    t0 = time.time()
    n_cpu,parameters = generate_parameters\
                       (p=p,n_images=n_images_step_size,n_cpu=n_cpu)
    current_image_total += n_images_step_size
    current_image_bin += 1
    t1 = time.time()

    sge_mp = multiprocessing_params.fetch().extract()
    sge_mp.enable_multiprocessing = True
    sge_mp.method = "sge"

    # multiple jobs
    pool = Pool(sge_mp)
    result = pool.map(map_correlate_gpu,parameters)
    del pool

    # single jobs
    ## result = [None for i in xrange(len(parameters))]
    ## for i in xrange(len(parameters)):
    ##   result[i] = map_correlate_gpu(parameters[i])

    # collect results
    t2 = time.time()
    if (len(parameters) == len(result)):
      for i in xrange(len(parameters)):
        for j in xrange(n_q):
          ac_average[j] += result[i].ac_average[j]
          ac_variance[j] += result[i].ac_variance[j]
        sc_average += result[i].sc_average
        sc_variance += result[i].sc_variance

      # average results
      hostname = socket.gethostname().split('.')[0]
      f = open(output_directory + 'ac_' + hostname + '.' +
               str(current_image_total),'w')
      sqrt_n = math.sqrt(current_image_total)
      for i in xrange(n_q):
        f.write('# q = %f\n'%ac.get_ring_q(i))
        current_average = ac_average[i]/current_image_total
        current_variance = ac_variance[i]/current_image_total -\
                           current_average*current_average
        dphi = 2.0*math.pi/len(current_average)
        for j in xrange(len(current_average)):
          if (current_variance[j] < 0):
            current_variance[j] = math.fabs(current_variance[j])
          f.write('%f %f %f\n'%(0.5*dphi + j*dphi,current_average[j],
                                math.sqrt(current_variance[j])/sqrt_n))
        f.write('&\n')
      f.close()

      f = open(output_directory + 'sc_' + hostname + '.' +
               str(current_image_total),'w')
      for i in xrange(n_q):
        current_average = sc_average[i]/current_image_total
        f.write('%f %f %f\n'%(ac.get_ring_q(i),current_average,
                              sc_variance[i]/current_image_total -\
                              current_average*current_average))
      f.close()

    else:
      print 'Warning: Jobs did not complete properly\n',\
            '         Submitted', len(parameters), 'jobs\n',\
            '         Received', len(result), 'results\n',\
            '         Re-trying frames',\
            current_image_total - n_images_step_size, 'to',\
            current_image_total, '\n'
      current_image_total -= n_images_step_size
      current_image_bin -= 1
    t3 = time.time()

    t_setup += t1 - t0
    t_processing += t2 - t1
    t_collecting += t3 - t2

  print t_setup/60, 'min (setup)'
  print t_processing/60, 'min (processing)'
  print t_collecting/60, 'min (collecting)'
