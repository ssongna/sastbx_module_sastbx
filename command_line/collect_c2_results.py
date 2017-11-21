# LIBTBX_SET_DISPATCHER_NAME sastbx.collect_c2_results

import cPickle, math, os, time
from scitbx.array_family import flex

# =============================================================================
if (__name__ == '__main__'):

  t0 = time.time()

  # file locations
  input_directory = './raw_data/'
  output_directory = './correlations/'
  if (os.path.exists(output_directory) is False):
    os.mkdir(output_directory)

  # sort filenames
  all_filenames = sorted(os.walk(input_directory).next()[-1])
  dat_files = list()
  par_files = list()
  for current_filename in all_filenames:
    split_name = current_filename.split('.')
    if (len(split_name) >= 2):
      if (split_name[-1] == 'dat'):
        dat_files.append(current_filename)
      elif (split_name[-1] == 'par'):
        par_files.append(current_filename)

  # collect results
  f = open(input_directory + dat_files[0],'rb')
  result = cPickle.load(f)
  q = result.ring_q
  f.close()
  image_total = int(dat_files[0].split('_')[0])
  ac_average = result.ac_average
  ac_variance = result.ac_variance
  sc_average = result.sc_average
  sc_variance = result.sc_variance
  for filename in dat_files[1:]:
    current_image_number = int(filename.split('_')[0])
    # output results after all chunks have been read for current image_total
    # or last filename is read
    if ( (current_image_number != image_total) or
         (filename == dat_files[-1]) ):
      if (filename == dat_files[-1]):
        image_total = current_image_number
      f = open(output_directory + str(image_total).zfill(8) + '.ac','w')
      sqrt_n = math.sqrt(image_total)
      for i in xrange(len(q)):
        f.write('#q = %f\n'%q[i])
        current_average = (ac_average[i].deep_copy())/image_total
        if (type(ac_variance[i]) == type(current_average)):
          current_variance = (ac_variance[i].deep_copy())/image_total -\
                             current_average*current_average
        else:
          current_variance = flex.double(len(current_average),0.0)
        dphi = 2.0*math.pi/len(current_average)
        for j in xrange(len(current_average)):
          f.write('%f %g %g\n'%(0.5*dphi + j*dphi,current_average[j],
                                math.sqrt(current_variance[j])/sqrt_n))
        f.write('&\n')
      f.close()
      f = open(output_directory + str(image_total).zfill(8) + '.sc','w')
      for i in xrange(len(q)):
        current_average = sc_average[i]/image_total
        current_variance = sc_variance[i]/image_total -\
                           current_average*current_average
        f.write('%f %g %g\n'%(q[i],current_average,
                              math.sqrt(current_variance)/sqrt_n))
      f.close()
      image_total = current_image_number
    # add chunk to results
    f = open(input_directory + filename,'rb')
    result = cPickle.load(f)
    f.close()
    for i in xrange(len(q)):
      ac_average[i] += result.ac_average[i]
      ac_variance[i] += result.ac_variance[i]
    sc_average += result.sc_average
    sc_variance += result.sc_variance

  # check for missing files
  f = open('missing.dat','w')
  for p in par_files:
    missing = True
    p = p.split('.')[0]
    for d in dat_files:
      d = d.split('.')[0]
      if (p == d):
        missing = False
        break
    if (missing):
      f.write('%s\n'%(input_directory + p + '.par'))
  f.close()

  t1 = time.time()
  t_collecting = t1 - t0

  print t_collecting/60, 'min (collecting)'
