# LIBTBX_SET_DISPATCHER_NAME sastbx.run_c2_jobs

import cPickle, sys, time

from sastbx.fXS.gpu_c2_calculation import map_correlate_gpu

# =============================================================================
if (__name__ == '__main__'):

  if (len(sys.argv) != 2):
    print 'usage: sastbx.run_c2_jobs <parameter file>\n'
    exit()

  t0 = time.time()

  filename = sys.argv[1]
  print filename

  f = open(filename,'rb')
  p = cPickle.load(f)
  f.close()

  result = map_correlate_gpu(p)

  f = open(filename.split('.')[0] + '.dat','wb')
  cPickle.dump(result,f,cPickle.HIGHEST_PROTOCOL)
  f.close()

  t1 = time.time()
  t_processing = t1 - t0

  print t_processing/60, 'min (processing)'
