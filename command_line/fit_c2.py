# LIBTBX_SET_DISPATCHER_NAME sastbx.fit_c2

import copy, math, os, sys
from multiprocessing import Pool

from libtbx.introspection import number_of_processors
from sastbx.fXS import fourier_legendre_series
from scitbx.array_family import flex

from sastbx.fXS.fit_data import fit_data, fit_parameters, read_ac, read_sc

if (__name__ == '__main__'):

  # user-defined parameters
  # ---------------------------------------------------------------------------
  fls_data = '/home/useraccts/mcf-staff/bkpoon/work/data/fls.dat'
  denormalize = True
  minimize = True
  standardize = False
  wavelength = 1.0
  fit_order = 40
  ring_indices = [] #range(0,120,10)
  # ---------------------------------------------------------------------------

  f = open(sys.argv[1])
  fl = f.readlines()
  f.close()

  # base parameter object
  p = fit_parameters()
  p.fls_data = fls_data
  p.minimize = minimize
  p.standardize = standardize
  p.wavelength = wavelength
  p.fit_order = fit_order

  # loop over files
  result = list()
  for filename in fl:
    q, x, ac, v = read_ac(filename=filename.rstrip())

    if ((ring_indices is None) or (len(ring_indices) == 0)):
      ring_indices = range(len(q))

    if (denormalize):
      q_sc, sc = read_sc(filename=filename.rstrip())
      assert(len(q) == q_sc.size())
      for i in ring_indices:
        ac[i] = sc[i] * sc[i] * ac[i]

    parameters = [copy.deepcopy(p) for i in xrange(len(ring_indices))]
    for i in xrange(len(ring_indices)):
      parameters[i].q = float(q[ring_indices[i]].split()[-1])
      parameters[i].x = x[ring_indices[i]]
      parameters[i].ac = ac[ring_indices[i]]
      parameters[i].v = v[ring_indices[i]]

    # fit each q independently
    np = number_of_processors()
    if (len(ring_indices) < np):
      np = len(ring_indices)
    pool = Pool(processes=np)
    try:
      result.append(pool.map(fit_data,parameters))
      pool.close()
    except Exception:
      pool.terminate()
    finally:
      pool.join()

  output_directory = './fit_' + str(p.fit_order) + '/'
  if (os.path.exists(output_directory) is False):
    os.mkdir(output_directory)

  # check fit for last file
  # q, x, ac for the last file will exist from the end of earlier loop over files
  # result[file_index][q_index][0 = coefficients, 1 = variances, 2 = scales]
  fls = fourier_legendre_series()
  fls.read_polynomials(p.fls_data)
  q_count = 0
  for q_i in ring_indices:
    f = open(output_directory + 'f_' + str(q_i).zfill(2) + '.dat','w')
    f2 = open(output_directory + 'e_' + str(q_i).zfill(2) + '.dat','w')

    cos_sq = float(q[q_i].split()[-1])*p.wavelength/(4.0*math.pi)
    cos_sq = cos_sq*cos_sq
    sin_sq = 1.0 - cos_sq
    ac_new = fls.compute_function(result[-1][q_count][0], cos_sq +
                                  sin_sq*flex.cos(x[q_i]))

    if (p.standardize):
      ac[q_i] = ac[q_i]/result[-1][q_count][2][0] - 1.0
      ac[q_i] = 1.0/result[-1][q_count][2][1] * ac[q_i] + 1.0
      ac_new = ac_new/result[-1][q_count][2][0] - 1.0
      ac_new = 1.0/result[-1][q_count][2][1] * ac_new + 1.0

    ac_error = 100.0*flex.fabs(ac_new - ac[q_i])/ac_new

    # write output
    for i in xrange(x[q_i].size()):
      f.write('%f %e\n'%(x[q_i][i], ac[q_i][i]))
    f.write('&\n')
    for i in xrange(x[q_i].size()):
      f.write('%f %e\n'%(x[q_i][i], ac_new[i]))
      f2.write('%f %e\n'%(x[q_i][i], ac_error[i]))

    f.close()
    f2.close()
    q_count += 1

    print flex.mean(ac_error), flex.median(ac_error)

  q_count = 0
  f = open(output_directory + 'b.dat','w')
  for i in ring_indices:
    f.write('%s '%q[i].split()[-1])
    for j in xrange(0,len(result[-1][q_count][0]),2):
      f.write('%f '%result[-1][q_count][0][j])
    f.write('\n')
    q_count += 1
  f.close()

  # output coefficients
  q_count = 0
  for q_i in ring_indices:
    f = open(output_directory + str(q_i).zfill(2) + '.dat','w')
    for c_i in xrange(fit_order + 1):
      for i in xrange(len(result)):
        f.write('%i %e %e\n'%(i, result[i][q_count][0][c_i],\
                              math.sqrt(result[i][q_count][1][c_i])))
      f.write('&\n')
    f.close()
    q_count += 1
