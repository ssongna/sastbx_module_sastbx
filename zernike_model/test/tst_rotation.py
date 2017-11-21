import sys, os, time
from scitbx.array_family import flex
from sastbx.zernike_model import correlation
from stdlib import math as smath
from libtbx import easy_pickle
from scitbx import math, fftpack


def tst_rotation(args):
  filename = args[0]
  filename2 = args[1]
  beta = float( args[2] )
  ngrid=40
  nmax = 20
  nlm_array=math.nlm_array(nmax)
  coefs = easy_pickle.load(filename)
  nlm_array.load_coefs(nlm_array.nlm(), coefs)

  this_nlm_array=math.nlm_array(nmax)
  coefs = easy_pickle.load(filename2)
  this_nlm_array.load_coefs(nlm_array.nlm(), coefs)

  beta=smath.pi*beta
  cc_obj = correlation( nlm_array, this_nlm_array, nmax, beta)


  fft_input = flex.complex_double( flex.grid(2*nmax+1, 2*nmax+1) )
  count = 0
  radian = 180.0/smath.pi
  out = open("scan_"+str(beta)+".dat", 'w')
  for ii in range(ngrid+1):
    for jj in range(ngrid+1):
      alpha = ii*smath.pi*2.0/ngrid
      gama  = jj*smath.pi*2.0/ngrid
      cc=cc_obj.calc_correlation( alpha, beta, gama)
      fft_input[count] = cc
      count = count + 1
      print>>out, alpha*radian, gama*radian, abs(cc)
    print>>out
  out.close()
  fft = fftpack.complex_to_complex_2d( 2*nmax+1, 2*nmax+1 )

  result = fft.forward( fft_input )
  #return result

  result = fft.backward( result )

  out = open("fft_fake_"+str(beta)+".dat", 'w')
  for ii in range( 2*nmax+1 ):
    for jj in range( 2*nmax+1 ):
      print>>out, ii*9, jj*9, abs(result[(jj,ii)] )
    print>>out
  out.close()

def tst_rotation_fft( args ):
  filename = args[0]
  filename2 = args[1]
  beta = float( args[2] )
  nmax = 20
  nlm_array=math.nlm_array(nmax)
  coefs = easy_pickle.load(filename)
  nlm_array.load_coefs(nlm_array.nlm(), coefs)

  this_nlm_array=math.nlm_array(nmax)
  coefs = easy_pickle.load(filename2)
  this_nlm_array.load_coefs(nlm_array.nlm(), coefs)

  beta=smath.pi*beta
  cc_obj = correlation( nlm_array, this_nlm_array, nmax, beta)
  mm= cc_obj.mm_coef(0)
  fft_input = flex.complex_double( flex.grid(2*nmax+1, 2*nmax+1) )
  fft_r = tst_rotation( args )
  for ii in range( mm.size() ):
    fft_input[ii] = mm[ii]
  #  print ii, mm[ii], fft_r[ii]/1681.0

  fft = fftpack.complex_to_complex_2d( 2*nmax+1, 2*nmax+1 )

  result = fft.backward( fft_input )

  out = open("fft_"+str(beta)+".dat", 'w')
  for ii in range( 2*nmax+1 ):
    for jj in range( 2*nmax+1 ):
      print>>out, ii*9, jj*9, abs(result[(ii,jj)] )
    print>>out
  out.close()



if __name__ == "__main__":
  t1 = time.time()
  args=sys.argv[1:]
  if(len(args) < 3):
    print "Usage: tst_rotation.py filename filename2 beta"
    exit()
  #tst_rotation(args)
  tst_rotation_fft(args)
  t2 = time.time()
  print "#OK","total time: ", t2-t1
