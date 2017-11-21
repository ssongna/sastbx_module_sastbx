from scitbx.array_family import flex
from iotbx import pdb
from scitbx import math, fftpack
from sastbx.fXS import projection
from sastbx.fXS import image_tools
from sastbx import fXS
from stdlib import math as smath
import os, sys


def from_I_to_C2( I_image, Nq, Nphi, N=100, further_reduct=True, smoothen=False):
  nx=N
  ny=nx
  image_cart = fXS.image_cartesian( nx, ny, I_image )
  image_cart.calc_c2_array( Nq, Nphi, further_reduct, smoothen )
  c2_image=image_cart.get_c2_array()
  normalized_sp = image_cart.get_normalized_sp()
  return c2_image, normalized_sp



def tst_zernike_expansion(pdbfile, nmax=10):
  projection_obj = projection.projection( pdbfile )
  image = projection_obj.project()
  ori_image=image.deep_copy()
  output=open('ori.dat', 'w')

  # get moments
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )
  # use moments to reconstruct image
  reconst = image_tools.generate_image2(nmax,zernike_2d_mom.moments(),ori_image,N)

  for pp, qq in zip(image,reconst):
    print pp[0],pp[1],pp[2],qq[2]



def tst_fft_proj(pdbfile, nmax=20):
  dx=1.0
  projection_obj = projection.projection( pdbfile, dx=dx, fraction=0.5)
  dq = projection_obj.get_dq()
  image = projection_obj.project()
  output=open('prj1.dat', 'w')
  for pt in image:
    print >> output, pt[0],pt[1],pt[2]
  output.close()

  values=projection_obj.get_value()
  np = projection_obj.get_np()

  flex_grid=flex.grid(np,np)
  fft_input=flex.complex_double(flex_grid)
  for ii in range( fft_input.size() ):
    fft_input[ii] = complex( values[ii],0 )
  fft_obj = fftpack.complex_to_complex_2d( (np,np) )
  result = fft_obj.forward( fft_input )

  sp_image = flex.vec3_double()

  for ii in range(np):
    for jj in range(np):
      kk = ii
      ll = jj
      if kk > np/2: kk=kk-np
      if ll > np/2: ll=ll-np
      sp_image.append([kk+np/2,ll+np/2,abs( result[(ii,jj)] )**2.0 ] )

  
  Nq=np/2
  Nphi=Nq
  c2_image, normalized_sp = from_I_to_C2(sp_image, Nq, Nphi, N=np, further_reduct=True)

  image = normalized_sp
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  sp_n = N
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

  sp_moments = zernike_2d_mom.moments()

#  sp_reconst =  image_tools.generate_image2(nmax, sp_moments, normalized_sp, sp_n)
#  c2_image, normalized_sp = from_I_to_C2(sp_reconst, Nq, Nphi, N=np)

  sp_out = open('sp_image.dat', 'w')
  for pt in sp_image:
    print>>sp_out, (pt[0]-sp_n)*dq, (pt[1]-sp_n)*dq, pt[2]
  sp_out.close()


  c2_output=open('c2.dat', 'w')
  for pt in c2_image:
    print >> c2_output, pt[0],pt[1],pt[2]
  c2_output.close()


  image = c2_image
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  c2_n = N
#  grid_2d = math.two_d_grid(N, nmax)
#  grid_2d.clean_space( image )
#  grid_2d.construct_space_sum()
#  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

#  c2_moments = zernike_2d_mom.moments()

#  coefs = ( c2_moments )
  nls = math.nl_array(nmax).nl()
#  for nl, c2m in zip( nls, coefs ):
#    print nl, c2m


  sp_moments[0]=0
  nmax4c2 = 2*nmax
  Inm = math.nl_c_array(nmax4c2)
  coef_table=image_tools.calc_integral_triple_zernike2d(nmax4c2)
  Inm.load_coefs(nls, sp_moments)
  cnm = image_tools.calc_Cnm_from_Inm( Inm, coef_table, nmax4c2 )
  cnm_coef = cnm.coefs()

#  for nl, mm,mmm in zip( nls, c2_moments, cnm_coef ):
#    print abs(mm), abs(mmm)


  #cnm_coef = c2_moments
  c2_reconst = image_tools.generate_image2(nmax,cnm_coef, c2_image, c2_n)

  c2_r=open('c2r.dat', 'w')
  for pt in c2_reconst:
    print >> c2_r, pt[0],pt[1],pt[2]
  c2_r.close()

  exit()


if __name__ == "__main__":
  pdbfile=sys.argv[1]
  nmax=int(sys.argv[2])
  tst_fft_proj( pdbfile, nmax=nmax )
  #tst_projection( pdbfile )
  #tst_zernike_expansion(pdbfile)
