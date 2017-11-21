from scitbx.array_family import flex
from iotbx import pdb
from scitbx import math, fftpack
from sastbx.fXS import projection
from sastbx.fXS import image_tools
from sastbx import fXS
from stdlib import math as smath
import random
import os, sys


def analytical_c2(Nq,Nphi,nls,moments):
  image=flex.vec3_double()
  for xx in range( 0, Nq+1 ):
    for yy in range( 0, Nq+1 ):
      r=smath.sqrt(xx*xx+yy*yy)/Nq
      if(r <= 1):
        phi = smath.atan2(yy,xx)
        if(phi<0): phi = phi + two_pi
        c2 = evaluate_c2( nls, moments, r, phi )
      else:
        c2 = 0
      image.append( [Nq+xx,Nq+yy, c2] )
      if( xx>0 and yy>0):
        image.append( [Nq-xx,Nq-yy, c2] )
      if( yy>0):
        image.append( [Nq+xx,Nq-yy, c2] )
      if( xx>0):
        image.append( [Nq-xx,Nq+yy, c2] )
  return image


def evaluate_c2( nls, moments, r, phi ):
  c2=0
  two_pi = smath.pi*2.0
  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    if(l/2*2 != l): continue
    if (l>0): c=c*2.0
    rap = math.zernike_2d_polynome(n,l)
    for nl1,c1 in zip( nls, moments):
      n1=nl1[0]
      l1=nl1[1]
      if( l != l1 ): continue
      if (l1>0): c1=c1*2.0
      rap1 = math.zernike_2d_polynome(n1,l1)
      tmp_c2 = rap.f(r,phi/2)      # using phi/2 because it will add up later
      tmp_c2 = tmp_c2*rap1.f(r,phi/2)
      c2 = c2+tmp_c2.real*c*c1
  return c2*two_pi


def build_grid(Nq,Nphi,nls, moments):
  two_pi = smath.pi*2.0
  dq=1
  dphi=two_pi/Nphi
  qps = flex.vec2_double()
  xys = flex.vec2_double()
  image=flex.vec3_double()
  for xx in range( 0, Nq+1 ):
    for yy in range( 0, Nq+1 ):
      r=smath.sqrt(xx*xx+yy*yy)/Nq
      if(r <= 1):
        phi = smath.atan2(yy,xx)
        if(phi<0): phi = phi + two_pi
        qps.append( [r, phi] )
        xys.append( [xx,yy] )
        c2 = compute_c2( nls, moments, r, phi )
      else:
        c2 = 0
      image.append( [Nq+xx,Nq+yy, c2] )
      if( xx>0 and yy>0):
        image.append( [Nq-xx,Nq-yy, c2] )
      if( yy>0):
        image.append( [Nq+xx,Nq-yy, c2] )
      if( xx>0):
        image.append( [Nq-xx,Nq+yy, c2] )
  return image

def compute_c2( nls, moments, r, phi, Nsample=100 ):
  c2=0
  two_pi = smath.pi*2.0
  sp1 = flex.double(Nsample,0)
  sp2 = flex.double(Nsample,0)
  theta = flex.random_double(Nsample)*two_pi
  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    if(l/2*2 != l): continue
    if (l>0): c=c*2.0
    rap = math.zernike_2d_polynome(n,l)
    for ii in range(Nsample):
      value1 = rap.f(r,theta[ii]+phi).real*c
      value2 = rap.f(r,theta[ii]).real*c
      sp1[ii] = sp1[ii] + value1
      sp2[ii] = sp2[ii] + value2
  return flex.sum(sp1*sp2)/Nsample






def from_I_to_C2( I_image, Nq, Nphi, N=100, further_reduct=True, smoothen=True):
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


def tst_line(nmax=20, width=2,N=201):

  nl_array=math.nl_array(nmax)
  nls = nl_array.nl()
  moments = flex.random_double(nls.size())*0
  nl_array.load_coefs(nls, moments)
  nl_array.set_coef(8,6,1.0)
  nl_array.set_coef(6,2,1.0)
  nl_array.set_coef(4,4,1.0)
  moments = nl_array.coefs()
  this_nls = [[8,6],[6,2],[4,4]]
#  this_nls = [[6,2]]
  this_moments = [1.0,1.0,1.0]

  N=N
  half_N = N/2
  Nq=half_N
  Nphi=Nq

  ## The following two methods should give consistent results ##
  #c2_image = build_grid(Nq, Nphi, this_nls, this_moments )
  c2_image = analytical_c2(Nq, Nphi, this_nls, this_moments )

  c2_output=open('c2_raw.dat', 'w')
  for pt in c2_image:
    print >> c2_output, pt[0],pt[1],pt[2]
  c2_output.close()

  #sp_out = open('sp_image.dat', 'w')
  #for pt in raw_image:
  #  print>>sp_out, pt[0], pt[1], pt[2]
  #sp_out.close()

  ## get cnm from sp_moments ##

  sp_moments=flex.complex_double(moments, moments*0)
  #sp_moments[0]=0
  nmax4c2 = nmax*2
  Inm = math.nl_c_array(nmax4c2)
  coef_table=image_tools.calc_integral_triple_zernike2d(nmax4c2)
  Inm.load_coefs(nls, sp_moments)
  cnm = image_tools.calc_Cnm_from_Inm( Inm, coef_table, nmax4c2 )
  cnm_coef = cnm.coefs()

  c2_reconst = image_tools.generate_image2(nmax4c2,cnm_coef, c2_image, half_N)

  c2_output=open('c2_reconst.dat', 'w')
  for pt in c2_reconst:
    print >> c2_output, pt[0],pt[1],pt[2]
  c2_output.close()


  ## get zm for sp_image ##
  image = c2_image
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  sp_n = N
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )
  sp_moments = zernike_2d_mom.moments()

  print "#C2"
  for nl, m1,m2 in zip( nls,cnm_coef, sp_moments):
    print nl, m1, m2.real, m1/(m2.real+1e-23)



if __name__ == "__main__":
  #pdbfile=sys.argv[1]
  nmax=int(sys.argv[1])
  tst_line(nmax=nmax)
  #tst_fft_proj( pdbfile, nmax=nmax )
  #tst_projection( pdbfile )
  #tst_zernike_expansion(pdbfile)
