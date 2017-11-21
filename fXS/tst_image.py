from scitbx import math
from scitbx import differential_evolution as de
from scitbx import simplex
from scitbx import lbfgs
from scitbx import direct_search_simulated_annealing as dssa
from scitbx.array_family import flex
from scitbx.array_family import shared
from stdlib import math as smath
import random
import time, os, sys
from fractions import Fraction
from sastbx import fXS

def read_data(filename):
  file=open(filename, 'r')
  data=flex.vec3_double()
  for line in file:
    keys=line.split()
    if(len(keys)==3):
      x=float(keys[0])
      y=float(keys[1])
      z=float(keys[2])
      data.append([x,y,z])
  file.close()
  return data

def generate_image(n,moments, N=100):
  nmax=n

  image = flex.vec3_double()
  NP=2*N+1
  reconst=flex.complex_double(NP**2, 0)
  nl_array = math.nl_array( nmax )
  nls = nl_array.nl()

  r_theta = flex.vec2_double()

  for x in range(0,NP):
    x=x-N
    for y in range(0,NP):
      y=y-N
      rr = smath.sqrt(x*x+y*y)/N
      if rr>1.0:
        tt=0.0
      else:
        tt = smath.atan2(y,x)
      r_theta.append( [rr, tt] )


  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    #if(l/2*2 != l): continue
    if (l>0): c=c*2.0

    rap = math.zernike_2d_polynome(n,l)
    i=0
    for pt in r_theta:
      rr = pt[0]
      if rr>1.0:
        value=0.0
      else:
        tt = pt[1]
        value = rap.f(rr,tt)
      reconst[i]=reconst[i]+value*c
      i=i+1


  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      image.append([x,y,value])
      i=i+1
  return image

def generate_image2(n,moments, template_image, np):
  nmax=n

  image = flex.vec3_double()
  np_tot = template_image.size()
  reconst=flex.complex_double(np_tot, 0)
  nl_array = math.nl_array( nmax )
  nls = nl_array.nl()

  r_theta = flex.vec2_double()
  for pt in template_image:
    x=pt[0]-np
    y=pt[1]-np
    rr = smath.sqrt(x*x+y*y)/np
    if rr>1.0:
      tt=0.0
    else:
      tt = smath.atan2(y,x)
    r_theta.append( [rr, tt] )

  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    #if(l/2*2 != l): continue
    if (l>0): c=c*2.0
    rap = math.zernike_2d_polynome(n,l)
    i=0
    for pt in r_theta:
      rr = pt[0]
      if rr>1.0:
        value=0.0
      else:
        tt = pt[1]
        value = rap.f(rr,tt)
      reconst[i]=reconst[i]+value*c
      i=i+1

  for i in range(np_tot):
    x = template_image[i][0]
    y = template_image[i][1]
    value=reconst[i].real
    image.append([x,y, value])
  return image



def from_I_to_C2( nmax, I_image, Nq, Nphi, N=100 ):
  nx=N*2+1
  ny=nx
  image_cart = fXS.image_cartesian( nx, ny, I_image )
  image_cart.calc_c2_array( Nq, Nphi )
  c2_image=image_cart.get_c2_array()
  return c2_image

def tst_image(n, N=100, filename=None):
  nmax = n
  nmax0=8
  nl_array0 = math.nl_array( nmax0 )
  nls0 = nl_array0.nl()

  nl_array = math.nl_array( nmax )
  nls = nl_array.nl()
  if(filename is not None):
    image=read_data(filename)
  else:
    moments=flex.random_double(nls0.size() )
#    moments=flex.double(nls.size(),0 )
#    moments[3] = 1.0
#    moments[7] = 1.0
    image=generate_image(n,moments)
    orig=open('original.dat','w')
    for pt in image:
      print >>orig, pt[0],pt[1],pt[2]
    orig.close()

  Nq=100
  Nphi=100
  c2_image = from_I_to_C2(nmax, image, Nq, Nphi)

  c2_output=open('c2.dat', 'w')
  for pt in c2_image:
    print >> c2_output, pt[0],pt[1],pt[2]
  c2_output.close()

  coef_table=tst_integral_triple_zernike2d(nmax)
  moments[0]=0
  new_mom = moments.concatenate( flex.double(nls.size()-nls0.size(),0) )
  Inm = nl_array
  Inm.load_coefs(nls, new_mom)
  cnm = calc_Cnm_from_Inm( Inm, coef_table, nmax )
  cnm_coef = cnm.coefs()
  print "#cnm[0]", cnm_coef[0]
#  cnm_coef[0]=0

  ### calculate 2d zernike moments for C2_image ###
  image = c2_image
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

  c2_moments = zernike_2d_mom.moments()

  #coefs = flex.real( c2_moments )
  #cnm_coef  = flex.real( c2_moments )
  cnm_coef  = c2_moments
  c2_reconst = generate_image2(n,cnm_coef, c2_image,Nq)


  c2_r=open('c2r.dat', 'w')
  for pt in c2_reconst:
    print >> c2_r, pt[0],pt[1],pt[2]
  c2_r.close()

  ls_score = 0
  np_tot = c2_image.size()
  for p1, p2 in zip( c2_image, c2_reconst ):
    ls_score += (p1[2]-p2[2])**2.0

  print nmax, nls.size(), ls_score, np_tot*smath.log(ls_score/np_tot), "SUM"

  for nl, c2m in zip( nls, cnm_coef ):
    print nl, c2m

  exit() #

  ### calculate 2d zernike moments for C2_image ###
  image = c2_image
  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  tt2=time.time()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

  c2_moments = zernike_2d_mom.moments()

  #coefs = flex.real( c2_moments )
  coefs = ( c2_moments )
  for nl, c2m, c2mm in zip( nls, cnm_coef, coefs ):
    if( nl[0]/2*2 == nl[0]):
      print c2m, c2mm

  coefs = flex.real( moments )
  nl_array.load_coefs( nls, coefs )

  for nl, c in zip( nls, moments):
    if(abs(c)<1e-3):
      c=0
    print nl, c
  print

  reconst=flex.complex_double(NP**2, 0)
  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    if(l>0):
      c=c*2.0
    rap = math.zernike_2d_polynome(n,l)
    i=0
    for x in range(0,NP):
      x=x-N
      for y in range(0,NP):
        y=y-N
        rr = smath.sqrt(x*x+y*y)/N
        if rr>1.0:
          value=0.0
        else:
          tt = smath.atan2(y,x)
          value = rap.f(rr,tt)
        reconst[i]=reconst[i]+value*c
        i=i+1

  rebuilt=open('rebuilt.dat','w')
  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      print>>rebuilt, x,y,image[i][2],value
      i=i+1
  rebuilt.close()


def integrate_triple_zernike2d(n1,n2,n3,m, Bnmk_obj):
  value=Fraction(0)
  temp = long(0)
  ck = [long(0)]*(n1+n2+n3+1)
  for k1 in range(m,n1+1,2):
    for k2 in range(m,n2+1,2):
      for k3 in range(m,n3+1,2):
       # value = value+Bnmk_obj.get_coef(n1,m,k1)*Bnmk_obj.get_coef(n2,m,k2)*Bnmk_obj.get_coef(n3,m,k3)/(k1+k2+k3+2.0)
        temp = Bnmk_obj.get_coef(n1,m,k1)*Bnmk_obj.get_coef(n2,m,k2)*Bnmk_obj.get_coef(n3,m,k3)
        ck[k1+k2+k3] = ck[k1+k2+k3]+temp

  for kk in range(3*m,n1+n2+n3+1,2):
  #  print "%4d, %30d"%(kk,ck[kk])
    value = value + Fraction( ck[kk],(kk+2))
  return float(value)

class Bnmk (object):
  "Bnmk coefficient object hold 2d zernike expansion coefs"
  def __init__(self, nmax):
    self.nmax=nmax
    self.Bnmk=math.nmk_array(nmax)
    self.initialize_bnmk()

  def initialize_bnmk(self):
    for n in range(self.nmax, -1,-1):
      self.Bnmk.set_coef(n,n,n,1.0)
      for m in range(n-2,-1,-2):
        value = self.Bnmk.get_coef(n,m+2,n)*(n+m+2.0)/(n-m)
        self.Bnmk.set_coef(n,m,n,value)
        for k in range(n-2,m-1,-2):
          value = -self.Bnmk.get_coef(n,m,k+2)*(k+m+2.0)*(k+2.0-m)/(k+2.0+n)/(n-k)
          self.Bnmk.set_coef(n,m,k,value)

  def get_coef(self,n,m,k):
    return int(self.Bnmk.get_coef(n,m,k).real)

  def print_bnmk(self):
    for n in range(self.nmax+1):
      for m in range(n,-1,-2):
        for k in range(m,n+1,2):
          print n,m,k,self.get_coef(n,m,k)


def tst_integral_triple_zernike2d(nmax):
  Bnmk_obj = Bnmk(nmax)
  #Bnmk_obj.print_bnmk()
  coef_table = []

  for m in range(nmax+1):
    C_m_3n = math.nmk_array(nmax)
    for n1 in range(m,nmax+1,2):
      for n2 in range(m,n1+1,2):
        for n3 in range(m,n2+1,2):
          value = integrate_triple_zernike2d(n1,n2,n3,m,Bnmk_obj)
          C_m_3n.set_coef(n1,n2,n3,value)
        #  print m,n1,n2,n3,value.real
    coef_table.append( C_m_3n )
  return coef_table

def calc_Cnm_from_Inm( Inm, coef_table, nmax ):
  two_pi = smath.pi*2.0
  Cnm=math.nl_array(nmax)
  for n in range( 0,nmax+1,2 ): # only even number (n,m)
    for m in range(2,n+1,2 ):  # when m=0, the In0 is subtracted
      temp = 0
      for n1 in range( m,nmax+1,2 ):
        for n2 in range( m,n1,2 ):
          i,j,k=sorted([n,n1,n2],reverse=True)
          temp = temp + coef_table[m].get_coef(i,j,k).real*Inm.get_coef(n1,m)*Inm.get_coef(n2,m)

        i,j,k=sorted([n,n1,n1],reverse=True)
        temp = temp + coef_table[m].get_coef(i,j,k).real*Inm.get_coef(n1,m)**2.0/2.0
      Cnm.set_coef(n,m,temp*2.0*two_pi)
  return Cnm





if __name__ == "__main__":
  args = sys.argv[1:]
  t1 = time.time()
  args = sys.argv[1:]
  if( len(args) > 0 ):
    nmax = int( args[0] )
  else:
    nmax=10

  tst_image(nmax)


  t2 = time.time()
  print "time used: ", t2-t1
  print "OK"
