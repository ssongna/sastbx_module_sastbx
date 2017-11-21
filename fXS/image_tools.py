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
  # Generate images from the zernike moments, N is the number of points from 0 to 1
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


  return reconst

  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      image.append([x,y,value])
      i=i+1
  return image

def generate_image2(n,moments, template_image, np):
  # Generate images from the zernike moments, the grid indices are copied from template_image
  nmax=n
  image = flex.vec3_double()
  np_tot = template_image.size()
  reconst=flex.double(np_tot, 0)
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
      reconst[i]=reconst[i]+(value*c).real
      i=i+1

  for i in range(np_tot):
    x = template_image[i][0]
    y = template_image[i][1]
    value=reconst[i]
    image.append([x,y, value])
  return image


class image_generator(object):
  def __init__(self, base_image,nmax):
    self.bi = base_image.deep_copy()
    self.np = self.bi.focus()[0]
    self.nmax = nmax
    self.n = int(self.np/2)
    print self.bi.focus(), self.np, self.n
    print "Building mom->img grid"
    t1 = time.time()
    self.z2dg = math.zernike_grid_2d(self.n, self.nmax)
    t2 = time.time()
    print "time used: %d (seconds)"%(t2-t1)
    print "Build img->mom grid"
    self.grid_2d = math.two_d_grid(self.n, self.nmax)

  def get_coefs(self,image):
    # get moments
    nimg = flex.vec3_double()
    for ii in range(self.np):
      for jj in range(self.np):
        nimg.append( (ii,jj,image[(ii,jj)]) )
    ori_image = nimg.deep_copy()
    self.grid_2d.clean_space( nimg )
    self.grid_2d.construct_space_sum()
    zernike_2d_mom = math.two_d_zernike_moments( self.grid_2d, self.nmax )
    return zernike_2d_mom.nm(), zernike_2d_mom.moments()

  def build_image(self, nm,cnm):
    print nm, cnm
    self.z2dg.load_coefs(nm,cnm)
    f = flex.real(self.z2dg.f())
    f.reshape( flex.grid(self.np,self.np) )
    return f

def dummy(n1,n2,n3,m):
  return 0.0

def num_int_triple_int(n1,n2,n3,m, np=1e5):
  r = flex.double( range(int(np+1)) )/np
  f1 = math.zernike_2d_radial(n1,m).f_array(r)
  f2 = math.zernike_2d_radial(n2,m).f_array(r)
  f3 = math.zernike_2d_radial(n3,m).f_array(r)
  result = flex.sum(f1*f2*f3*r)/np
  return result

class fast_tic(object):
  def __init__(self,nmax):
    self.nmax = nmax
    self.Bnmk_obj = Bnmk(nmax)

  def ennnm(self,n1,n2,n3,m):
    value = integrate_triple_zernike2d(n1,n2,n3,m,self.Bnmk_obj)
    return value

class triple_zernike_integral_manager(object):
  def __init__(self,nmax, method=None):
    self.nmax = nmax
    self.nnnm = {}
    self.coefs = flex.double()
    self.method = method
    if self.method is None:
      self.method = dummy

    count=0
    for n1 in range(self.nmax+1):
      for n2 in range(n1,self.nmax+1):
        if self.parity_check(n1,n2):
          for n3 in range(n2,self.nmax+1):
            if self.parity_check(n2,n3):
              no1,no2,no3 = self.sort_n( n1,n2,n3 )
              for mm in range(no1+1):
                if self.parity_check(n1,mm):
                  #print no1,no2,no3, mm,
                  self.nnnm[ (no1,no2,no3,mm) ] = count
                  tmp = self.get_integral(no1,no2,no3,mm)
                  self.coefs.append( tmp )
                  #print tmp
                  count += 1

  def parity_check(self, a,b):
    d = abs(a-b)
    dd = 2*int(d/2)
    if dd==d:
      return True
    return False


  def get_integral(self, n1,n2,n3,m):
     return self.method( n1,n2,n3,m)

  def sort_n(self,n1,n2,n3):
    # sort low to high
    n = flex.int( [n1,n2,n3] )
    o = flex.sort_permutation( n )
    return n[o[0]],n[o[1]],n[o[2]]





def tst_img_gen():
  nmax=20
  ft = fast_tic(nmax)
  tzim = triple_zernike_integral_manager(nmax,ft.ennnm)
  print "(0,0)",tzim.get_integral(0,4,4,0)
  print "(2,0)",tzim.get_integral(2,4,4,0)
  print "(4,0)",tzim.get_integral(4,4,4,0)*5
  print "(6,0)",tzim.get_integral(4,4,6,0)
  print "(8,0)",tzim.get_integral(4,4,8,0)*9
  print "(10,0)",tzim.get_integral(4,4,10,0)
  #tzim = triple_zernike_integral_manager(5,num_int_triple_int)

  print "Building ball"
  import basic_shapes
  rod = basic_shapes.ball(261,nmax)
  ig = image_generator(rod,nmax)
  print "Buiding IG"
  nm,cnm = ig.get_coefs( rod )

  # test single coefficient
  this_nm = (4,0)
  cnm = cnm*0.0
  for iii,nnn in enumerate(nm):
    if nnn[0]==this_nm[0]:
      if nnn[1]==this_nm[1]:
        cnm[iii]=1.0
  f = ig.build_image(nm,cnm)
  import zernike_phaser
  io = zernike_phaser.integration_object(261)
  #c2 = io.build_c2_function( f )
  c2 = f*f
  zernike_phaser.write_generic_image(c2, "sp_c2.dat")
  zernike_phaser.write_generic_image(f, "sp.dat")

  nm,cnm = ig.get_coefs( c2 )
  for n,nn in zip(nm,cnm):
    if(abs(nn)>1e-6):
      print n,nn #*130**2.0
  g = ig.build_image(nm,cnm)
  cnm = cnm*0.0
  for iii,nnn in enumerate(nm):
    if (nnn[1]==0):
      if nnn[0]==0:
        cnm[iii]=0.1
      if (nnn[0]==4): 
        cnm[iii]=0.02857*5.0
      if (nnn[0]==8): 
        cnm[iii]=0.02857*9.0
  zernike_phaser.write_generic_image(g, "sp_c2_r.dat")

  return
  nm,nnm = ig.get_coefs( f )
  g = ig.build_image(nm,nnm)
  h = generate_image(10,cnm,100)
  h = flex.double(flex.real(h))
  h.reshape( flex.grid(201,201) )

  for c,i,j in zip(nm,cnm,nnm):
    print "#", c,i,j

  for ii in range(ig.np):
    for jj in range(ig.np):
      print ii,jj, f[(ii,jj)], g[(ii,jj)], h[ (ii,jj) ]
    print


  return
  #moms = ig.get_coefs( rod )
  f = ig.build_image(nm,cnm)
  for ii in range(ig.np):
    for jj in range(ig.np):
      print ii,jj, rod[ (ii,jj) ], f[(ii,jj)]
    print
  return


# from intensity image, the C2 is calculated
def from_I_to_C2( nmax, I_image, Nq, Nphi, N=100 ):
  nx=N*2+1
  ny=nx
  image_cart = fXS.image_cartesian( nx, ny, I_image )
  further_reduct=True
  image_cart.calc_c2_array( Nq, Nphi,further_reduct )
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

  coef_table=calc_integral_triple_zernike2d(nmax)
  moments[0]=0
  new_mom = moments.concatenate( flex.double(nls.size()-nls0.size(),0) )
  nmax4c2 = 2*nmax
  Inm = math.nl_c_array(nmax4c2)

  Inm.load_coefs(nls, new_mom)
  cnm = calc_Cnm_from_Inm( Inm, coef_table, nmax4c2 )
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


# integration of R_n1^m R_n2^m R_n3^m * r, from [0,1]
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
#    print "%4d, %30d"%(kk,ck[kk])
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
    return int(self.Bnmk.get_coef(n,m,k))

  def print_bnmk(self):
    for n in range(self.nmax+1):
      for m in range(n,-1,-2):
        for k in range(m,n+1,2):
          print n,m,k,self.get_coef(n,m,k)

# calculate triple integrals for all possible (n1,n2,n3,m) up to nmax
# the coefficients En1,n2,n3,m are saved in coef_table



def calc_integral_triple_zernike2d(nmax):
  Bnmk_obj = Bnmk(nmax)
  Bnmk_obj.print_bnmk()
  coef_table = []

  for m in range(nmax+1):
    C_m_3n = math.nmk_array(nmax)
    for n1 in range(m,nmax+1,2):
      for n2 in range(m,n1+1,2):
        for n3 in range(m,n2+1,2):
          value = integrate_triple_zernike2d(n1,n2,n3,m,Bnmk_obj)
          C_m_3n.set_coef(n1,n2,n3,value)
          #print n1,n2,n3, m, C_m_3n.get_coef(n1,n2,n3), value
    coef_table.append( C_m_3n )
  return coef_table

# Given Inm (expansion coef from intensity), Cnm can be calculated
# using formula involving E(n1,n2,n3,m)
def calc_Cnm_from_Inm( Inm, coef_table, nmax ):
  two_pi = smath.pi*2.0
  four_pi= 2.0*two_pi
  Cnm=math.nl_array(nmax)
  for n in range( 2,nmax+1,2 ): # only even number (n,m)
    for m in range(2,n+1,2 ):  # when m=0, the In0 is subtracted
      temp = 0
      for n1 in range( m,nmax+1,2 ):
        for n2 in range( m,nmax+1,2 ):
          i,j,k=sorted([n,n1,n2],reverse=True)
          temp = temp + coef_table[m].get_coef(i,j,k)*(Inm.get_coef(n1,m).conjugate()*Inm.get_coef(n2,m)).real

#        i,j,k=sorted([n,n1,n1],reverse=True)
#        temp = temp + coef_table[m].get_coef(i,j,k)*abs(Inm.get_coef(n1,m))**2.0/2.0
      Cnm.set_coef(n,m,temp*four_pi*(n+1))
  return Cnm



def tst_integral_triple_zernike2d(nmax):
  # E^m_(n1,n2,n3)
  # (0;0,0,0) = 1/2
  # (1;1,1,1) = 1/5
  # (0;2,0,0) = 0
  # (0;2,2,0) = 1/6
  # (2;2,2,2) = 1/8
  # (2;4,2,2) = 1/40
  # (2;4,4,2) = 7/120
  # (2;4,4,4) = -1/280
  # (4;4,4,4) = 1/14
  # (6;6,6,6) = 1/20

  coefs_table = calc_integral_triple_zernike2d(nmax)
  assert abs( coefs_table[0].get_coef(0,0,0) - 1/2.0 ) < 1e-3
  assert abs( coefs_table[1].get_coef(1,1,1) - 1/5.0 ) < 1e-3
  assert abs( coefs_table[0].get_coef(2,0,0) - 0.0 ) < 1e-3
  #print coefs_table[0].get_coef(2,2,0), 1/6.0
  assert abs( coefs_table[0].get_coef(2,2,0) - 1/6.0 ) < 1e-3
  assert abs( coefs_table[2].get_coef(2,2,2) - 1/8.0 ) < 1e-3
  assert abs( coefs_table[2].get_coef(4,2,2) - 1/40.0 ) < 1e-3
  assert abs( coefs_table[2].get_coef(4,4,2) - 7/120.0 ) < 1e-3
  assert abs( coefs_table[2].get_coef(4,4,4) + 1/280.0 ) < 1e-3
  assert abs( coefs_table[4].get_coef(4,4,4) - 1/14.0 ) < 1e-3
  assert abs( coefs_table[6].get_coef(6,6,6) - 1/20.0 ) < 1e-3

  scale=3.14*4
  for ii in range(2,6+1,2):
    print coefs_table[2].get_coef(6,6,ii)*scale*(ii+1)
  for ii in range(8,nmax+1,2):
    print coefs_table[2].get_coef(ii,6,6)*scale*(ii+1)

def tst_Inm_array(nmax):
  Inm = math.nl_c_array(nmax)
  nls = Inm.nl()
  moments = flex.double(range(nls.size()) )
  moments = flex.complex_double( moments, moments)
  Inm.load_coefs(nls, moments)
  for ii, nl in zip( range( nls.size()), nls):
    assert(Inm.get_coef(nl[0],nl[1]).real==ii)

if __name__ == "__main__":
  #nmax=12
  #tst_integral_triple_zernike2d(nmax)
  #tst_Inm_array(nmax)
  tst_img_gen()
  print "OK"
