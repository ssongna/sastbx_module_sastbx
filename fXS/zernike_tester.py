from sastbx.fXS import projection, image_tools
from scitbx.array_family import flex
from iotbx import pdb
from scitbx import math, fftpack
import os, sys
from stdlib import math as smath
from libtbx.utils import n_dim_index_from_one_dim 
from scitbx.python_utils import random_transform as rt



def write_image(data,outobj):
  for pp in data:
    print >> outobj, pp[0],pp[1],pp[2]

def write_generic_image(data,filename):
    pp = open(filename,'w')
    is_list = False
    m = None
    np=None
    if type(data) == type([]):
      m = len(data)
      np,np = data[0].focus()
    else:
      np,np = data.focus()
  
    for ii in range(np):
      for jj in range(np):
          print >> pp, ii, jj, 
          if m is None:
            print >>pp, data[ (ii,jj) ]
          else:
            for mm in range(m):
              print >> pp, data[mm][ (ii,jj) ],
            print >> pp
      print >> pp

def pad_around(data,n):
  m,m=data.focus()
  d=n+m+n
  new_data = flex.double( flex.grid(d,d), 0.0 )
  for ii in range(m):
    for jj in range(m):
      new_data[ (n+ii,n+jj) ] = data[ (ii,jj) ]
  return new_data
  

def enlarge(data,factor,sigma=0.1,full=False):
  n,n = data.focus()
  m = int(n*factor)
  x    = flex.double()
  y    = flex.double()
  vals = flex.double()
  sigma=sigma*factor 
  new_data = flex.double( flex.grid(m,m), -9 )
  visited = flex.bool( flex.grid(m,m), False )
  oo = n/2.0
  for ii in range(n):
    for jj in range(n):
      dd = smath.sqrt( (ii-oo)**2.0 + (jj-oo)**2.0 )
      if dd <= oo:
        nx = ii*factor
        ny = jj*factor
        x.append( nx )
        y.append( ny )
        vals.append( data[ (ii,jj) ] )
        new_data[ (int(nx), int(ny)) ] = data[ (ii,jj) ]
        if not full:
          visited[ (int(nx), int(ny)) ] = True

  
  not_visited = ~visited
  not_visited = not_visited.iselection()
  # now we need to loop over all non-visited pixel values
  for pixel in not_visited:
        nv = -9
        index =  n_dim_index_from_one_dim(pixel, [m,m] )
        nvx  = index[1]
        nvy  = index[0]
        dx   = x-nvx
        dy   = y-nvy
        dd   = (dx*dx+dy*dy)
        ss   = flex.exp( -dd/(sigma*sigma) )
        nv   = flex.sum(ss*vals)/(1e-12+flex.sum(ss))
        new_data[ (nvx,nvy) ] = nv
        visited[  (nvx,nvy) ] = True
        #print nvx, nvy, nv
  not_visited = ~visited
  not_visited = not_visited.iselection()
  #print not_visited.size()
  return new_data
  oo=m/2.0
  for ii in range(m):
    for jj in range(m):
      dd = smath.sqrt( (ii-oo)**2.0 + (jj-oo)**2.0 )
      new_data[ (ii,jj) ] = new_data[ (ii,jj) ]/(1+smath.sqrt(dd))
  return new_data
      

       




def cut_data(data,np,max_radius):
   mid_point = int(np/2)
   nii = mid_point-max_radius
   njj = mid_point-max_radius
   data = data.matrix_copy_block(nii,njj,max_radius*2+1,max_radius*2+1)
   # set things to zero for r values above radius
   n,n = data.focus()
   for ii in range(n):
     for jj in range(n):
       dx = ii-n/2.0
       dy = jj - n/2.0
       rr = smath.sqrt(dx*dx + dy*dy) 
       if rr >= max_radius:
         data[ (ii,jj) ]= 1e-12
   return data

def reorder_and_cut_data(data,np,max_radius=None):
   if max_radius is None:
     max_radius = np/4
   new_data = flex.double(flex.grid(np,np),0.0 )
   for ii in range(np):
      for jj in range(np):
         nope=False
         kk = ii
         ll = jj
         if kk > np/2: kk=kk-np
         if ll > np/2: ll=ll-np
         
         #radius is now needed
         rr = (kk)**2.0 + (ll)**2.0
         if rr > max_radius*max_radius:
           new_val= 1e-15
         else:
           new_val = float( data[(ii,jj)] )
         kok = int(kk+np/2)
         kuk = int(ll+np/2)
         if kok < 0: kok = 0
         if kok > np-1 : kok = np-1
         if kuk < 0: kuk = 0
         if kuk > np-1 : kuk = np-1
         new_data[ ( kok,kuk ) ] = new_val

   #now we need to redefine the left upper corner
   new_data = cut_data(new_data,np,max_radius)
   return new_data

def write_scattering_pattern(data,np,outobj):
   for ii in range(np):
      for jj in range(np):
         print >> outobj, ii,jj,smath.log(data[(ii,jj)] )
      print >> outobj


class integration_object(object):
  def __init__(self,np):
    self.np = np
    self.grid = flex.grid( np,np)
    self.mask = flex.double(self.grid,1.0)
    self.r_image = flex.double(self.grid,0)
    self.t_image = flex.double(self.grid,0)
    self.saxs_image = flex.double(self.grid,0)
    self.ox, self.oy = np/2.0, np/2.0

    self.build_mask()
    self.build_images()
    self.build_saxs_template(self.np)
    self.c2_pixels, self.c2_pixels_angle_bin_nos, self.c2_angle_bins = self.build_angle_list(15)
    self.n_sets = 0.0

  def build_mask(self, rat=0.95):
    max_d = rat*min(self.np-self.ox, self.np-self.oy)
    max_d = (max_d**2.0)*rat
    new_mask = flex.double(self.grid,1.0)
    # we only want / need the circular part of the image
    for ii in range(self.np):
      for jj in range(self.np):
        dd = (ii-self.ox)**2 + (jj-self.oy)**2
        if dd > max_d:
          new_mask[ (ii,jj) ] = 0
    self.mask = new_mask

  def build_images(self):
    for ii in range(self.np):
      for jj in range(self.np):
        rr = smath.sqrt( (ii+0.5-self.ox)**2.0 + (jj+0.5-self.oy)**2.0 )
        self.r_image[ (ii,jj) ] = rr
        tt = smath.atan2( jj+0.5-self.oy,ii+0.5-self.ox )
        self.t_image[ (ii,jj) ] = tt

  def build_saxs_template(self, n_slots=20):
    print "Build SAXS template"
    # we need to get mean values for the intensity and build an image of that
    min_dist,max_dist = flex.min( self.r_image ), flex.max( self.r_image )
    bin_step = (max_dist-min_dist)/(n_slots-1)

    # prebin
    self.saxs_bin_image = flex.int( self.grid, -9 ) 
    self.bin_array = []
    self.r_array = []
    for ii in range(n_slots):
      self.bin_array.append( [] )
      self.r_array.append( ii )

    # make a binned image
    for ii,dd in enumerate(self.r_image):
      this_d = dd
      this_ii = int(this_d/bin_step+0.5)
      self.saxs_bin_image[ii] = this_ii
      if self.mask[ii] > 0.1:
        self.bin_array[this_ii].append( ii )       

  def build_angle_list(self, bins=100):
    # the angles we are after will lie inbetween 0 and 90 ( smath.pi/2.0)
    # the rest we can generate via symmetry
    these_pixels = []
    these_angles = flex.double()
    these_bin_nos = []
    angle_bins = flex.double( range(bins) )*smath.pi*0.5/(bins-1)
    da = smath.pi/2.0
    da = da/(bins-1)
    self.angle_bin_image = flex.int( self.grid,-9 )
    for ii in range(self.np):
      for jj in range(self.np):
        theta = self.nice_angle( self.t_image[ (ii,jj) ] )
        ntheta = self.angle_syms(theta)
        self.angle_bin_image[ (ii,jj) ] = int(ntheta/da+0.5)
        if theta <= smath.pi/2.0:
          # we want to keep this
          these_pixels.append( (ii,jj) )
          these_angles.append(  theta  )
          these_bin_nos.append( int(theta/da+0.5) )
    return these_pixels, these_bin_nos, angle_bins          


  def c2_at_angle(self, data, angle):
    # take the data and rotate it
    rot_data = self.rotate_image(data,angle*180.0/smath.pi)
    # multiply the lot
    multi = rot_data*data
    # now we need a saxs curve
    c2p = self.get_saxs_data(multi)    
    return c2p

  def angle_syms(self,x, return_all=False):
    x0 = self.nice_angle(x)
    x1 = self.nice_angle(smath.pi+x)
    x2 = self.nice_angle(-x)
    x3 = self.nice_angle(smath.pi-x)
    angs = flex.double( [x0,x1,x2,x3] )
    if return_all:
      return angs
    else:
      return flex.min(angs)


  def build_c2_function(self,data):
    partial_c2_curves = []
    for phi in self.c2_angle_bins:
      print phi
      c2p = self.c2_at_angle(data,phi)
      partial_c2_curves.append( c2p )

    # we have tabluated all entries for each possible angle and each possible r value, please fill in the image
    c2 = flex.double(self.grid,0)

    for ix in range(self.np):
      for iy in range(self.np):
        this_r = self.r_image[ (ix,iy) ]
        this_t = self.t_image[ (ix,iy) ]
        # map this t back to 0/90 domain
        this_t = self.angle_syms(this_t)
        
        # please get back the r bin number
        r_bin = self.saxs_bin_image[ (ix,iy) ] 
        # please get the angle bin id
        t_bin = self.angle_bin_image[ (ix,iy) ]
        # retrieve entry please
        val = partial_c2_curves[t_bin][r_bin]
        c2[ (ix,iy) ] = val        

    self.write_generic_image( c2, 'c2.dat' )
    return c2

  def nice_angle(self,x):
    while x < 0:
      x = x+smath.pi*2.0
    x=x%(smath.pi*2.0)
    return x

  def get_saxs_data(self,data):
    saxs_data = flex.double( len(self.bin_array) ,0.0 )
    for ii,this_array in enumerate(self.bin_array):
      for jj in this_array:
        saxs_data[ii]+=data[jj]
    for ii,ss in enumerate(saxs_data):
      saxs_data[ii] = (ss+1e-12)/(1e-12+len(self.bin_array[ii]))
    return saxs_data

 
  def accumulate_saxs_data(self, data):
    self.n_sets += 1.0 
    #first get the saxs data
    saxs_data = self.get_saxs_data(data)
    # please loop over the sax image and fill in values
    for ii in range(self.np**2):
      self.saxs_image[ii] += saxs_data[ self.saxs_bin_image[ii] ] 

    self.saxs_image = self.saxs_image*self.mask 
   

     

  def rotate_image(self,data,angle_deg,sigma=1.0):
    # we take an image, rotate, interpolate it back onto a regular grid and return it
    #
    # There is however the problem that some pixels will never be seen. those pixels will be filled up by interpolation
    data = data*self.mask 
    cost = smath.cos( angle_deg*smath.pi/180.0 )
    sint = smath.sin( angle_deg*smath.pi/180.0 )
    new_image = flex.double(self.grid,0)
    visited   = flex.bool(self.grid,False)
    nx_array  = flex.double()
    ny_array  = flex.double()
    val_array = flex.double()
    for ix in range(self.np):
      for iy in range(self.np):
        x  = ix - self.ox
        y  = iy - self.oy
        nx = cost*x - sint*y
        ny = sint*x + cost*y
        jx = int(nx+self.ox+0.5)
        jy = int(ny+self.oy+0.5)
        if jx >= self.np : jx = self.np-1
        if jy >= self.np : jy = self.np-1
        if jx<0: jx=0
        if jy<0: jy=0
        new_image[ (jx,jy) ] = data[ (ix,iy) ]
        visited[ (jx,jy) ] = True
        nx_array.append( nx+self.ox )
        ny_array.append( ny+self.oy )
        val_array.append( data[ (ix,iy) ] )
    
    assert nx_array.size() == self.mask.size()
    not_visited = ~visited
    not_visited = not_visited.iselection()
    # now we need to loop over all non-visited pixel values
    for pixel in not_visited:
      #for ix in range(self.np):
      # for iy in range(self.np):
      nv = 0
      if self.mask[(ix,iy)]>-0.1:
        nv = -9
        index =  n_dim_index_from_one_dim(pixel, [self.np,self.np] )
        nvx  = index[1]
        nvy  = index[0]
        dx   = nx_array-nvx
        dy   = ny_array-nvy
        dd   = (dx*dx+dy*dy)
        ss   = flex.exp( -dd/(sigma*sigma) )*self.mask
        nv   = flex.sum(ss*val_array)/(1e-12+flex.sum(ss))
        new_image[ (ix,iy) ] = nv
    new_image = new_image*self.mask
    return new_image
   
    
  def write_saxs_image(self, file_name):
    pp= open(file_name,'w')
    for ii in range(self.np):
      for jj in range(self.np):
        print >> pp, ii, jj, smath.log( self.saxs_image[ (ii,jj) ]+1e-12 )
      print >> pp

  def write_generic_image(self,data,filename):
    pp = open(filename,'w')
    for ii in range(self.np):
      for jj in range(self.np):
        if self.mask[ (ii,jj) ] > 0.1:
          print >> pp, ii, jj, data[ (ii,jj) ] 
        else:
          print >> pp, ii, jj, 1e-6
      print >> pp

  def build_relative_fluctuations(self,data):
    tmp = (data)/ (1e-12+self.saxs_image)
    self.write_generic_image(tmp,'flucts.dat')
    return tmp

  def single_correlation(self,data):
    tmp = (data)/ (1e-12+self.saxs_image)
    c2 = self.accumulate_correlation_data(tmp)
    self.write_generic_image(c2,'c2.dat')
 


def build_hook(np):
  img = flex.double( flex.grid(2*np+1,2*np+1), 0 )
  for ii in range(2*np+1):
    img[ (np+0,ii) ] = 1.0
    img[ (np-1,ii) ] = 1.0
    img[ (np+1,ii) ] = 1.0

    img[ (ii,np+0) ] = 1.0
    img[ (ii,np-1) ] = 1.0
    img[ (ii,np+1) ] = 1.0
  return img


def zernike_expansion(image, nmax=30):
  ori_image=image.deep_copy()
  output=open('ori_c2.dat', 'w')

  # get moments
  NP=int(smath.sqrt( image.size() ))
  nimg = flex.vec3_double()
  ori_image = nimg.deep_copy()
  for ii in range(NP):
    for jj in range(NP):
      nimg.append( (ii,jj,image[(ii,jj)]) )
  ori_image = nimg.deep_copy()
  N=NP/2
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( nimg )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )
  # use moments to reconstruct image
  reconst = image_tools.generate_image2(nmax,zernike_2d_mom.moments(),ori_image,N)
  for pp,qq, in zip( ori_image, reconst):
    print pp[0], pp[1], pp[2], qq[2]

def add_error(data, fraction=0.1):
  n = data.size()
  e = rt.normal_variate(mu=0.0,sigma=fraction,N=n)
  e = e*data
  data = data+e
  return data



def build_scat_pat(data):
  print "DO FFT"
  np,np = data.focus()
  flex_grid=flex.grid(np,np)
  fft_input=flex.complex_double(flex_grid)
  for ii in range( fft_input.size() ):
    fft_input[ii] = complex( data[ii],0 )
  fft_obj = fftpack.complex_to_complex_2d( (np,np) )
  result = fft_obj.forward( fft_input )
  result = flex.abs( result )
  result = result*result
  result = reorder_and_cut_data(result,np,np/10)
  return result


class scat_pat_lib(object):
  def __init__(self, data, N=1000,fraction=0.15):
    print "Building library of rotated scattering patterns"
    self.data = data
    self.np,self.np = self.data.focus()
    self.N = N
    self.angles = flex.random_double(N)*smath.pi
    self.rot_img = []
    # build an object that allows image rotation
    self.io = integration_object(self.np)
    for ang in self.angles:
      print "Rotating Image"
      nimg = self.io.rotate_image( self.data, ang, sigma=0.5 )
      nimg = add_error(nimg,fraction)
      self.rot_img.append( nimg )

  def random_image(self):
    ind = flex.sort_permutation( flex.random_double(self.N) )[0]
    return self.rot_img[ ind ]


class experiment_gen(object):
  def __init__(self, libs):
    self.libs = libs

  def generate_random_image(self, stochoimetries,n_part=1):
    result = self.libs[0].random_image().deep_copy()*0.0
    parts = stochoimetries*n_part
    for this_component,part in enumerate(parts):
      part = int(part)
      for ii in range(part):
        result += self.libs[this_component].random_image()
    return result

class mixture_prep(object):
  def __init__(self, scat_lib, dro=None,  Ntot=50, navg=10):
    self.scat_lib = scat_lib
    self.Ntot=Ntot
    self.dro = dro
    self.stored_series = []
    self.navg = navg
    self.n = 0


  def do_all_series(self):
    for ii in range(self.navg):
      self.do_series()
    for ii in range(len(self.stored_series)):
      self.stored_series[ii]=self.stored_series[ii]/self.navg

    # do an svd analyses to find out what is going on
    mean = self.stored_series[0]*0.0
    new_data = flex.double()
    n = len(self.stored_series)
    for ii in range(n):
      mean = mean+self.stored_series[ii]
    mean = mean/n
    for ii in range(n):
      new_data.extend( (self.stored_series[ii]-mean).as_1d() )
    new_data.reshape( self.grid(n,mean.size() ) )
    
    from scitbx.linalg import svd
    svd_obj = svd.real( new_data.deep_copy(),True,True )
    for ii,ss in enumerate(svd_obj.sigma):
      print ii,ss, ss/flex.sum(svd_obj.sigma) , "MEAN EDUCED EV's"




    
    # get accurate fractions / kinetic data
    self.m = flex.double()
    for jj in range(len(self.stored_series)):
      print self.vals[jj], 'series vals'
      for oo in self.vals[jj]:
        self.m.append(oo/self.Ntot)
    self.m.reshape( flex.grid(len(self.stored_series),3) )

  def do_series(self):
    self.stochies = []
    self.start=1000 
    self.vals = [ [0,0,1000],[0,1000,0], [900,200,100] ]
    # build appropriate mixtures please
    for ii in range(3,100):
      val = [ self.start-ii*5, ii*2*5, ii*5 ]
      self.vals.append( val )
    self.mix = []
    for ijk in range(100):
      c1 = self.vals[ijk][0]
      c2 = self.vals[ijk][1]
      c3 = self.vals[ijk][2]
      print c1,c2,c3
      mixx = flex.int()
      for ii in range(c1):
        mixx.append(0)
      for ii in range(c2):
        mixx.append(1)
      for ii in range(c3):
        mixx.append(2)

      permut = flex.sort_permutation( flex.random_double(mixx.size()) )
      tm = mixx.select(permut)[0:self.Ntot]
      stch = []
      stch.append( tm.count(0) )
      stch.append( tm.count(1) )
      stch.append( tm.count(2) )
      self.mix.append( stch )
      print ijk, list(stch)

    self.images = []
    self.c2_images = []
    for stc in self.vals[:100]:
      img = self.scat_lib.generate_random_image( stc )
      self.images.append(img)
      c2_img = self.dro.build_c2_function(img)
      self.c2_images.append(c2_img )

    if len(self.stored_series)==0:
      self.stored_series = self.c2_images
    else:
      for ii in range(len(self.stored_series)):
        self.stored_series[ii] += self.c2_images[ii]
         




  def do_svd(self):
    # make data array
    data = flex.double()
    for set in self.stored_series: 
      data.extend( set.as_1d() )

    cols = self.stored_series[0].size()
    rows = len(self.stored_series) 

    data.reshape( flex.grid(rows,cols) )
    from scitbx.linalg import svd
    svd_obj = svd.real( data.deep_copy(),True,True )
    for ii,ss in enumerate(svd_obj.sigma):
      print ii,ss, ss/flex.sum(svd_obj.sigma) 

   
    from scitbx.linalg import svd
    svd_obj = svd.real( self.m.deep_copy(),True,True )
    m_pseudo_inverse = svd_obj.pseudo_inverse()
    
    cols = len(self.stored_series)
    rows = 3
    print rows,cols, m_pseudo_inverse.focus() 
    deconvol = []
    for ith in range(rows):
      result = self.stored_series[0]*0.0
      for jset in range(cols):
        result += m_pseudo_inverse[ (ith,jset) ]*self.stored_series[jset]
      deconvol.append( result )
    write_generic_image(deconvol,'deconvol_series.data')

    write_generic_image(self.stored_series, 'series.dat')




def tst_shapes():
  N=80
  length = 70
  width=10
  radius = 15
  from sastbx.fXS import basic_shapes
  rod = basic_shapes.rod(N,length,width,1.5)
  ball = basic_shapes.ball(N,radius,1.0)
  sball = basic_shapes.single_ball(N,radius,1.0)
  combo = rod+ball
  print "Building Rod"
  rod   = pad_around( enlarge(rod,1,sigma=3.0,full=True),   N  )
  print "Building Balls" 
  ball  = pad_around( enlarge(sball,1,sigma=3.0,full=True),  N  )
  print "Building Combo"
  combo = pad_around( enlarge(combo,1,sigma=3.0,full=True), N  )

  
  write_generic_image(rod,'rod.dat')
  write_generic_image(ball,'ball.dat')
  write_generic_image(combo,'combo.dat')

  frod   = build_scat_pat(rod) 
  fball  = build_scat_pat(ball)
  fcombo = build_scat_pat(combo)

  rod_lib   =  scat_pat_lib(frod,N=45)
  ball_lib  =  scat_pat_lib(fball,N=45)
  combo_lib =  scat_pat_lib(fcombo,N=45)


  #lets build c2 functiosn for each of these buggers
  np,np = fcombo.focus() 
  int_obj   = integration_object( np )
  c2_fball  = int_obj.build_c2_function(fball)
  c2_frod   = int_obj.build_c2_function(frod)
  c2_fcombo = int_obj.build_c2_function(fcombo)
  write_generic_image(c2_fball,  "c2_ball.dat")
  write_generic_image(c2_frod,   "c2_rod.dat")
  write_generic_image(c2_fcombo, "c2_combo.dat")


  t0_s = flex.double( [4, 5, 5] )
  t1_s = flex.double( [1, 11, 7] )
  t2_s = flex.double( [0, 9,  0] )
  t3_s = flex.double( [0,  0, 9] )

  eg = experiment_gen( [combo_lib,ball_lib,rod_lib] )

  mp = mixture_prep(eg, int_obj)
  mp.do_all_series()
  mp.do_svd()


  sys.exit()

  # build synthetic images
  t0_c2 =  c2_fcombo*t0_s[0] + c2_fball*t0_s[1] + c2_frod*t0_s[2]
  t1_c2 =  c2_fcombo*t1_s[0] + c2_fball*t1_s[1] + c2_frod*t1_s[2]
  t2_c2 =  c2_fcombo*t2_s[0] + c2_fball*t2_s[1] + c2_frod*t2_s[2]
  t3_c2 =  c2_fcombo*t3_s[0] + c2_fball*t3_s[1] + c2_frod*t3_s[2]
  a = flex.double( [4, 5, 5,1, 11, 7,0, 9,  0,0,  0, 9] )

  




  a.reshape( flex.grid(4,3) )
  from scitbx.linalg import svd
  svd_obj = svd.real( a.deep_copy(),True,True )
  a_pseudo_inverse = svd_obj.pseudo_inverse()
  print list( a_pseudo_inverse )
  print list(svd_obj.sigma) 

  rows,cols =  a_pseudo_inverse.focus()
  print rows, cols
  cal_back = []
  obs_sets = [ t0_c2 ,  t1_c2 ,  t2_c2,  t3_c2 ]
  deconvol = []
  for ith in range(rows):
    result = t0_c2*0.0
    for jset in range(cols):
      result += a_pseudo_inverse[ (ith,jset) ]*obs_sets[jset]
    deconvol.append( result )
  write_generic_image(deconvol,'deconvol.data')  


  MM = 10
  # DO T=0

  all_images = []

  res = c2_fball.deep_copy()*0.0
  synth = c2_fcombo*t0_s[0] + c2_fball*t0_s[1] + c2_frod*t0_s[2]
  write_generic_image(synth,"t0_calc.dat")

  for ii in range(MM):
    print "Do image ", ii+1
    ri =  eg.generate_random_image(t0_s,1.0)
    pc2 = int_obj.build_c2_function( ri )
    res += pc2
  res = res/MM
  t0_c2_obs = res.deep_copy()
  write_generic_image(res,"t0_obs.dat") 
    

  # Do T=1
  res = c2_fball.deep_copy()*0.0
  synth = c2_fcombo*t1_s[0] + c2_fball*t1_s[1] + c2_frod*t1_s[2]
  write_generic_image(synth,"t1_calc.dat")
  for ii in range(MM):
    print "Do image ", ii+1
    ri =  eg.generate_random_image(t1_s,1.0)
    res += int_obj.build_c2_function( ri )
  res = res/MM
  t1_c2_obs = res.deep_copy()
  write_generic_image(res,"t1_obs.dat") 

  # DO T=2
  res = c2_fball.deep_copy()*0.0
  synth = c2_fcombo*t2_s[0] + c2_fball*t2_s[1] + c2_frod*t2_s[2]
  write_generic_image(synth,"t2_calc.dat")
  for ii in range(MM):
    print "Do image ", ii+1
    ri =  eg.generate_random_image(t2_s,1.0)
    res += int_obj.build_c2_function( ri )
  res = res/MM
  t2_c2_obs = res.deep_copy()
  write_generic_image(res,"t2_obs.dat") 


  # DO T=3
  res = c2_fball.deep_copy()*0.0
  synth = c2_fcombo*t3_s[0] + c2_fball*t3_s[1] + c2_frod*t3_s[2]
  write_generic_image(synth,"t3_calc.dat")
  for ii in range(MM):
    print "Do image ", ii+1
    ri =  eg.generate_random_image(t3_s,1.0)
    res += int_obj.build_c2_function( ri )
  res = res/MM
  t3_c2_obs = res.deep_copy()
  write_generic_image(res,"t3_obs.dat")

  cal_back = []
  obs_sets = [ t0_c2_obs ,  t1_c2_obs ,  t2_c2_obs,  t3_c2_obs ]
  #deconvol = []
  for ith in range(rows):
    result = t0_c2*0.0
    for jset in range(cols):
      result += a_pseudo_inverse[ (ith,jset) ]*obs_sets[jset]
    deconvol.append( result )
  write_generic_image(deconvol,'deconvol_obs.data')



def tst_zernike_expansion():
  N=80
  length = 60
  width=10
  radius=15
  from sastbx.fXS import basic_shapes
  rod = enlarge(basic_shapes.rod(N,length,width,1.0),3,3,True)
  #ball = enlarge(basic_shapes.ball(N,radius,1.0),3,3,True)
   
  zernike_expansion(rod,5)
  write_generic_image(rod,'this_rod.dat')



if __name__ == "__main__":
  #run(sys.argv[1],sys.argv[2])
  tst_shapes()
  #tst_zernike_expansion()




