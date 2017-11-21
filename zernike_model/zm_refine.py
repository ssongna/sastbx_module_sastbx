import sys, os, random
from stdlib import math as smath
from scitbx.array_family import flex
import iotbx.phil
from iotbx import pdb
from sastbx.data_reduction import saxs_read_write
from sastbx import zernike_model as zm
from sastbx.zernike_model.search_pdb import reduce_raw_data, get_mean_sigma
from sastbx.zernike_model import model_interface
from scitbx.math import zernike_align_fft as fft_align
import time

from scitbx import math
from scitbx.golden_section_search import gss
from iotbx.xplor import map as xplor_map
from cctbx import uctbx

from libtbx import easy_pickle

from sastbx.interface import get_input


master_params = iotbx.phil.parse("""\
zrefine{
  target = None
  .type=path
  .help="the experimental intensity profile"

  start = None
  .multiple=True
  .type=path
  .help="starting model in xplor format or pickle file of coefs"

  pdb = None
  .type=path
  .help="PDB model to be compared"

  rmax = None
  .type=float
  .help="estimated rmax of the molecule"

  qmax = 0.15
  .type=float
  .help="maximum q value where data beyond are disgarded"

  n_trial=5
  .type=int
  .help="number of refinement trials"

  nmax=10
  .type=int
  .help="maximum order of zernike expansion"

  splat_range=1
  .type=int
  .help="depth of alterable region from original surface"

  prefix='prefix'
  .type=path
  .help="prefix of the output files"
}
""")

banner = "-------------------Shape Refinement-------------------"

def xplor_map_type(m,N,radius,file_name='map.xplor'):
  gridding = xplor_map.gridding( [N*2+1]*3, [0]*3, [2*N]*3)
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  m.reshape( grid )
  uc = uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90")
  xplor_map.writer( file_name, ['no title lines'],uc, gridding,m) # is_p1_cell=True)  # True)

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\n\t Usage: libtbx.python zm_refine.py <target=target.iq start=start.xplor rmax=rmax> [splat=splat prefix=prefix qmax=qmax]\n "


def run(args):
  params = get_input(args, master_params, "zrefine", banner, help)
  if params is None:
    return
  np=30
  nmax=params.zrefine.nmax
  start_file=params.zrefine.start
  target_file=params.zrefine.target
  rmax = params.zrefine.rmax
  qmax = params.zrefine.qmax
  prefix = params.zrefine.prefix
  splat_range = params.zrefine.splat_range
  pdb = params.zrefine.pdb
  n_trial = params.zrefine.n_trial
  if pdb is not None:
    pdb_obj = model_interface.container( pdbfile=pdb, rmax=rmax, nmax=nmax )
    pdb_nlm = pdb_obj.nlm_array
  else:
    pdb_nlm = None

  data = saxs_read_write.read_standard_ascii_qis( target_file )
  refine_obj = zm_xplor_refine( data, start_file, rmax, qmax=qmax, nmax=nmax, np_on_grid=np, prefix=prefix, splat_range=splat_range, pdb_nlm=pdb_nlm, n_trial=n_trial )

  if pdb is not None:
    out_pdb_filename = pdb.split('.')[0]+'_shift.pdb'
    pdb_obj.write_pdb( rmax=refine_obj.rmax, filename=out_pdb_filename )


class zm_xplor_refine(object):
  def __init__(self, data, xplor_file, rmax, qmax=0.15, nmax=20, np_on_grid=30, prefix='prefix', splat_range=1, n_trial=5, fraction=0.9, pdb_nlm=None, firstn=1):
    self.raw_data=data
    self.rmax = rmax/fraction
    self.fraction = fraction
    self.qmax = qmax
    self.nmax = nmax
    self.pdb_nlm = pdb_nlm
    self.load_maps(xplor_file)
    self.nmodel = len(xplor_file)
    self.zga = math.zernike_grid( self.np_on_grid, self.nmax, False )
    self.ngp = self.raw_map.size()  # number of grid point in 3D box
    self.n = self.np_on_grid*2+1
    self.n2 = self.n**2
    self.prefix=prefix+'_'
    self.splat_range = splat_range
    self.start_model = None
    neighbors = [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]
    self.neighbors = flex.int()
    for n in neighbors:
      self.neighbors.append( self.convert_indx_3_1( n ) )
    self.build_sphere_list()
    #### Labels for different regions ####
    self.solvent_label=0
    self.molecule_label=1
    self.surface_label=2
    self.nlm_array=math.nlm_array(nmax)
    self.nlm=self.nlm_array.nlm()
    self.counter = 0
    self.n_accept = 0
    self.bandwidth = min( smath.pi/rmax/2.0, 0.01 )
    self.scale_2_expt = flex.mean( self.raw_data.i[0:firstn] )
    self.data = reduce_raw_data(self.raw_data,self.qmax, self.bandwidth,level=0.00001 )
    self.data.i = self.data.i/self.scale_2_expt
    self.data.s = self.data.i

    self.build_starting_model()
    self.mark_mod_core_region(splat_range)
    self.n_trial = n_trial
    self.finals = []
    for ii in range( n_trial ):
      self.refine(ii)
      self.finals.append( self.best_nlm_coefs.deep_copy() )
    self.pair_align()

  def pair_align(self):
    ms = flex.double()
    ss = flex.double()
    tmp_nlm_array = math.nlm_array( self.nmax )
    for coef in self.finals:
      mean = abs( coef[0] )
      var = flex.sum( flex.norm( coef ) )
      sigma = smath.sqrt( var - mean*mean )
      ms.append( mean )
      ss.append( sigma)

    grids = flex.grid(self.n_trial, self.n_trial)
    self.cc_array=flex.double( grids, 1.0 )
    for ii in range( self.n_trial ):
      self.nlm_array.load_coefs( self.nlm, self.finals[ii] )
      for jj in range( ii ):
        tmp_nlm_array.load_coefs( self.nlm, self.finals[jj] )
        cc = fft_align.align( self.nlm_array, tmp_nlm_array, nmax=self.nmax, refine=True ).best_score
        cc = (cc-ms[ii]*ms[jj])/(ss[ii]*ss[jj])
        self.cc_array[(ii,jj)]=cc
        self.cc_array[(jj,ii)]=cc

    outfile = self.prefix+"pair.cc"
    comment = "# electron density correlation coefficient, < rho_1(r)*rho_2(r) >"
    out=open(outfile, 'w')
    print>>out, comment
    for ii in range(1,self.n_trial+1):
      print>>out,"%6d"%ii,
    print>>out, "   average"

    for ii in range(self.n_trial):
      for jj in range(self.n_trial):
        print>>out,"%6.3f"%self.cc_array[(ii,jj)],
      print>>out, flex.mean( self.cc_array[ii*self.n_trial:(ii+1)*self.n_trial] )
    out.close()




  def load_maps(self, files): # take about 0.5s to load one map, not big deal
    xplor_file = files[0]
    this_xplor = xplor_map.reader(xplor_file)
    self.raw_map = flex.double( this_xplor.data.size(), 0)
    threshold = flex.max( this_xplor.data )/3.0
    for ii in range( this_xplor.data.size() ):
      if( this_xplor.data[ii] > threshold ):
        self.raw_map[ii] = 1

    self.np_on_grid = (this_xplor.gridding.n[0]-1 ) /2

    for xplor_file in files[1:]:
      this_xplor = xplor_map.reader(xplor_file)
      threshold = flex.max( this_xplor.data )/3.0
      for ii in range( this_xplor.data.size() ):
        if( this_xplor.data[ii] > threshold ):
          self.raw_map[ii] += 1

  def refine(self, trial):
    print "--------------Trial %d-----------------"%trial, time.ctime()
    self.working_model = self.start_model.deep_copy()
    self.nlm_coefs = self.start_nlm_coefs.deep_copy()
    self.best_nlm_coefs = self.start_nlm_coefs.deep_copy()
    self.best_i = self.start_i.deep_copy()
    self.lowest_score = self.start_score
    init_scores = flex.double()
    for ii in range(10):
      self.modify()
      init_scores.append( self.target() )
    mean = flex.mean( init_scores )
    self.deltaS = smath.sqrt( flex.sum(flex.pow2(init_scores-mean) )/10.0 )
    self.T = self.deltaS * 100
    self.nsteps = 200
    self.score = mean
    self.working_model = self.start_model.deep_copy()
    while( self.T > self.deltaS/2.0):
      self.n_reject = 0
      for ii in range( self.nsteps ):
        self.move()
      print "Number of Accept/Reject: %d/%d"%(self.nsteps-self.n_reject, self.n_reject)
      if( self.n_reject > self.nsteps*0.9 ):
        print "Too Many rejections (%d), quit at temperature (%f)"%(self.n_reject, self.T)
        break
      self.T = self.T*0.9

    out = open(self.prefix+str(trial)+'_final.iq', 'w')
    self.nlm_array.load_coefs( self.nlm, self.best_nlm_coefs )
    best_i = self.zm.calc_intensity_nlm( self.nlm_array )
    best_i = best_i/best_i[0]*self.scale_2_expt
    for qq,ic,io in zip( self.data.q, best_i, self.data.i*self.scale_2_expt):
      print>>out, qq, ic, io
    out.close()
    print "total number of moves %d"%self.counter
    print "total number of accepted moves %d"%self.n_accept

    if (self.pdb_nlm is not None):
      align_obj = fft_align.align(self.pdb_nlm, self.nlm_array, nmax=self.nmax, refine=True)
      mean = abs( self.best_nlm_coefs[0] )
      var = flex.sum( flex.norm( self.best_nlm_coefs ) )
      sigma = smath.sqrt( var - mean*mean )
      cc = align_obj.best_score
      cc = ( cc - mean*self.pdb_m ) / ( sigma*self.pdb_s )
      print "C.C. (PDB, trial%6d) = %8.5f, Score = %8.5f"%(trial, cc, self.lowest_score)
      self.best_nlm_coefs = align_obj.moving_nlm.coefs()

    reconst_model = self.reconst_model( self.best_nlm_coefs )
    xplor_map_type( reconst_model, self.np_on_grid, self.rmax, file_name=self.prefix+str(trial)+'_final_rbt.xplor')
    xplor_map_type( self.best_model, self.np_on_grid, self.rmax, file_name=self.prefix+str(trial)+'_final.xplor')
    print "-----------End of Trial %d--------------"%trial, time.ctime()


  def reconst_model( self,coefs ):
    self.zga.load_coefs( self.nlm, coefs )
    map = flex.abs( self.zga.f() )
    return map

  def optimize_rmax( self ):
    rmax_min = self.rmax/2.0
    rmax_max = self.rmax*2.0
    rmax = gss( self.rmax_target, rmax_min, rmax_max, eps=0.5 )
    return rmax

  def rmax_target( self, rmax ):
    self.zm = zm.zernike_model( self.moments, self.data.q, rmax, self.nmax )
    calc_i = self.zm.calc_intensity( self.nn )
    calc_i = calc_i / calc_i[0]
    chi_score= ( (calc_i-self.data.i)/self.data.s ).norm()
    return chi_score


  def build_starting_model(self):
    grids = flex.grid([self.np_on_grid*2+1]*3)
    self.start_model = flex.double(grids,0.0)
    ####   BUILD STARTING MODEL  #####
    max_map = flex.max( self.raw_map )
    self.molecule = flex.int()
    self.mod_list_1d = flex.int()
    self.mod_list = flex.vec3_double()
    cutoff = self.nmodel/2.0

    distance_indx = flex.vec3_double()
    for ii in range(self.ngp ):
      if(self.raw_map[ii] >= cutoff):
        self.start_model[ii]=1.0
        self.molecule.append( ii )  ## Need to be moved to mark region


    self.working_model= self.start_model.deep_copy()  # Make a working model
    self.best_model = self.start_model.deep_copy()    # Make a starting model

    #### Reusable Objects for moments calculation ####
    self.grid_obj=math.sphere_grid(self.np_on_grid, self.nmax)
    self.grid_obj.construct_space_sum_via_list( self.molecule, self.start_model.as_1d() )
    self.moment_obj = math.zernike_moments( self.grid_obj, self.nmax )
    self.moments = self.moment_obj.moments()
    self.nn = self.moment_obj.fnn()

    print "Starting Rmax is %5.2f"%(self.rmax)
    self.rmax = self.optimize_rmax()
    print "Optimized Rmax is %5.2f"%(self.rmax)

    self.zm = zm.zernike_model( self.moments, self.data.q, self.rmax, self.nmax)
    self.calc_i = self.zm.calc_intensity(self.nn)
    self.calc_i = self.calc_i / self.calc_i[0]
    self.start_i = self.calc_i.deep_copy()
    self.start_nlm_coefs = self.moments.coefs().deep_copy()
    self.best_nlm_coefs = self.start_nlm_coefs.deep_copy()
    self.start_score= ( (self.calc_i-self.data.i)/self.data.s ).norm()
    out = open(self.prefix+'start.iq', 'w')
    for qq,ic,io in zip( self.data.q, self.calc_i*self.scale_2_expt, self.data.i*self.scale_2_expt):
      print>>out, qq, ic, io
    out.close()

    self.nlm_array.load_coefs( self.nlm, self.start_nlm_coefs )
    self.start_m, self.start_s = get_mean_sigma( self.nlm_array )

    if (self.pdb_nlm is not None):
      self.pdb_m, self.pdb_s = get_mean_sigma( self.pdb_nlm )
      align_obj = fft_align.align(self.pdb_nlm, self.nlm_array, nmax=self.nmax, refine=True)
      cc = align_obj.best_score
      cc = ( cc - self.start_m*self.pdb_m ) / ( self.start_s*self.pdb_s )
      print "C.C. (PDB, Start) = %8.5f, Score = %8.5f"%(cc, self.start_score)

    xplor_map_type( self.raw_map, self.np_on_grid, self.rmax, file_name=self.prefix+'start.xplor')
    print "Fraction: ", flex.sum(self.start_model)/(self.np_on_grid**3.0)/8.0

  def build_sphere_list(self):
    self.indx_list=flex.int()
    indx_range = range(self.np_on_grid*2+1)
    np = self.np_on_grid
    np2 = np**2
    for ix in indx_range:
      for iy in indx_range:
        for iz in indx_range:
          if( (ix-np)**2 + (iy-np)**2 + (iz-np)**2 < np2 ):
            self.indx_list.append( self.convert_indx_3_1( (ix,iy,iz) ) )
    return

  def build_boundary(self, label):
    if( label == 0):
      mark_label = 1
    else:
      mark_label = 0
    new_surf_list = flex.int()
    if (self.surf_list.size() == 0 ):
      go_through_list = self.indx_list
    else:
      go_through_list = self.surf_list

    for pt in go_through_list:
          if( self.mask[pt] == label ):
            if( self.mark_surface( pt, mark_label) ):
              new_surf_list.append( (pt) )

    self.surf_list=new_surf_list.deep_copy()
    for indx in new_surf_list:
      self.mask[ indx ] = mark_label
    return


  def mark_surface( self, indx, mark_label):
    for n in self.neighbors:
      t= indx + n
      if (self.mask[t] == mark_label):
        return True
    return False

  def mark_mod_core_region(self, depth):
    self.mask = self.start_model.deep_copy()
    self.surf_list=flex.int()  # clean the surface list
    for dd in range(1, depth+1 ):
      self.erode(1)

    self.mod_list = flex.vec3_double()
    self.mod_list_1d = flex.int()
    self.mask = self.start_model-self.mask
    for indx in self.indx_list:
      if( self.mask[indx] > 0 ):
        self.mod_list_1d.append( indx )
        self.mod_list.append( self.convert_indx_1_3( indx ) )

    self.mask = self.start_model.deep_copy()
    self.surf_list=flex.int()
    for dd in range(1, depth+1 ):
      self.dilate(0)

    self.mask = self.mask - self.start_model
    for indx in self.indx_list:
      if( self.mask[indx] > 0 ):
        self.mod_list_1d.append( indx )
        self.mod_list.append( self.convert_indx_1_3( indx ) )

    print "Size of Changable Region", self.mod_list.size()

    #xplor_map_type( self.working_model, self.np_on_grid, self.rmax, file_name=self.prefix+'_mod.xplor')

  def erode(self, d):
    self.build_boundary(d)
    self.working_model = self.working_model + self.mask

  def dilate(self,d):
    self.build_boundary(d)
    self.working_model = self.working_model + self.mask


  def modify(self):
    self.space_alter_list=flex.int()

    indx = int( random.random()*self.mod_list.size() )
    center_1d  = self.mod_list_1d[ indx ]
    center_indx  = self.mod_list[ indx ]
    radius = self.splat_range + 2

    if( self.working_model[ center_1d ] == 0 ):
      self.label = 0
      self.add = False   ### shrink
    else:
      self.label = 1
      self.add = True   ### grow

    dist = (self.mod_list - center_indx ).norms()
    for d, m1d in zip( dist, self.mod_list_1d ):
      if(d< radius and self.working_model[m1d] != self.label):
        self.space_alter_list.append( m1d )

    space_sum = self.grid_obj.construct_space_sum_via_list( self.space_alter_list )
    self.moment_obj.calc_moments( space_sum.as_1d() )
    self.delta_nlm_coefs = self.moment_obj.moments().coefs()
    if( self.add ):
      self.nlm_coefs = self.nlm_coefs + self.delta_nlm_coefs
    else:
      self.nlm_coefs = self.nlm_coefs - self.delta_nlm_coefs
    self.nlm_array.load_coefs( self.nlm_array.nlm(), self.nlm_coefs )
    self.counter += 1
    return


  def convert_indx_3_1(self, this_indx):
    return this_indx[0]*self.n2+this_indx[1]*self.n+this_indx[2]

  def convert_indx_1_3(self, this_indx):
    a = this_indx / self.n2
    tmp = this_indx - a* self.n2
    b = tmp / self.n
    c = tmp - b * self.n
    return ( a, b, c )

  def commit_changes(self):
    for m in self.space_alter_list:
      self.working_model[ m ] = self.label
    return

  def revoke_changes(self):
    if( self.add ):
      self.nlm_coefs = self.nlm_coefs - self.delta_nlm_coefs
    else:
      self.nlm_coefs = self.nlm_coefs + self.delta_nlm_coefs

  def target(self):
    self.calc_i = self.zm.calc_intensity_nlm( self.nlm_array )
    self.calc_i = self.calc_i/self.calc_i[0]
    score = ((self.calc_i-self.data.i)/self.data.s).norm()
    return score

  def move(self):
    self.modify()
    score = self.target() # metropolis criterion
    if( score < self.score ):
      self.score = score
      self.n_accept += 1
      self.commit_changes()
      if( score < self.lowest_score):
        self.lowest_score = self.score
        self.best_model = self.working_model.deep_copy()
        self.best_i = self.calc_i.deep_copy()
        self.best_nlm_coefs = self.nlm_coefs.deep_copy()
    else:
      delta = (self.score - score ) / self.T
      delta = smath.exp( delta )
      if( random.random() < delta ):
        #print "Delta, T", delta, self.T, score, self.score
        self.revoke_changes()
        self.n_reject += 1
      else:
        self.score = score
        self.n_accept += 1
        self.commit_changes()
    return


if __name__ == "__main__":
   t1 = time.time()
   args=sys.argv[1:]
   run(args)
   t2 = time.time()
   print "Time spent: ", t2-t1
