from scitbx import math
from stdlib import math as smath
from scitbx.array_family import flex
from sastbx import zernike_model as zm
from sastbx.intensity import sas_library
import os,sys
from iotbx import pdb
import time

def linear_fit(x,y,s):
  var = s*s
  sum_x2 = flex.sum( x*x/var )
  sum_xy = flex.sum( x*y/var )
  sum_x  = flex.sum( x / var )
  sum_y  = flex.sum( y / var )
  N = x.size()
  sum_inv_var = flex.sum(1.0/var)
  det = sum_inv_var * sum_x2 - sum_x * sum_x
  scale = (sum_inv_var * sum_xy - sum_x*sum_y ) / det
  offset = (sum_x2*sum_y - sum_x * sum_xy) /det
  return scale, offset


def get_chi_score( x, y, s):
  scale, offset = linear_fit( x, y, s )
  x = x*scale + offset
  score = flex.mean_sq( (y-x)/s )
  return score, x


class xyz2znk(object):
  def __init__(self, xyz, abs_I0,nmax=30, density=None, external_rmax=-1, np=50, splat_range=1, uniform=False, fix_dx=True, default_dx=0.5, fraction=0.7, protein=0.44, sol_layer=0.03, solvent=0.334, layer=4):
    self.nmax=nmax
    self.abs_I0 = abs_I0
    self.density=flex.double(xyz.size(), 1.0)
    self.xyz = xyz
    self.external_rmax=external_rmax
    self.default_dx = default_dx
    self.zernike_mom=math.nlm_array(nmax)
    self.protein=protein
    self.solvent=solvent*1.06 ## empirical parameter
    self.sol_layer=sol_layer/2.0 ## empirical parameter
    self.layer_thick=layer
    self.calc_mom(np, splat_range, uniform, fix_dx, default_dx, fraction)

  def calc_mom(self, np, splat_range, uniform, fix_dx, default_dx, fraction, pdb_out=False):
    if(self.xyz.size() == 0 ):
      print "no atom found, double check the model"
      return
    self.voxel_obj = math.sphere_voxel(np,splat_range+self.layer_thick,uniform,fix_dx,self.external_rmax, default_dx, fraction,self.xyz, self.density)
    self.volume_unit=self.default_dx**3.0
    self.total_volume=self.voxel_obj.occupied_sites()*self.volume_unit
    np=self.voxel_obj.np()
    self.rmax=self.voxel_obj.rmax()/fraction
    dx_used = self.rmax/np
    new_layer_thick=int( self.rmax*(1-fraction)/dx_used)-1
    if( new_layer_thick < self.layer_thick ):
      print "solvent layer extended too much, and thickness is reset to ", new_layer_thick
      self.layer_thick=new_layer_thick
    self.grid_obj = math.sphere_grid(np, self.nmax)

    ## displaced solvent ##
    border_list=self.voxel_obj.border(self.layer_thick-1)  ## only inner section left now
    inner_sites=self.voxel_obj.occupied_sites()
    self.grid_obj.clean_space(self.voxel_obj, pdb_out)
    self.grid_obj.construct_space_sum()
    self.mom_obj_s = math.zernike_moments(self.grid_obj, self.nmax)
    self.zernike_mom_s=self.mom_obj_s.moments()

    ## protein ##
    self.voxel_obj = math.sphere_voxel(np,splat_range,uniform,fix_dx,self.external_rmax, default_dx, fraction,self.xyz, self.density)
    real_inner_sites=self.voxel_obj.occupied_sites()
    self.grid_obj.clean_space(self.voxel_obj, pdb_out) ## output protein
    self.grid_obj.construct_space_sum()
    self.mom_obj_p = math.zernike_moments(self.grid_obj, self.nmax)
    self.zernike_mom_p=self.mom_obj_p.moments()

    ## border layer ##
    self.inner_volume=inner_sites*self.volume_unit
    self.grid_obj.construct_space_sum_via_list(border_list)
    self.mom_obj_b = math.zernike_moments(self.grid_obj, self.nmax)
    self.zernike_mom_b=self.mom_obj_b.moments()
    ## scale the densities ##
    self.density_scale=self.voxel_obj.weight_sum()
    print self.density_scale
    #self.density_scale=float(real_inner_sites)/inner_sites
    self.sol_layer=self.sol_layer/self.density_scale
    #self.solvent = self.solvent*self.density_scale
    ##self.protein=self.protein*(1.0+0.1*smath.exp(-inner_sites/134000.0) )
    ### adding moment vectors ##
    self.zernike_mom_coef  = self.zernike_mom_p.coefs()*self.protein #-self.solvent)
    self.zernike_mom_coef -= self.zernike_mom_s.coefs()*self.solvent
    self.zernike_mom_coef += self.zernike_mom_b.coefs()*self.sol_layer
    self.nlm = self.zernike_mom.nlm()
    self.zernike_mom.load_coefs( self.nlm, self.zernike_mom_coef )
    return

  def calc_intensity(self, q_array=None, rmax=None):
    return self.calc_old_intensity(q_array, rmax)

  def calc_new_intensity(self, q_array=None, rmax=None):
    if rmax is None:
      rmax=self.rmax
    if q_array is None:
      q_array=flex.double( range(51) )/100.0
    self.intensity = q_array*0.0
    self.intensity_vac = q_array*0.0
    self.intensity_sol = q_array*0.0
    self.intensity_lay = q_array*0.0
    scale_p = flex.exp(-0.23*q_array*q_array)
    scale_s = flex.exp(0.5*q_array*q_array)
    # Get the coefs #
    self.zp_coef=self.zernike_mom_p.coefs()*self.protein
    self.zs_coef=self.zernike_mom_s.coefs()*self.solvent
    self.zb_coef=self.zernike_mom_b.coefs()*self.sol_layer

    ii = 0 # starting with the first q
    for qq, sp, ss in zip( q_array, scale_p, scale_s):
      #Make a zernike model for particular q value qq #
      this_z_model=zm.zernike_model(self.zernike_mom, flex.double([qq]), rmax, self.nmax)
      self.mom_update(sp,ss)
      this_intensity = this_z_model.calc_intensity_nlm(self.zernike_mom)
      this_intensity_vac = this_z_model.calc_intensity_nlm(self.zernike_mom_p)
      this_intensity_lay = this_z_model.calc_intensity_nlm(self.zernike_mom_b)
      this_intensity_sol = this_z_model.calc_intensity_nlm(self.zernike_mom_s)
      self.intensity_vac[ii]=( this_intensity_vac[0] )
      self.intensity[ii]=( this_intensity[0] )
      self.intensity_sol[ii]=( this_intensity_sol[0] )
      self.intensity_lay[ii]=( this_intensity_lay[0] )
      ii = ii + 1
    ## Normalized to Absolute I(0)
    self.i_znk_0 = self.intensity_vac[0]
    self.intensity_vac = self.intensity_vac*(self.abs_I0/self.intensity_vac[0])
    self.intensity = self.intensity*(self.abs_I0/self.i_znk_0)
    self.intensity_lay = self.intensity_lay*(self.abs_I0/self.i_znk_0)
    self.intensity_sol = self.intensity_sol*(self.abs_I0/self.i_znk_0)
    return self.intensity, self.intensity_vac, self.intensity_sol, self.intensity_lay

  def mom_update(self, sp,ss):
    this_mom_p = self.zp_coef*sp
    this_mom_s = self.zs_coef*ss
    this_mom_b = self.zb_coef*ss
    this_mom = this_mom_p - this_mom_s + this_mom_b
    self.zernike_mom.load_coefs( self.nlm, this_mom)
    self.zernike_mom_p.load_coefs(self.nlm, this_mom_p)
    self.zernike_mom_s.load_coefs(self.nlm, this_mom_s)
    self.zernike_mom_b.load_coefs(self.nlm, this_mom_b)
    return

  def calc_old_intensity(self, q_array=None, rmax=None):
    if rmax is None:
      rmax=self.rmax
    if q_array is None:
      q_array=flex.double( range(51) )/100.0
    #scales = flex.exp(-0.23*q_array*q_array)
    #scales = scales * scales   # directly scaling intensity, not amplitudes
    scales = flex.double(q_array.size(), 1)
    self.z_model=zm.zernike_model(self.zernike_mom, q_array, rmax, self.nmax)
    self.intensity = self.z_model.calc_intensity_nlm(self.zernike_mom)
    self.intensity_vac = self.z_model.calc_intensity_nlm(self.zernike_mom_p)
    self.intensity_lay = self.z_model.calc_intensity_nlm(self.zernike_mom_b)*self.sol_layer**2.0
    self.intensity_sol = self.z_model.calc_intensity_nlm(self.zernike_mom_s)*self.solvent**2.0
    ## Normalized to Absolute I(0)
    self.i_znk_0 = self.intensity_vac[0]*self.protein*self.protein
    self.intensity_vac = self.intensity_vac*(self.abs_I0/self.intensity_vac[0]) * scales
    self.intensity = self.intensity*(self.abs_I0/self.i_znk_0) * scales
    self.intensity_lay = self.intensity_lay*(self.abs_I0/self.i_znk_0) * scales
    self.intensity_sol = self.intensity_sol*(self.abs_I0/self.i_znk_0) * scales
    return self.intensity, self.intensity_vac, self.intensity_sol, self.intensity_lay

  def calc_intensity_no_scale(self, q_array):
    rmax=self.rmax
    self.z_model=zm.zernike_model(self.zernike_mom, q_array, rmax, self.nmax)
    self.intensity = self.z_model.calc_intensity_nlm(self.zernike_mom)
    return self.intensity

  def optimize_solvent(self,data):
    self.optimize_old_solvent(data)
    return

  def optimize_new_solvent(self,data):
    scales=flex.double(range(5))*self.sol_layer
    self.best_score=1e8
    self.best_scale=0.0
    self.best_i_calc=None
    for s in scales:
      self.sol_layer=s
      self.calc_intensity(data.q)
      score, i_calc = get_chi_score( self.intensity, data.i, data.s )
      if (self.best_score > score ):
        self.best_score = score
        self.best_scale = s
        self.best_i_calc = i_calc.deep_copy()
    self.sol_layer=self.best_scale
    return


  def optimize_old_solvent(self, data):
    scales=flex.double(range(21))/5.0*self.sol_layer
    p_coef=self.zernike_mom_p.coefs()*self.protein - self.zernike_mom_s.coefs()*self.solvent
    b_coef=self.zernike_mom_b.coefs()
    self.z_model=zm.zernike_model(self.zernike_mom, data.q, self.rmax, self.nmax)
    self.best_score=1e8
    self.best_scale=0.0
    self.best_i_calc=None
    for s in scales:
      this_coef=p_coef+b_coef*s
      self.zernike_mom.load_coefs( self.nlm, this_coef )
      i_calc = self.z_model.calc_intensity_nlm( self.zernike_mom )
      score, i_calc = get_chi_score( i_calc, data.i, data.s )
      if (self.best_score > score ):
        self.best_score = score
        self.best_scale = s
        self.best_i_calc = i_calc.deep_copy()
    self.sol_layer=self.best_scale
    return



  def summary(self):
    output='\n'
    output=output+"  Protein Radius is  %11.3f\n"%self.rmax
    output=output+"  Grid resolution is %11.3f\n"%self.default_dx
    output=output+"  Electron densities:\n"
    output=output+"    Protein  ------- %11.3f\n"%self.protein
    output=output+"    Solvent  ------- %11.3f\n"%(self.solvent/self.density_scale)
    output=output+"    Layer Contrast - %11.3f\n"%(self.sol_layer*self.density_scale)
    output=output+"  Protein Volume is  %11.3f\n"%self.inner_volume
    output=output+"  Layer Volume is    %11.3f\n"%(self.total_volume-self.inner_volume)
    output=output+"  Total Volume is    %11.3f\n"%self.total_volume
    output=output+"\n"
    return output

def get_density(atoms):
  elements=atoms.extract_element()
  density = flex.double( elements.size(),0 )
  eles=[' C', ' N', ' O', ' S', ' P', ' H']
  nes =[6, 7, 8, 16, 15, 1]
  for ele, nelectron in zip( eles, nes ):
    sel=flex.bool(elements==ele)
    density.set_selected(sel, nelectron)
  #print list(elements[1:10])
  #print list(density[1:10])
  return density

def count_electron(atoms):
  elements=atoms.extract_element()
  eles=[' C', ' N', ' O', ' S', ' P', ' H']
  nes =[6, 7, 8, 16, 15, 1]
  tot_n_e = 0
  tot_atm = 0
  for ele, nelectron in zip( eles, nes ):
    sel=flex.bool(elements==ele)
    isel=sel.iselection()
    tot_n_e = tot_n_e + isel.size()*nelectron
    tot_atm = tot_atm + isel.size()
  return tot_n_e, tot_atm

def calc_abs_Io(atoms, explicit_H=False):
  n_electron, n_atm = count_electron(atoms)
  if(not explicit_H):
    n_electron = n_electron + n_atm
  abs_Io = n_electron*n_electron
  return abs_Io



def tst(pdb_file):
  pdbi = pdb.hierarchy.input(file_name=pdb_file)
  if(len( pdbi.hierarchy.models() ) == 0):
    return None,None,None
  atoms = pdbi.hierarchy.models()[0].atoms()
  xyz = flex.vec3_double()
  for atom in atoms:
    if(not atom.hetero):
      xyz.append( atom.xyz )
  density=flex.double(xyz.size(),1.0)
  nmax=20
  abs_Io=calc_abs_Io(atoms)
  znk_model=xyz2znk(xyz,abs_Io,nmax,density=density)
  q_array = flex.double( range(51) )/100.0
  calc_i, calc_i_vac, i_sol, i_lay=znk_model.calc_intensity(q_array)
  for qq,ii,iv in zip(q_array, calc_i, calc_i_vac):
    print qq,ii,iv


if __name__ == "__main__":
  pdb_file=sys.argv[1]
  t1=time.time()
  tst(pdb_file)
  print time.time() - t1
