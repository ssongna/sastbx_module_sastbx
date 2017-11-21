from scitbx.array_family import flex
from mmtbx.monomer_library import server, pdb_interpretation
from iotbx import pdb
from sastbx.intensity import sas_library
from sastbx import intensity
from stdlib import math as smath
import libtbx.env_config
from libtbx.test_utils import run_command
import os, sys, time


class she(object):
  def __init__(self,pdb_file,q_array=None, rho=0.334, drho=0.03, max_L=15, max_i=13, implicit_hydrogens=True):
    if(q_array is not None):
      q_values = q_array
    else:
      q_values = flex.double( range(51) )/100.0
    # read in pdb file
    pdbi = pdb.hierarchy.input(file_name=pdb_file)
    atoms = pdbi.hierarchy.atoms()
    # predefine some arrays we will need
    atom_types = flex.std_string()
    dummy_atom_types = flex.std_string()
    radius= flex.double()
    b_values = flex.double()
    occs = flex.double()
    xyz = flex.vec3_double()
    # keep track of the atom types we have encountered
    dummy_at_collection = []
    for atom in atoms:
      b_values.append( atom.b )
      occs.append( atom.occ )
      xyz.append( atom.xyz )

    self.natom=xyz.size()

    dummy_ats= sas_library.read_dummy_type(file_name=pdb_file)
    for at in dummy_ats:
      if at not in dummy_at_collection:
        dummy_at_collection.append( at )

    # Hydrogen controls whether H is treated explicitly or implicitly
    #Hydrogen = True
    Hydrogen = not implicit_hydrogens

    radius_dict={}
    ener_lib=server.ener_lib()
    for dummy in dummy_at_collection:
      if(Hydrogen):
        radius_dict[dummy]=ener_lib.lib_atom[dummy].vdw_radius
      else:
        radius_dict[dummy]=ener_lib.lib_atom[dummy].vdwh_radius
      if(radius_dict[dummy] is None):
        radius_dict[dummy]=ener_lib.lib_atom[dummy].vdw_radius
    if(radius_dict[dummy] is None):
      print  "********************* WARNING WARNING  *************************"
      print  "Did not find atom type: ", dummy, "default value 1.58 A was used"
      print  "****************************************************************"
      radius_dict[dummy]=1.58

    for at in dummy_ats:
      dummy_atom_types.append( at)
      radius.append(radius_dict[at])

    Scaling_factors=sas_library.load_scaling_factor()

    #------------------
    #
    B_factor_on=False
    f_step= 0.8
    q_step= 0.01
    solvent_radius_scale=0.91
    protein_radius_scale=1.2
    delta=3.0
    #------------------

    scat_lib_dummy =  sas_library.build_scattering_library( dummy_at_collection,
                                                q_values,
                                                radius_dict,
                                                solvent_radius_scale,
                                                Hydrogen,
                                                Scaling_factors)

    model=intensity.model(xyz,
                          radius*protein_radius_scale,
                          b_values,
                          occs,
                          dummy_ats,
                          scat_lib_dummy,
                          B_factor_on)

    max_z_eps=0.02
    max_z=model.get_max_radius()*(q_values[-1]+max_z_eps)
    self.engine = intensity.she_engine( model, scat_lib_dummy,max_i,max_L,f_step, q_step,max_z, delta,rho,drho )
    self.engine.update_solvent_params(rho,drho)

  def update_coord(self, coords, indices):
    self.engine.update_coord(coords, indices)

  def get_intensity( self ):
    return self.engine.I()

  def get_all_blq( self ):
    return self.engine.get_all_coefs()

def linear_fit(x,y,s): # this is to fit y = s* x + o, not y = (x+o)*s
  var = s*s
  sum_x2 = flex.sum( x*x/var )
  #sum_y2 = flex.sum( y*y/var )
  sum_xy = flex.sum( x*y/var )
  sum_x  = flex.sum( x / var )
  sum_y  = flex.sum( y / var )
  N = x.size()
  sum_inv_var = flex.sum(1.0/var)
  det = sum_inv_var * sum_x2 - sum_x * sum_x
  scale = (sum_inv_var * sum_xy - sum_x*sum_y ) / det
  offset = (sum_x2*sum_y - sum_x * sum_xy) /det

  return scale, offset


class solvent_parameter_optimisation(object):
  def __init__(self, she_object, observed_data):
    # we'll only optimize the scale factor and form factor of excluded solvent
    self.rm = 1.62
    self.rho=0.334
    self.drho=0.03
    self.obs = observed_data
    self.she_object = she_object
    self.rm_fluct_scale = -(4.0*smath.pi/3.0)**1.5*smath.pi*flex.pow(self.obs.q,2.0)*self.rm**2.0
    ### setup the scan range ###
    self.default_a = 1.0
    self.a_range=flex.double(range(-10,11))/50.0+self.default_a
    self.drho_range = (flex.double(range(-10,21))/10.0+1.0)*self.drho

    self.scan()


  def scan(self):
    self.best_score, self.best_i_calc, self.best_scale, self.best_offset=self.target(self.drho, self.default_a)
    self.best_drho = self.drho
    self.best_a = self.default_a
    for drho in self.drho_range:
      for a in self.a_range:
        score, i_calc, s, off = self.target(drho, a)
        if(score < self.best_score ):
          self.best_score = score
          self.best_drho = drho
          self.best_a = a
          self.best_i_calc = i_calc.deep_copy()
          self.best_scale = s
          self.best_offset = off
    return

  def compute_scale_array(self, a):
    scales = a**3.0 * flex.exp(self.rm_fluct_scale*(a**2.0 - 1))
    return scales


  def get_scaled_data(self):
    return self.best_i_calc

  def get_scales(self):
    return self.best_scale, self.best_offset, self.best_drho, self.best_a*self.rm

  def target(self, drho, a):
    self.she_object.update_solvent_params(self.rho,drho)
    this_scale = self.compute_scale_array( a )
    i_calc = self.she_object.Iscale(this_scale)

    s, off = linear_fit( i_calc, self.obs.i, self.obs.s )
    i_calc = s*i_calc + off
    result = flex.sum(flex.pow2( (i_calc-self.obs.i)/self.obs.s ) )
    return result, i_calc, s, off
