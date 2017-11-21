import math

from cctbx.eltbx.xray_scattering import gaussian, wk1995
from cctbx.xray import scattering_type_registry
from scitbx.array_family import flex

from sastbx.fXS import solvent_accessible_area

class solvent_model(object):
  """
  Nucleic Acids Research, 2010, Vol. 38

  FoXS web server

  f_i(q) = f_v(q) - c_1 f_s(q) + c_2 s_i f_w(q)          (eq 2)

  """
  def __init__(self):
    self.interpreted_pdb = None
    self.xyz = None
    self.scattering_types = None
    self.element_types = None
    self.explicit_h = False
    self.bulk_solvent_scale = 0.334
    self.boundary_layer_scale = 0.1
    self.probe_radius = 1.4
    self.default_radius = 1.7
    self.fudge_factor = 0.75

  def add_bulk_solvent(self,reference_registry):
    assert (self.interpreted_pdb is not None)
    assert (self.xyz is not None)

    # determine scattering types
    self.scattering_types = self.interpreted_pdb.all_chain_proxies.\
                            nonbonded_energy_type_registry.symbols
    self.element_types = self.interpreted_pdb.all_chain_proxies.\
                         scattering_type_registry.symbols
    assert(len(self.scattering_types) == len(self.element_types))
    for i in xrange(len(self.scattering_types)):
      if (len(self.scattering_types[i]) == 0):
        self.scattering_types[i] = self.element_types[i]

    atom_library = self.interpreted_pdb.ener_lib.lib_atom
    new_registry = scattering_type_registry()
    for s_type,e_type in zip(self.scattering_types,self.element_types):
      if (new_registry.has_key(s_type) is False):
        new_registry.process(s_type)

        # old parameters
        g = reference_registry.gaussian(e_type)
        a = list(g.array_of_a())
        b = list(g.array_of_b())
        c = g.c()
        use_c = g.use_c()

        # add new parameters for bulk solvent
        if (s_type in atom_library.keys()):
          if (self.explicit_h):
            vdw_radius = atom_library[s_type].vdw_radius
          else:
            vdw_radius = atom_library[s_type].vdwh_radius
        else:
          vdw_radius = self.default_radius
        if (vdw_radius is None):
          vdw_radius = self.default_radius
        v = 4.0/3.0*math.pi*math.pow(vdw_radius,3) * self.fudge_factor
        a.append(-self.bulk_solvent_scale*v)
        b.append(math.pi*math.pow(v,2.0/3.0))

        # finalize form factor
        new_g = gaussian(a,b,c,use_c)
        new_registry.assign(s_type,new_g)
      else:
        new_registry.process(s_type)

    ## # add O for boundary layer solvent
    ## o_g = reference_registry.gaussian('O')
    ## new_registry.process('HOH')
    ## new_registry.assign('HOH',o_g)

    return new_registry

  def add_boundary_layer_solvent(self,reference_registry):

    assert(len(self.scattering_types) == len(self.xyz))

    # get van der Waals radii
    atom_library = self.interpreted_pdb.ener_lib.lib_atom
    vdw_radii = flex.double(len(self.xyz))
    for i in xrange(len(self.xyz)):
      s_type = self.scattering_types[i]
      if (s_type in atom_library.keys()):
        if (self.explicit_h):
          vdw_radius = atom_library[s_type].vdw_radius
        else:
          vdw_radius = atom_library[s_type].vdwh_radius
        if (vdw_radius is None):
          vdw_radii[i] = self.default_radius
        else:
          vdw_radii[i] = vdw_radius
      else:
        vdw_radii[i] = self.default_radius

    # get solvent accessible areas
    n_points = 1000
    indices = flex.int(len(self.xyz))
    for i in xrange(len(self.xyz)):
      indices[i] = i
    areas = solvent_accessible_area(self.xyz,vdw_radii,indices,
                                    self.probe_radius,n_points)

    # construct array of scaling factors
    boundary_layer_scaling_factors = flex.double(len(self.xyz),0.0)
    for i in xrange(len(self.xyz)):
      boundary_layer_scaling_factors[i] = self.boundary_layer_scale*\
                                          areas[i]/n_points

    return boundary_layer_scaling_factors

    ## # form factor for water (O)
    ## g_w = wk1995('O',True).fetch()
    ## a_w = g_w.array_of_a()
    ## b_w = g_w.array_of_b()
    ## c_w = g_w.c()
    ## use_c_w = g_w.use_c()

    ## new_registry = scattering_type_registry()
    ## new_scattering_types = flex.std_string(len(self.scattering_types))
    ## for i in xrange(len(self.xyz)):
    ##   reference_s_type = self.scattering_types[i]
    ##   s_type = self.scattering_types[i] + '_' + str(areas[i])
    ##   new_scattering_types[i] = s_type
    ##   if (new_registry.has_key(s_type) is False):
    ##     new_registry.process(s_type)

    ##     # old parameters
    ##     g = reference_registry.gaussian(reference_s_type)
    ##     a = list(g.array_of_a())
    ##     b = list(g.array_of_b())
    ##     c = g.c()
    ##     use_c = g.use_c()

    ##     # add new parameters for boundary layer
    ##     if (areas[i] != 0):
    ##       scale = self.boundary_layer_scale*areas[i]/n_points
    ##       for i in xrange(len(a_w)):
    ##         a.append(scale)
    ##         b.append(b_w[i])
    ##       if (use_c_w):
    ##         if (use_c):
    ##           c += scale * c_w
    ##         else:
    ##           use_c = use_c_w
    ##           c = scale * c_w

    ##     # finalize form factor
    ##     new_g = gaussian(a,b,c,use_c)
    ##     new_registry.assign(s_type,new_g)
    ##   else:
    ##     new_registry.process(s_type)

    ## self.scattering_types = new_scattering_types
    ## return new_registry
