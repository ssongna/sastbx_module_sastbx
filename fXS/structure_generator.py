import math,random

from cctbx.xray import scattering_type_registry
from libtbx.utils import Sorry
from mmtbx.monomer_library import server, pdb_interpretation
from scitbx.array_family import flex
from scitbx.math import basic_statistics

from sastbx.fXS.solvent_model import solvent_model

# =============================================================================
class species(object):
  """
  =============================================================================
  Molecular Species Class
  Data structure for structure_generator

  Useful accessible data members:
    n_copies - the number of copies of the structure
    pdb_input - the iotbx_pdb_ext.input object for the structure
    n_electrons - the sum of the scattering factors at diffraction angle 0
    radius - the minimum bounding radius of the structure

  -----------------------------------------------------------------------------
  """
  def __init__(self):
    self.n_copies = None
    self.pdb_input = None
    self.scattering_type_registry = None
    self.radius = None
    self.n_electrons = None
    self.xyz = None
    self.scattering_types = None
    self.boundary_layer_scaling_factors = None

# =============================================================================
class structure_generator(object):
  """
  =============================================================================
  Model Generator Class
  Generates copies of models with random rotations and translations

  Arguments:
    None

  Useful accessible methods:
    __len__() - returns number of species
    add_species(iotbx_pdb_ext.input,int) - adds n_copies of a model
    randomize() - generates random rotations and translations for added species

  Useful accessible data members:
    max_attempts - maximum number of attempts at placing molecule
    min_separation - minimum distance between molecules
    box_size - the size of the experimental box (cube)

  Notes:
    All distances (min_separation, box_size, radius) are in Angstroms
    All species should be added before anything else

  -----------------------------------------------------------------------------
  """
  def __init__(self):

    # data storage
    self.max_attempts = 100
    self.min_separation = 100.0
    self.box_size = 10000.0
    self.use_solvent = True
    self.species = list()
    self.rotations = list()
    self.translations = list()
    self.random = random.Random()
    self.total_electrons = None
    self.neighbors = dict()
    self.gridsize = 2.0 * self.min_separation

  def __len__(self):
    return len(self.species)

  def add_species(self,pdb_input=None,n_copies=1,form_factor_table='WK1995'):
    if (n_copies <= 0):
      raise Sorry('n_copies has to be greater than or equal to 1')

    # construct entry
    s = species()
    s.n_copies = n_copies
    s.pdb_input = pdb_input

    # find center and radius of sphere enclosing molecule
    atoms = s.pdb_input.atoms()
    x = flex.double(len(atoms))
    y = flex.double(len(atoms))
    z = flex.double(len(atoms))
    for i in xrange(len(atoms)):
      x[i] = atoms[i].xyz[0]
      y[i] = atoms[i].xyz[1]
      z[i] = atoms[i].xyz[2]
    x_stats = basic_statistics(x)
    y_stats = basic_statistics(y)
    z_stats = basic_statistics(z)
    r = flex.double( ( (x_stats.max - x_stats.min)/2.0,
                       (y_stats.max - y_stats.min)/2.0,
                       (z_stats.max - z_stats.min)/2.0 ) )
    center = (x_stats.mean,y_stats.mean,z_stats.mean)
    s.radius = r.norm()

    # center model at origin
    center = flex.double(center)
    s.xyz = flex.vec3_double(len(atoms))
    for i in xrange(len(atoms)):
      s.xyz[i] = tuple(flex.double(atoms[i].xyz) - center)

    # determine scattering types
    mon_lib_srv = server.server()
    ener_lib = server.ener_lib()
    interpreted_pdb = pdb_interpretation.process\
                        (mon_lib_srv=mon_lib_srv,ener_lib=ener_lib,
                         pdb_inp=s.pdb_input)
    s.scattering_types = interpreted_pdb.all_chain_proxies.\
                         scattering_type_registry.symbols
    s.scattering_type_registry = scattering_type_registry()
    for f in s.scattering_types:
      s.scattering_type_registry.process(f)
    s.scattering_type_registry.assign_from_table(form_factor_table)
    s.n_electrons = s.scattering_type_registry.\
                    sum_of_scattering_factors_at_diffraction_angle_0()

    # apply solvent model
    if (self.use_solvent):
      sm = solvent_model()
      sm.interpreted_pdb = interpreted_pdb
      sm.xyz = s.xyz
      sm.fudge_factor = 0.6
      s.scattering_type_registry = sm.add_bulk_solvent(s.scattering_type_registry)
      s.scattering_types = sm.scattering_types
      sm.boundary_layer_scale = 0.6
      s.boundary_layer_scaling_factors = sm.add_boundary_layer_solvent\
                                         (s.scattering_type_registry)
    else:
      s.boundary_layer_scaling_factors = flex.double(len(s.xyz),0.0)

    # finalize entry
    self.species.append(s)

  def randomize_rotations(self):
    hypersphere_point = flex.double(4)
    self.rotations = [ None for i in xrange(len(self.species)) ]
    for i in xrange(len(self.species)):
      self.rotations[i] = [ None for j in xrange(self.species[i].n_copies) ]
      for j in xrange(self.species[i].n_copies):

        # pick random point on a 3-sphere
        for k in xrange(hypersphere_point.size()):
          hypersphere_point[k] = self.random.normalvariate(0.0,1.0)
        hypersphere_point = hypersphere_point /\
                            math.sqrt(flex.sum_sq(hypersphere_point))
        a = hypersphere_point[0]
        b = hypersphere_point[1]
        c = hypersphere_point[2]
        d = hypersphere_point[3]

        # convert quaternion to rotation matrix
        self.rotations[i][j] =\
          (a*a + b*b - c*c - d*d,     2*b*c - 2*a*d,         2*b*d + 2*a*c,
               2*b*c + 2*a*d,     a*a - b*b + c*c - d*d,     2*c*d - 2*a*b,
               2*b*d - 2*a*c,         2*c*d + 2*a*b    , a*a - b*b - c*c + d*d)

  def valid(self,species=None,translation=None):
    # since the reference structures are all at the origin, each translation
    # is actually the center of each molecule in the box
    # check that molecule is in the box
    min_buffer = -0.5*self.box_size + self.species[species].radius
    max_buffer = 0.5*self.box_size - self.species[species].radius
    for ti in translation:
      if ( (ti < min_buffer) or (ti > max_buffer) ):
        return False

    # compute grid location for translation
    # assumes min_separation >> radius of largest species
    gridpoint = ( int(math.floor(translation[0]/self.gridsize)),
                  int(math.floor(translation[1]/self.gridsize)),
                  int(math.floor(translation[2]/self.gridsize)) )

    # check for overlap in boxes around grid location
    t = flex.double(translation)
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          test_gridpoint = (gridpoint[0]+i, gridpoint[1]+j, gridpoint[2]+k)
          if (self.neighbors.has_key(test_gridpoint)):
            for test_point in self.neighbors[test_gridpoint]:
              min_distance = self.species[species].radius +\
                             self.min_separation +\
                             self.species[test_point[0]].radius
              norm = (t - test_point[1]).norm()
              if (norm < min_distance):
                return False

    # add to neighbor list if translation is valid
    if (self.neighbors.has_key(gridpoint) is False):
      self.neighbors[gridpoint] = list()
    self.neighbors[gridpoint].append((species,t))

    return True

  def randomize_translations(self):
    # reset neighbor list
    self.neighbors = dict()

    # define random translation vectors
    self.translations = [ list() for i in xrange(len(self.species)) ]
    shift = 0.5*self.box_size
    for i in xrange(len(self.species)):
      for j in xrange(self.species[i].n_copies):
        current_attempt = 0
        while (current_attempt < self.max_attempts):
          translation = [self.box_size*self.random.random() - shift,
                         self.box_size*self.random.random() - shift,
                         self.box_size*self.random.random() - shift]
          if (self.valid(species=i,translation=translation)):
            self.translations[i].append(translation)
            break
          current_attempt += 1
        if (current_attempt == self.max_attempts):
          raise Sorry('unable to place molecule in box without overlap')

  def randomize(self):
    max_radius = 0.0
    for i in xrange(len(self.species)):
      if (self.species[i].radius > max_radius):
        max_radius = self.species[i].radius
    self.gridsize = self.min_separation + 2*max_radius
    self.randomize_rotations()
    self.randomize_translations()
    self.total_electrons = 0
    for s in self.species:
      self.total_electrons += s.n_electrons * s.n_copies
