import random

from iotbx import pdb
from libtbx.test_utils import approx_equal
from mmtbx.monomer_library import server, pdb_interpretation
from scitbx import matrix
from scitbx.array_family import flex

from sastbx.fXS import ewald_sphere, direct_sum_structure_factors,\
     bilinear_interpolation, nearest_neighbors, multiple_poisson,\
     solvent_accessible_area
from sastbx.fXS.solvent_model import solvent_model
from sastbx.fXS.structure_generator import structure_generator

test_pdb = "\
CRYST1   90.000   90.000   90.000  90.00  90.00  90.00 P 1           8\n\
ATOM      1  N   LYS A   1       3.287  10.092  10.329  1.00  5.89           N\n\
ATOM      2  CA  LYS A   1       2.445  10.457   9.182  1.00  8.16           C\n\
ATOM      3  C   LYS A   1       2.500  11.978   9.038  1.00  8.04           C\n\
ATOM      4  O   LYS A   1       2.588  12.719  10.041  1.00  7.07           O\n\
ATOM      5  CB  LYS A   1       1.006   9.995   9.385  1.00  3.88           C\n\
ATOM      6  CG  LYS A   1       0.016  10.546   8.377  1.00  3.81           C\n\
ATOM      7  CD  LYS A   1      -1.404  10.093   8.699  1.00  2.91           C\n\
ATOM      8  CE  LYS A   1      -2.269  10.030   7.451  1.00  3.97           C\n\
ATOM      9  NZ  LYS A   1      -3.559   9.362   7.735  1.00  2.08           N\n\
END\n\
"

# =============================================================================
def test_nearest_neighbors():
  xyz = flex.vec3_double(2)
  xyz[0] = (0.0,0.0,0.0)
  xyz[1] = (5.0,5.0,5.0)
  radius = flex.double(2)
  radius[0] = 1.0
  radius[1] = 1.0
  probe_radius = 1.0
  neighbors = nearest_neighbors(xyz,radius,0,probe_radius)
  assert(len(neighbors) == 0)

  xyz[1] = (2.0,2.0,2.0)
  neighbors = nearest_neighbors(xyz,radius,0,probe_radius)
  assert(len(neighbors) == 1)

# =============================================================================
def test_solvent_accessible_area():
  xyz = flex.vec3_double(2)
  xyz[0] = (0.0,0.0,0.0)
  xyz[1] = (5.0,5.0,5.0)
  radius = flex.double(2)
  radius[0] = 1.0
  radius[1] = 1.0
  probe_radius = 1.0
  indices = flex.int(1)
  indices[0] = 0
  n_points = 1000
  sas = solvent_accessible_area(xyz,radius,indices,probe_radius,n_points)
  assert(approx_equal(sas[0],n_points))

  xyz[1] = (1.0,1.0,1.0)
  sas = solvent_accessible_area(xyz,radius,indices,probe_radius,n_points)
  assert(approx_equal(sas[0],714))

# =============================================================================
def test_bilinear_interpolation():
  f = flex.vec3_double(4)
  f[0] = (0,0,10)
  f[1] = (10,0,10)
  f[2] = (0,10,20)
  f[3] = (10,10,20)

  for i in xrange(len(f)):
    assert(approx_equal(bilinear_interpolation(f[i][0:2],f),f[i][2]))
  assert(approx_equal(bilinear_interpolation((5,5),f),15.0))

# =============================================================================
def test_multiple_poisson():
  maxint = 2147483647
  t = 10.0
  l = flex.double(1000000,t)
  p = multiple_poisson(l,random.randint(0,maxint)).as_double()
  m = flex.mean(p)
  v = flex.mean(p*p) - m*m
  assert approx_equal(m,t,eps=0.05)
  assert approx_equal(v,t,eps=0.05)

# =============================================================================
def test_direct_sum():

  # correct values
  p = pdb.input(source_info='string',lines=test_pdb)
  x = p.xray_structure_simple()
  for s in x.scatterers():
    s.set_use_u(False,False)
  fc = x.structure_factors(anomalous_flag=False,d_min=2.0,algorithm='direct').f_calc()
  fcd = fc.data()
  indices = fc.indices()

  # test values
  xyz = x.sites_cart()
  h = flex.vec3_double(len(indices))
  fm = matrix.sqr(p.crystal_symmetry().unit_cell().fractionalization_matrix())
  for i in xrange(len(indices)):
    h[i] = fm * indices[i]
  sr = x.scattering_type_registry()
  st = x.scattering_types()
  bls = flex.double(len(xyz),0.0)
  amplitudes = direct_sum_structure_factors(st,xyz,bls,h,sr)

  for i in xrange(len(indices)):
    assert(approx_equal(fcd[i],amplitudes[i]))

# =============================================================================
def test_solvent_model():

  # correct values
  p = pdb.input(source_info='string',lines=test_pdb)
  sr = p.xray_structure_simple().scattering_type_registry()

  # test values
  p1 = pdb.input(source_info='string',lines=test_pdb)
  mls = server.server()
  el = server.ener_lib()
  ip = pdb_interpretation.process(mon_lib_srv=mls,ener_lib=el,pdb_inp=p1)
  sg = structure_generator()
  sg.add_species(p,1)
  sm = solvent_model()
  sm.interpreted_pdb = ip
  sm.xyz = sg.species[0].xyz

  new_sr = sm.add_bulk_solvent(sr)
  assert(new_sr.sum_of_scattering_factors_at_diffraction_angle_0() <
         sr.sum_of_scattering_factors_at_diffraction_angle_0())
  sm.bulk_solvent_scale = 0.0
  new_sr = sm.add_bulk_solvent(sr)
  assert(approx_equal(sr.sum_of_scattering_factors_at_diffraction_angle_0(),
                      new_sr.sum_of_scattering_factors_at_diffraction_angle_0()))

  sm.boundary_layer_scale = 0.0
  new_sr = sm.add_boundary_layer_solvent(new_sr)
  assert(approx_equal(sr.sum_of_scattering_factors_at_diffraction_angle_0(),
                      new_sr.sum_of_scattering_factors_at_diffraction_angle_0()))
  sm.boundary_layer_scale = 0.6
  new_sr = sm.add_boundary_layer_solvent(new_sr)
  assert(new_sr.sum_of_scattering_factors_at_diffraction_angle_0() >
         sr.sum_of_scattering_factors_at_diffraction_angle_0())

# =============================================================================
def test_ewald_sphere():
  es = ewald_sphere()
  es.set_wavelength(1.0)
  es.set_distance(10.0)
  xy = (0.0,0.0)
  h = es.get_h(xy)
  q = es.h_to_q(h)
  assert(approx_equal(h,(0.0,0.0,0.0)))
  assert(approx_equal(q,0.0))
  xy = (10.0,10.0)
  h = es.get_h(xy)
  q = es.h_to_q(h)
  assert(approx_equal(h,(0.57735026919, 0.57735026919, -0.42264973081)))
  assert(approx_equal(q,5.77677116966))

# =============================================================================
def test_structure_generator():
  p = pdb.input(source_info='string',lines=test_pdb)
  sg = structure_generator()
  sg.add_species(p,100)
  sg.randomize()

  t = sg.translations[0]
  d = flex.double()
  for i in xrange(len(t)):
    for j in xrange(i+1,len(t)):
      d.append( (flex.double(t[i]) - flex.double(t[j])).norm() )
  assert ( flex.min(d) > (sg.min_separation + 2.0*sg.species[0].radius) )

# =============================================================================
if (__name__ == '__main__'):
  test_nearest_neighbors()
  test_solvent_accessible_area()
  test_bilinear_interpolation()
  test_multiple_poisson()
  test_direct_sum()
  #test_solvent_model()
  test_ewald_sphere()
  test_structure_generator()
  print 'Ok'
