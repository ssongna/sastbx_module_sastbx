from libtbx import easy_pickle
from stdlib import math as smath
from scitbx.array_family import flex
import os, sys
from sastbx.data_reduction import saxs_read_write
from sastbx.zernike_model import zernike_model, pdb2zernike
from sastbx.basic_analysis import guinier_analyses
from scitbx import math, matrix
from iotbx import pdb

import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
import libtbx.phil
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.math import zernike_align_fft as fft_align
from iotbx.xplor import map as xplor_map
from cctbx import uctbx
from scitbx.golden_section_search import gss
import time
from cctbx.macro_mol import rotation_parameters

from sastbx.interface import get_input


master_params = iotbx.phil.parse("""\
model{
  pdb_file = None
  .type=path
  .help="model in PDB format"

  nlm_file = None
  .type=path
  .help="model in zernike moments (coefficients) format"

  map_file = None
  .type=path
  .help="model in X-plor map format"

  rmax = None
  .type=float
  .help="maxium radial distance to the C.o.M."

  nmax = 20
  .type=int
  .help="maximum order of zernike polynomial:fixed for the existing database"

  output = "converted"
  .type=path
  .help="output file name"

}
""")
banner = "-------------------Model Converter (PDB/Zernike/Xplor)-------------------"

def build_3d_grid( np_on_grid, rmax_over_fraction ):
  grid = flex.vec3_double()
  one_d = range( -np_on_grid, np_on_grid+1 )
  one_d_x = flex.double(one_d)/np_on_grid*rmax_over_fraction
  for x in one_d_x:
    for y in one_d_x:
      for z in one_d_x:
        grid.append([x,y,z])
  return grid


def get_moments_for_scaled_model( map, np_on_grid, grid, nmax, rmax, external_rmax ):
  ### default params for moments calculation ###
  splat_range=1
  uniform=True
  fix_dx=False
  default_dx=0.7
  fraction=0.9
  ### end of default params for moments calculation ###

  threshold = flex.max( map )/3.0
  select = flex.bool( map.as_1d()>=threshold )
  xyz = grid.select( select )
  xyz_norms=xyz.norms()
  select = flex.bool ( xyz_norms<=rmax )
  xyz = xyz.select( select )
  density=flex.double(xyz.size(),1.0)
  vox_obj = math.sphere_voxel( np_on_grid, splat_range, uniform, fix_dx, external_rmax, default_dx, fraction, xyz,density )
  grid_obj = math.sphere_grid( np_on_grid, nmax )
  grid_obj.clean_space(vox_obj, False)
  grid_obj.construct_space_sum()
  moment_obj = math.zernike_moments( grid_obj, nmax)
  nlm_array=moment_obj.moments()
  return nlm_array.coefs()


def build_model( filename, type, nmax, rmax ):
  if( type == 'pdb' ):
    model = container( pdbfile=filename, nmax=nmax, rmax=rmax, id=filename )
  elif( type == 'nlm' ):
    nlm_coefs = easy_pickle.load( filename )
    id = filename.split('.')[0]

    model = container( nlm_coefs=nlm_coefs, nmax=nmax, rmax=rmax, id=id )
  elif( type == 'map' ):
    model = container( mapfile=filename, nmax=nmax, rmax=rmax, id=filename )
  else:
    print "model type is not identified, it should be 'pdb', 'nlm', or 'map' "
    model = None
  return model

class container(object):
  def __init__(self, pdbfile=None, nlm_coefs=None, mapfile=None, rmax=None, nmax=10, id='from_nlm', np_on_grid=30, fraction=0.9):
    self.nmax=nmax
    self.np_on_grid = np_on_grid
    self.id = id
    self.rmax=rmax
    self.vox_obj = None
    self.pdb_inp = None
    if(self.rmax is not None):
      self.rmax = self.rmax/fraction
    self.fraction=fraction    ### Attens: The Rmax is defined as length equivalent to scaled [0,1], not [0, fraction],
                              ###         Where [0, fraction] is equivalent to the Rmax of molecule
    self.map = None

    if( pdbfile is not None):
      self.build_with_pdb(pdbfile)
    elif( (nlm_coefs is not None) ):
      self.build_with_nlm(nlm_coefs)
    elif( mapfile is not None):
      self.build_with_map(mapfile)

  def build_with_pdb(self, pdbfile):
    mom_obj, vox_obj, pdb_inp = pdb2zernike.zernike_moments(pdbfile, nmax=self.nmax, fix_dx=True, coef_out=False, calc_intensity=False)
    if( mom_obj is None):
      print "The file does not have coordinate info, please check", pdbfile
    else:
      self.nlm_array = mom_obj.moments()
      if( self.rmax is None ):
        self.rmax = vox_obj.rmax() / self.fraction
      self.vox_obj = vox_obj
      self.pdb_inp = pdb_inp
      oned_np = self.vox_obj.np()*2 + 1
      self.volume = self.rmax**3.0*8.0*self.vox_obj.occupied_sites()/( oned_np**3.0 )
    self.id = pdbfile.split('.')[0]

  def build_with_nlm(self, nlm_coefs):
    self.nlm_array=math.nlm_array( self.nmax )
    self.nlm_total = self.nlm_array.coefs().size()
    if(nlm_coefs.size() < self.nlm_total):
      print "The nmax is bigger than the nmax used to build the database"
    else:
      self.nlm_coefs = nlm_coefs[0:self.nlm_total]
      self.nlm_array.load_coefs( self.nlm_array.nlm(), self.nlm_coefs )
      if( self.rmax is None):
        self.rmax=50
        print "WARNING: rmax was not specified, and default value ( 50.0 ) was used"

  def build_with_map(self, mapfile):
    this_xplor = xplor_map.reader(mapfile)
    self.np_on_grid = (this_xplor.gridding.n[0]-1 ) /2  ## this is the np_on_grid stored in xplor map
    self.np_on_grid = int(self.np_on_grid)
    self.raw_map = this_xplor.data
    self.map = flex.double( this_xplor.data.size(), 0)
    self.id = mapfile.split('.')[0]
    this_rmax = this_xplor.unit_cell.parameters()[0]/2.0
    if(self.rmax is not None and ( self.rmax != this_rmax ) ): # do the scaling
      grid=build_3d_grid(self.np_on_grid, this_rmax)
      self.nlm_coefs = get_moments_for_scaled_model( self.raw_map, self.np_on_grid, grid, self.nmax, this_rmax, self.rmax)
      self.nlm_array = math.nlm_array(self.nmax)
      self.nlm_array.load_coefs(self.nlm_array.nlm(), self.nlm_coefs )
    else:
      self.rmax = this_rmax
      print self.rmax,"test"
      threshold = flex.max( self.raw_map )/3.0
      select = flex.bool( this_xplor.data.as_1d() >= threshold )
      self.map.set_selected( select, 1)
      self.ngp = self.map.size()  # number of grid point in 3D box
      self.all_indx = flex.int(range( self.ngp ))
      self.molecule = self.all_indx.select( select )
      self.grid_obj=math.sphere_grid(self.np_on_grid, self.nmax)
      ss = self.grid_obj.construct_space_sum_via_list( self.molecule.as_1d(), self.map.as_1d() )
      self.moment_obj = math.zernike_moments( self.grid_obj, self.nmax )
      self.moment_obj.calc_moments( ss.as_1d() )
      self.nlm_array = self.moment_obj.moments()
      self.nlm_coefs = self.nlm_array.coefs()

  def write_map(self, filename='converted.map'):
    if(self.map is None):
      hex=False
      self.zga =  math.zernike_grid( self.np_on_grid, self.nmax, hex )
      self.zga.load_coefs( self.nlm_array.nlm(), self.nlm_array.coefs() )
      self.map = flex.abs( self.zga.f() )
    xplor_map_type( self.map, self.np_on_grid, self.rmax, file_name=filename )

  def write_pdb( self, rmax=None, rotation=None, filename='shifted.pdb' ):
     if(rmax is None):
       rmax = self.rmax
     shift = (rmax, rmax, rmax)
     if( rotation is not None ):
       new_xyz = self.vox_obj.rotate( (-rotation[0], rotation[1], -rotation[2]), False ) + shift
     else:
       new_xyz = self.vox_obj.xyz() + shift
     for a, xyz in zip( self.pdb_inp.hierarchy.atoms(), new_xyz ):
       a.set_xyz( new_xyz = xyz )
     self.pdb_inp.hierarchy.write_pdb_file( file_name=filename, open_append=False )

  def write_bead( self, filename='converted.pdb'):
    self.n_bds = int(self.rmax/(3.6))+1
    hex = True
    self.zga_hex = math.zernike_grid( self.n_bds, self.nmax, hex)
    self.zga_hex.load_coefs( self.nlm_array.nlm(), self.nlm_array.coefs() )
    self.map_hex = flex.abs( self.zga_hex.f() )
    max_value = flex.max( self.map_hex )
    threshold = max_value*0.3
    select = flex.bool( self.map_hex.as_1d() >= threshold )

    xyz = self.zga_hex.xyz().select( select )
    xyz = (xyz+(1.0,1.0,1.0)) * self.rmax
    res_id = 1
    out = open( filename, 'w')
    for x,y,z in xyz:
      print>>out, "ATOM  %5d  CA  ASP  %4d    %8.3f%8.3f%8.3f  1.00  1.00"%(res_id, res_id, x, y, z)
      res_id = res_id + 1
    out.close()



def xplor_map_type(m,N,radius,file_name='map.xplor'):
  gridding = xplor_map.gridding( [N*2+1]*3, [0]*3, [2*N]*3)
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  m.reshape( grid )
  uc = uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90")
  xplor_map.writer( file_name, ['no title lines'],uc, gridding,m)



def run(args):
  params = get_input(args, master_params, "model", banner, help)
  if params is None:
    return
  pdb_file = params.model.pdb_file
  map_file = params.model.map_file
  nlm_file = params.model.nlm_file
  nmax = params.model.nmax
  rmax = params.model.rmax
  output = params.model.output

  nlm_coefs = None
  if( nlm_file is not None):
    nlm_coefs = easy_pickle.load( nlm_file )

  model = container( pdbfile=pdb_file, mapfile=map_file, nlm_coefs=nlm_coefs, rmax=rmax, nmax=nmax )
  model.write_map( output+'.xplor' )
  model.write_bead(output+'.pdb' )



def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\n\tUsage: sastbx.buildmap <pdb_file=pdb_file | map_file=map_file | nlm_file=nlm_file rmax=rmax >\n"

if __name__=="__main__":
  args = sys.argv[1:]
  t1 = time.time()
  run(args)
  t2 = time.time()
  print "total time used: ", t2-t1
