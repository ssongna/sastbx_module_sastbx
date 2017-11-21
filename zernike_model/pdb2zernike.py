import sys, os
from stdlib import math as smath
from scitbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import pdb, ccp4_map
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from cctbx.eltbx import xray_scattering
from sastbx.data_reduction import saxs_read_write
import time
from mmtbx.monomer_library import server, pdb_interpretation

from scitbx import math
from sastbx import zernike_model as zm
from iotbx.xplor import map as xplor_map
from cctbx import uctbx, sgtbx

from libtbx import easy_pickle
from sastbx.interface import get_input


master_params = iotbx.phil.parse("""\
zernike{
  pdbfile = None
  .type=path
  .help="the pdb file"

  qmax = 0.3
  .type=float
  .help="maximum q value, for which the intensity to be evaluated"

  nmax=20
  .type=int
  .help="maximum order of zernike expansion"

  np=50
  .type=int
  .help="number of point covering [0,1]"

  fix_dx = True
  .type=bool
  .help="Whether the dx will be adjusted to match the size of molecule"

  buildmap=False
  .type=bool
  .help="Whether xplor map will be constructed or not"

  shift = False
  .type=bool
  .help="Whether the pdb coordates will be shifted or not"

  coef_out=True
  .type=bool
  .help="Whether dump zernike moments to picle files"
}
""")


banner = "--------------Zernike Moments Calculation and Map Construction----------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "Usage: libtbx.python pdb2zernike.py pdbfile=pdbfile nmax=nmax buildmap=True/False shift=True/False np=np_of_point"

def ccp4_map_type(map, N, radius,file_name='map.ccp4'):
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  map.reshape( grid )
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90"),
      space_group=sgtbx.space_group_info("P1").group(),
      gridding_first=(0,0,0),
      gridding_last=(N*2, N*2, N*2),
      map_data=map,
      labels=flex.std_string(["generated from zernike moments"]))



def xplor_map_type(m,N,radius,file_name='map.xplor'):
  gridding = xplor_map.gridding( [N*2+1]*3, [0]*3, [2*N]*3)
  grid = flex.grid(N*2+1, N*2+1,N*2+1)
  m.reshape( grid )
  uc = uctbx.unit_cell(" %s"%(radius*2.0)*3+"90 90 90")
  xplor_map.writer( file_name, ['no title lines'],uc, gridding,m) # is_p1_cell=True)  # True)


def zernike_moments(pdbfile, nmax=20, np=50, fix_dx=False, np_on_grid=20, shift=False, buildmap=False, coef_out=True, calc_intensity=True, external_rmax=-1):
  base = pdbfile.split('.')[0]
  splat_range = 0
  fraction = 0.9
  default_dx = 0.7
  uniform = True
  #iotbx/pdb/hierarchy.py-->class input
  pdbi = pdb.hierarchy.input(file_name=pdbfile)
  if(len( pdbi.hierarchy.models() ) == 0):
    return None,None,None

  atoms = pdbi.hierarchy.models()[0].atoms()
  # predefine some arrays we will need
  atom_types = flex.std_string()
  radius= flex.double()
  b_values = flex.double()
  occs = flex.double()
  xyz = flex.vec3_double()
  # keep track of the atom types we have encountered
  for atom in atoms:
    if(not atom.hetero):
      xyz.append( atom.xyz )
#    b_values.append( atom.b )
#    occs.append( atom.occ )
  #print time.time()

  if(xyz.size() == 0):
    return None,None,None
  density=flex.double(xyz.size(),1.0)
  voxel_obj = math.sphere_voxel(np,splat_range,uniform,fix_dx,external_rmax, default_dx, fraction,xyz,density)
  np = voxel_obj.np()
  rmax=voxel_obj.rmax()/fraction

  #print base, "RMAX: ", voxel_obj.rmax()

  #print time.time()
  grid_obj = math.sphere_grid(np, nmax)
  pdb_out = False
  grid_obj.clean_space(voxel_obj, pdb_out)
  #print time.time(), "END GRID CONST"
  grid_obj.construct_space_sum()

  #print time.time(), "SS CALC"
  mom_obj = math.zernike_moments(grid_obj, nmax)
  #print time.time(), "END MOM CALC"

  moments = mom_obj.moments()
  if(coef_out):
    nn = mom_obj.fnn()
    easy_pickle.dump(base+'.nlm.pickle', moments.coefs() )
    easy_pickle.dump(base+'.nn.pickle', nn.coefs() )

  #
####### The following will be optional ###
  if(shift):
    shift = [rmax, rmax, rmax]
    centered_xyz = voxel_obj.xyz() + shift
    out_pdb_name=base+'_centered.pdb'
    for a,xyz in zip( atoms, centered_xyz):
      a.set_xyz( new_xyz=xyz)
    pdbi.hierarchy.write_pdb_file( file_name=out_pdb_name, open_append=False)

    original_map = voxel_obj.map()
    xplor_map_type( original_map, np, rmax, file_name=base+'_pdb.xplor')


######## Calculate Intnesity Profile ############
  if(calc_intensity):
    q_array = flex.double( range(51) )/100.0
    z_model = zm.zernike_model( moments, q_array, rmax, nmax)
    nn = mom_obj.fnn()
    intensity = z_model.calc_intensity(nn)
    #iq_file = open(base+"_"+str(nmax)+"_"+str(int(rmax*fraction))+".zi", 'w')
    iq_file = open(base+"_"+str(nmax)+".zi", 'w')
    for qq, ii in zip(q_array, intensity):
      print>>iq_file, qq,ii

    iq_file.close()
####### END of Intensity Calculation ###########

  if(buildmap):
    print "grid setting up..."
    zga = math.zernike_grid(np_on_grid, nmax, False)
    print "grid setting up...Done"
    zga.load_coefs(moments.nlm(), moments.coefs() )
    print "start reconstruction"
    map=flex.abs(zga.f() )
    print "finished reconstruction"
    xplor_map_type( map, np_on_grid, rmax, file_name=base+'.xplor')
    ccp4_map_type( map, np_on_grid, rmax, file_name=base+'.ccp4' )

  return mom_obj, voxel_obj, pdbi

def disp_inv(inva):
  for indx, inv in zip(inva.nl(), inva.coefs() ):
    print indx, inv

def display(moments):
  for indx, mom in zip(moments.nlm(), moments.coefs()):
    print indx, mom

def run(args):
  params = get_input(args, master_params, "zernike", banner, help)
  if params is None:
    return
  pdbfile = params.zernike.pdbfile
  nmax=params.zernike.nmax
  np = params.zernike.np
  fix_dx = params.zernike.fix_dx
  shift = params.zernike.shift
  buildmap = params.zernike.buildmap
  coef_out = params.zernike.coef_out

  zernike_moments(pdbfile, nmax, np, fix_dx=fix_dx, shift=shift, buildmap=buildmap, coef_out=coef_out, calc_intensity=True)


if __name__ == "__main__":
   args=sys.argv[1:]
   run(args)
