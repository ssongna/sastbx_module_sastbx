import sys, os, math
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import pdb
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
from libtbx import easy_pickle
from cctbx.eltbx import xray_scattering
from sastbx import intensity
from sastbx.data_reduction import saxs_read_write
from mmtbx.monomer_library import server, pdb_interpretation

def read_dummy_type(file_name=None,pdb_inp=None):
  assert ((file_name is not None) or (pdb_inp is not None))
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()

  if (pdb_inp is None):
    tmp_obj = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=file_name).all_chain_proxies
  else:
    tmp_obj = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      pdb_inp=pdb_inp).all_chain_proxies    

  ats = tmp_obj.nonbonded_energy_type_registry.symbols
  els = tmp_obj.scattering_type_registry.symbols

  for i in range(ats.size()):
    if ats[i]=="":
      ats[i]=els[i]

  return ats

def scattering_factor(a,b,c,q):
  result = q*0
  stol = q/(math.pi*4)
  for aa,bb in zip(a,b):
    result+= aa*flex.exp(-bb*stol*stol)
  result += c
  return result

def dummy_factor(v,q):
  result = q*0
  v2_3 = v**(2.0/3.0)
  stol = q/(math.pi*4.0)
  result = v * flex.exp(-math.pi*stol*stol*v2_3)
  return result

def load_scaling_factor():
  Scaling_factors={}
  #Carbon sp3
  Scaling_factors['CH1']=([1.16648272088 ,0.0274345176965 ])
  Scaling_factors['CH2']=([1.33301049942 ,0.06681996399 ])
  Scaling_factors['CH3']=([1.49986734604 ,-0.0154123735207])

  Scaling_factors['NH1']=([1.14278088926 ,-0.0726165370032 ])
  Scaling_factors['NH2']=([1.28593051297 ,0.0958307550963 ])
  Scaling_factors['NH3']=([1.42827623616,0.0615917265355])

  Scaling_factors['OH1']=([1.12467042278 ,0.0212677806006 ])
  Scaling_factors['SH1']=([1.06219008731 ,0.0172765796161 ])
  #carbon sp2 (may not be the same as sp3, but using sp3 for now!)
  Scaling_factors['C1']=Scaling_factors['CH1']
  Scaling_factors['C2']=Scaling_factors['CH2']
  Scaling_factors['CR1H']=Scaling_factors['CH1']
  Scaling_factors['CR15']=Scaling_factors['CH1']
  Scaling_factors['CR16']=Scaling_factors['CH1']

  Scaling_factors['CSP1']=Scaling_factors['CH1']

  Scaling_factors['NS1']=Scaling_factors['NH1']
  Scaling_factors['NC1']=Scaling_factors['NH1']
  Scaling_factors['NC2']=Scaling_factors['NH2']
  Scaling_factors['NR15']=Scaling_factors['NH1']
  Scaling_factors['NR16']=Scaling_factors['NH1']
  Scaling_factors['NT1']=Scaling_factors['NH1']
  Scaling_factors['NT2']=Scaling_factors['NH2']
  Scaling_factors['NT3']=Scaling_factors['NH3']
  return Scaling_factors

def build_scattering_library(atom_types, q_array, radii, radius_scale, Explicit_H, S_factor):
  scat_lib = intensity.scattering_library( q_array )
  ener_lib = server.ener_lib()
  for at in atom_types:
    element=ener_lib.lib_atom[at].element
    if(element is ''):
      element='C'
      print "Warning: unexpected atom found, and C element is used"
    val = xray_scattering.it1992(element, True).fetch()
    a = val.array_of_a()
    b = val.array_of_b()
    c = val.c()
    sf_result = scattering_factor(a,b,c,q_array)
    # getting displaced solvent volumn
    if(radii.__contains__(at)):
       r=radii[at]*radius_scale
       v=math.pi*4.0/3.0*r**3.0
    else:
       v=16.44
    dummy_sf = dummy_factor(v,q_array)

    if(not Explicit_H):
      if(S_factor.__contains__(at)):
        scale=S_factor[at]
        sf_result *= (scale[0]*flex.exp(-scale[1]*q_array*q_array))
        dummy_sf *= (flex.exp(-1.25*q_array*q_array))

    scat_lib.load_scattering_info( at, sf_result,dummy_sf )
  return scat_lib

def calc_abs_Io( atom_types, at_collection, S_factor, Explicit_H=False ):
  ener_lib = server.ener_lib()
  abs_Io = 0.0
  for at in at_collection:
    element=ener_lib.lib_atom[at].element
    if(element is ''):
      element='C'
      print "Warning: unexpected atom found, and C element is used"
    val = xray_scattering.it1992(element, True).fetch()
    a = val.array_of_a()
    c = val.c()
    this_sf0=c
    for aa in a: this_sf0 += aa
    if(not Explicit_H):
      if(S_factor.__contains__(at)):
        scale=S_factor[at]
        this_sf0 = this_sf0 * scale[0]
    sel=flex.bool(atom_types==at)
    isel=sel.iselection()
    n_at = isel.size()
#    print at, this_sf0
    abs_Io = abs_Io + n_at*this_sf0
  abs_Io = abs_Io**2.0
  return abs_Io

def calc_abs_Io_from_pdb( filename, Explicit_H=False ):
  atom_types = read_dummy_type( file_name=filename )
  at_collection=[]
  for at in atom_types:
    if at not in at_collection:
      at_collection.append( at )
  Scaling_factors = load_scaling_factor()
  abs_Io = calc_abs_Io( atom_types, at_collection, Scaling_factors, Explicit_H=Explicit_H )
  return abs_Io
