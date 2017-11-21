from scitbx.array_family import flex
from iotbx import pdb
from scitbx import math, fftpack
import os, sys
from stdlib import math as smath


class projection(object):
  """ This class takes pdb file (xyz) and project it to 2D plane
      The projection will be along given vector
  """
  def __init__(self, pdb_file, splat_range=1, dx=1.0, external_rmax=-1, fraction=0.3):
    self.pdb_file = pdb_file
    self.pdbi = pdb.hierarchy.input( file_name=pdb_file )
    self.atoms = self.pdbi.hierarchy.atoms()
    self.xyz = self.atoms.extract_xyz().deep_copy()

    self.splat_range=splat_range
    self.dx = dx
    self.external_rmax= external_rmax
    self.fraction= fraction

    self.voxel = None


  def rotate(self, vector):
    rotation_matrix = math.r3_rotation_vector_to_001( vector )
    self.new_xyz = self.xyz*rotation_matrix

  def mapping(self, xyz):
    self.voxel = math.two_d_voxel(self.splat_range, self.external_rmax, self.dx, self.fraction, xyz )
    return self.voxel.get_image()

  def project(self, vector=None):
    if(vector is not None):
      self.rotate( vector )
      xyz = self.new_xyz
    else:
      xyz = self.xyz
    projection = self.mapping(xyz)
    return projection

  def get_value(self):
    if( self.voxel is None ): self.project()
    return self.voxel.get_value()

  def get_np(self):
    if( self.voxel is None ): self.project()
    return self.voxel.np()*2+1
 
  def get_dq(self):
    if( self.voxel is None ): self.project()
    return 1.0/(self.dx*self.voxel.np())*smath.pi
    

def tst_projection(pdbfile):
  projection_obj = projection( pdbfile )
  image = projection_obj.project()
  vector=[0,1,0]
  new_image = projection_obj.project( vector )

  output=open('prj1.dat', 'w')
  for pt in image:
    print >> output, pt[0],pt[1],pt[2]
  output.close()

  output=open('prj2.dat', 'w')
  for pt in new_image:
    print >> output, pt[0],pt[1],pt[2]
  output.close()

def tst_fft_proj(pdbfile):
  projection_obj = projection( pdbfile )
  image = projection_obj.project()
  output=open('prj1.dat', 'w')
  for pt in image:
    print >> output, pt[0],pt[1],pt[2]
  output.close()

  values=projection_obj.get_value()
  np = projection_obj.get_np()

  flex_grid=flex.grid(np,np)
  fft_input=flex.complex_double(flex_grid)
  for ii in range( fft_input.size() ):
    fft_input[ii] = complex( values[ii],0 )
  fft_obj = fftpack.complex_to_complex_2d( (np,np) )
  result = fft_obj.forward( fft_input )

  for ii in range(np):
    for jj in range(np):
      kk = ii
      ll = jj
      if kk > np/2: kk=kk-np
      if ll > np/2: ll=ll-np

      print kk+np/2,ll+np/2, abs( result[(ii,jj)] )


if __name__ == "__main__":
  pdbfile=sys.argv[1]
  tst_fft_proj( pdbfile )
  #tst_projection( pdbfile )
