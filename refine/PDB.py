import os,sys
from scitbx.array_family import flex
from iotbx import pdb
from sastbx import refine
import time

class PDB(object):
  def __init__(self,file_name, method='ca', n_res_in_block=3, max_n_block=300 ):
    self.method=method
    self.n_res_in_block=n_res_in_block
    self.max_n_block=max_n_block
    self.block_start=flex.int()
    self.xyz=flex.vec3_double()
    self.natm=0
    self.readPDB(file_name)
    self.eigens=[]
    self.dxyz=self.xyz*0.0
    self.r = 0.0  # ca-pair distance restraint
    if(self.method == 'ca' ):
      self.r_list=self.build_restraint_list()

  def readPDB(self, file_name):
    self.pdbi = pdb.hierarchy.input(file_name=file_name)
    if(len( self.pdbi.hierarchy.models() ) == 0):
      return None

    self.atoms = self.pdbi.hierarchy.models()[0].atoms()

    if(self.method == 'ca'):
      # keep track of the atom types we have encountered
      for atom in self.atoms:
        #if(not atom.hetero):
          self.xyz.append( atom.xyz )
          if(atom.name in [" CA "," P  "," C5 "]):
            self.block_start.append(self.natm)
          self.natm=self.natm+1
      print "CA-based Elastic Network Model will be used"
      print "number of points used for ENM: ", self.block_start.size()
      self.n_block=self.block_start.size()
      self.ca_dxyz = flex.vec3_double( self.n_block, (0.0, 0.0, 0.0) )
    else:
      self.residue_groups=self.pdbi.hierarchy.models()[0].residue_groups()
      n_atm_in_res=flex.int()
      res_start=flex.int()
      self.xyz=self.atoms.extract_xyz()
      self.natm = self.xyz.size()
      self.atoms.reset_i_seq()
      for rg in self.residue_groups:
        rg_atoms = rg.atoms()
        n_atm_in_res.append( rg_atoms.size() )
        res_start.append( rg_atoms[0].i_seq )
      n_res=res_start.size()
      ###################20170617####################
      self.n_res_in_block=max(n_res/self.max_n_block, self.n_res_in_block)
      self.n_res_in_block = int(self.n_res_in_block)
      ### Merge the residues to block ###
      self.n_atm_in_block=flex.int()
     

      for ii in range(0,n_res,self.n_res_in_block):
      #  merge_n_atm = flex.sum( n_atm_in_res[ii:ii+self.n_res_in_block] )
      #  self.n_atm_in_block.append(merge_n_atm)
        self.block_start.append( res_start[ii] )
      self.n_block=self.block_start.size()
      ### is there any remaining residues ###
      nres_left = n_res - self.n_block*self.n_res_in_block
      if(nres_left>0):
      #  merge_n_atm = flex.sum(n_atm_in_res[-res_left:])
      #  self.n_atm_in_block.append(merge_n_atm)
        self.block_start.append(res_start[-res_left])
        self.n_block +=1
      self.ca_dxyz = flex.vec3_double( self.natm, (0.0, 0.0, 0.0) )  ##Because each block has 3*2 d.o.f.
      print "RTB-based Elastic Network Model will be used"
      print "number of blocks used for ENM: ", self.block_start.size()



  def Hessian(self,cutoff,Nmodes,sf,percentage=0.5):
    if(self.method == 'ca'):
      self.model=refine.elastic(self.xyz,self.block_start,cutoff,sf)
    else:
      self.model=refine.elastic_rtb(self.xyz,self.block_start,cutoff,sf)
    for modes in range(Nmodes):
      self.eigens.append(self.model.mode(modes+6))
    suggested_nmode = self.model.nmode(percentage, Nmodes+100)
    #print "Based on the Eigenvalues, the top %d modes cover %% %d motions"%(suggested_nmode, percentage*100)
    return suggested_nmode

  def eigenvalues(self):
    return self.model.eigenvalues()

  def writePDB(self, crd, file_name):
    self.atoms.set_xyz( crd )
    self.pdbi.hierarchy.write_pdb_file( file_name=file_name, open_append=False)

  def NMPerturb(self, modes, weights):
    self.dxyz = self.dxyz*0.0
    for m, w in zip(modes, weights):
      self.dxyz += self.eigens[m-7]*w
    return self.xyz+self.dxyz


  def build_restraint_list(self):
    n_restraint = self.n_block*0
    r_a = flex.random_size_t(n_restraint)%self.n_block
    r_b = flex.random_size_t(n_restraint)%self.n_block
    bead_pairs = flex.tiny_size_t_2()
    for a,b in zip(r_a, r_b):
      bead_pairs.append( (a,b) )
    return bead_pairs
