import os,sys,math
from sastbx.pr import pr_tools
from sastbx.rigid_body import rigidbody as rb
from scitbx.array_family import flex

class rb_engine(object):
  def __init__(self, rbs, n_slots, q_array=None, she_obj=None):
    self.rbs=[]
    self.contact_list=flex.vec2_double()
    self.intra_pr=flex.double(n_slots,0.0)
    self.pr=flex.double(n_slots,0.0)
    self.norm_pr=flex.double(n_slots,0.0)
    if q_array is None:
      self.q_array = flex.double( range(51) ) /100.0
    else:
      self.q_array = q_array
    
    self.she_obj = she_obj
    self.new_coord = flex.vec3_double()
    self.new_index= flex.int()

    self.total_n_pair=0
    self.clash = 0
    for rb in rbs:
      self.rbs.append(rb)
    self.nbody=len(rbs)
    self.integrator = pr_tools.fast_integrator(n_slots, self.q_array, div=1)
    self.setup_intra_hist()

  def setup_intra_hist(self):
    for rb in self.rbs:
      self.intra_pr += rb.get_hist()

  def calc_inter_hist(self):
    self.inter_pr=self.intra_pr*0
    self.var=self.intra_pr*0
    self.clash = 0
    for ii in range(self.nbody):
      for jj in range(ii):
        self.inter_pr += self.rbs[ii].that_hist(self.rbs[jj].get_crd() )
#	self.var += self.rbs[ii].get_hist_var()
	self.clash += self.rbs[ii].get_clash()

  def get_norm_pr(self):
    self.calc_inter_hist()
    #print list(self.inter_pr)
    self.pr = self.intra_pr + self.inter_pr
    if(self.total_n_pair==0):
      self.total_n_pair=sum(self.pr)
    self.norm_pr = self.pr / self.total_n_pair
    return self.norm_pr

  def get_norm_pr_with_var(self):
    self.calc_inter_hist()
    self.pr = self.intra_pr + self.inter_pr
    if(self.total_n_pair==0):
      self.total_n_pair=sum(self.pr)
    self.norm_pr = self.pr / self.total_n_pair
    self.var = flex.pow( self.var, 0.5)
    self.var = self.var/self.total_n_pair
    return self.norm_pr, self.var


  def get_norm_intra_pr(self):
    self.norm_intra_pr = self.intra_pr / self.total_n_pair
    return self.norm_intra_pr

  def get_norm_inter_pr(self):
  #  self.calc_inter_hist()
    self.norm_inter_pr = self.inter_pr / self.total_n_pair
    return self.norm_inter_pr


  def get_clash(self):
    return self.clash

  def get_intensity(self):
    if self.she_obj is not None:
      self.new_coord.clear()
      for ii in range(self.nbody):
        self.new_coord = self.new_coord.concatenate( self.rbs[ii].get_crd()) 
      if(self.new_index.size() == 0):
        self.new_index = flex.int( range( self.new_coord.size() ) )
#        self.new_index = flex.int( range( self.rbs[1].get_crd().size() )) + self.rbs[0].get_crd().size()
      
      self.she_obj.engine.update_coord( self.new_coord, self.new_index)
      i = self.she_obj.engine.I()
    else:
      i = self.integrator.get_intensity_from_pr_array( self.get_norm_pr()  )
    return i
#    return i/i[0]

  def update_contact_list(self):
    self.contact_list = self.rbs[0].update_contact_list(self.rbs[1].get_crd() )

  def get_restraint(self):
    # for two-body case now, testing purpose
    h = self.rbs[0].calc_contact_hist( self.rbs[1].get_crd(), self.contact_list )
    return h


  def rotate_translate(self,t,r,index):
    if(len(r) == 3):
      self.rbs[index].rotate_translate(list(t),list(r*0.05) )
    else:
      self.rbs[index].rotate_translatev(list(t),list(r[0:3]),r[3] )
      #self.rbs[index].rotate_translate(list(t),r ) # quaternion


class PDB(object):
  def __init__(self,file_name,scale=10):
    self.CA_indx=flex.int()
    self.label=[]
    self.xyz=flex.vec3_double()
    self.natm=0
    self.readPDB(file_name)

  def readPDB(self, file_name):
    ff = open( file_name, 'r')
    for line in ff:
      if(line[0:4] == 'ATOM'):
        self.label.append(line[0:30])
        if(line[13:15] == 'CA'):
          self.CA_indx.append(self.natm)
        x=float(line[30:38])
        y=float(line[38:46])
        z=float(line[46:54])
        self.xyz.append([x,y,z])
        self.natm=self.natm+1

  def writePDB(self, crd, file_name):
    ff = open( file_name, 'w')
    for ii in range(len(crd)):
      str=self.label[ self.CA_indx[ii] ]+"%8.3f%8.3f%8.3f" % (crd[ii][0],crd[ii][1],crd[ii][2])
      ff.write(str+"\n")
    ff.close()

def print_pr(r,pr):
  for rr,pp in zip(r,pr):
    if(pp > 0):
      print rr,pp

def test(args):
  rbs=[]
  dmax=300
  pdb=None
  for arg in args:
    pdb= PDB(arg)
    rbs.append(rb(pdb.xyz, pdb.CA_indx, dmax))   

  rb_eng = rb_engine(rbs,int(dmax+0.5) )
  rb_eng.rotate_translate([0,0,0],[math.pi,0,0],1)
  crd = rb_eng.rbs[1].get_crd()
  pdb.writePDB(crd,'rotated.pdb') 
  pr = rb_eng.get_norm_pr()
  r = range(dmax+1)
  print_pr(r,pr)

if __name__ == "__main__":
  test(sys.argv[1:])
