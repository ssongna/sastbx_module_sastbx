import os,sys
from scitbx.array_family import flex
from sastbx import refine
import time
import iotbx.phil
from libtbx.utils import Sorry, date_and_time, multi_out
import libtbx.phil
import libtbx.phil.command_line
from cStringIO import StringIO

master_params = iotbx.phil.parse("""\
pr{
  pdb_file= None
  .type=path
  .help="PDB filename to be calculated"

  Max = 200
  .type=float
  .help="upper limit of the histgram building"

  bin_size = 1.0
  .type=float
  .help="bin size for the histgram"

  type = "all"
  .type=str
  .help="the p(r) of atom selection: CA | all"

  output="pdb.pr"
  .type=path
  .help="Output file name for calculated P(r)"

}

""")


class PDB(object):
  def __init__(self,file_name,scale):
    self.CA_indx=flex.int()
    self.label=[]
    self.xyz=flex.vec3_double()
    self.natm=0
    self.readPDB(file_name)
    self.crd=flex.double()
    self.eigens=[]
    for xyz in self.xyz:
      for value in xyz:
        self.crd.append(value)

  def readPDB(self, file_name):
    ff = open( file_name, 'r') 
    for line in ff:
      if(line[0:4] == 'ATOM'):
#        self.label.append(line[0:30])
        if(line[13:15] == 'CA'):
          self.CA_indx.append(self.natm)
        x=float(line[30:38])
        y=float(line[38:46])
        z=float(line[46:54])
        self.xyz.append([x,y,z])
        self.natm=self.natm+1        
  def Hessian(self,cutoff,Nmodes,scale):
    self.model=refine.elastic(self.xyz,self.CA_indx)
 #   for modes in range(Nmodes):
 #      self.eigens.append(self.model.Project2All(modes+7))

  def writePDB(self, crd, file_name):
    ff = open( file_name, 'w')
    for ii in range(self.natm):
      str=self.label[ii]+"%8.3f%8.3f%8.3f" % (crd[ii][0],crd[ii][1],crd[ii][2])
      ff.write(str+"\n")

  def NMPerturb(self, modes, weights):
    new_crd = flex.double()
    new_crd = self.crd.deep_copy()
    for m, w in zip(modes, weights):
      print m,w, new_crd.size(), self.eigens[0].size()
      new_crd += self.eigens[m-7]*w
    new_xyz = flex.vec3_double(new_crd)
    return new_xyz

def run(params,log):
  file=params.pr.pdb_file
  t1 = time.time()
  scale = 1.0
  pdb=PDB(file,scale)
  cutoff = 14
  Nmodes = 10
  pdb.Hessian(cutoff,Nmodes,scale)
  if(params.pr.type == "all"):
    pdb.model.updateDistArrayAll(pdb.xyz)
  elif(params.pr.type == "CA"):
    pdb.model.updateDistArray(pdb.xyz)
  else:
    raise Sorry("wrong type: please specify either 'all' or 'CA'")

  dist_array=pdb.model.getDistArray()
  dMax=params.pr.Max
  bin_size = params.pr.bin_size
  n_slot=int(dMax/bin_size+0.5)
  pr=pdb.model.Histogram(dist_array,dMax,n_slot)
  rs=range(n_slot)

  output=open(params.pr.output,'w')
  for r,p in zip(rs,pr):
    if(p > 0 or (r < 10)):
      print >>output, r*bin_size,p/bin_size

  t2 = time.time()
  #print '&', time.ctime(t1)
  #print '&', time.ctime(t2)

def get_input(args):
  if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil=master_params,
      home_scope="pr")

    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log,expert_level=1)
    print >> log, "#phil __END__"
    run( params, log )

def print_help():
  print "\nUsage:\n sastbx.pr pdb_file=PDB_FILE Max=MAX type=type(\"all\" or \"CA\") output=outputfile\n"

    
if __name__ == "__main__":
  get_input(sys.argv[1:])
