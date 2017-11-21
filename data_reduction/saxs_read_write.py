import os,sys
from scitbx.array_family import flex
import curves

def read_standard_ascii_qis(file_name, sigma_simul=0.08):
   separators = [ '', ' ', ',',';', '->','&']
   file = open(file_name,'r')
   q = flex.double()
   i = flex.double()
   s = flex.double()
   comments = []
   for line in file:
    all_good=False
    try:
     if line[0]!="#":
       keys=line.split("\n")[0].split()
       new_keys = []
       for key in keys:
         if key not in separators:
           new_keys.append( key )
       keys = new_keys
       q.append( float(keys[0]) )
       i.append( float(keys[1]) )
       if len(keys)==3:
         s.append( float(keys[2]) )
       else:
         s.append( float(keys[1])*sigma_simul )
     else:
       comments.append( line[0:len(line)-1] )
     all_good=True
    except: pass
    if not all_good:
      print "WARNING TROUBLE READING THIS LINE:"
      print line
      print "-------"

   this_data = curves.simple_saxs_data(q,i,s,comments)
   return this_data



def write_standard_ascii_qis(data,file_name=None):
  set_ext=False
  if file_name is None:
    file_name = data.data_id
    set_ext=True
  if file_name is None:
    file_name = "did_you_forget_to_give_a_sample_id_in_your_input_file"
  if set_ext:
    file_name=file_name +".qis"
  output = open(file_name,'w')
  for q,i,s in zip(data.q, data.i, data.s):
    print >> output, "%8.3e    %8.3e    %8.3e"%(q,i,s)




def read_array_of_saxs_data(file_name_array):
  result = []
  for file_name in file_name_array:
    saxs_data = read_standard_ascii_qis( file_name )
    result.append( saxs_data )


def write_array_of_saxs_data(data_array):
  for data in data_array:
    write_standard_ascii_qis(data)
