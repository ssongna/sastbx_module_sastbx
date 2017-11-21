from scitbx.array_family import flex
from sastbx.fXS import zernike_moment_variants as znk_blq
import sys,os

class blq_data(object):
  """This is a data container for BLQ data structure, which has:
     .q: momentum transfer values (4*pi*sin(theta)/lambda)
     .blq: 1-D flex array, containing the BLQ expansion coefs
  """
  def __init__(self,qs,blq,lmax):
    self.q = qs
    self.blq = blq
    self.lmax = lmax
  def format_blq(self,data=None):
    if(data is None):
      data=self.blq.deep_copy()
    new_blq=[]
    n_l = int(self.lmax/2)+1  # lmax is even number, by assumption
    for ii in range(0,data.size(),n_l):
      new_blq.append(data[ii:ii+n_l])
    return new_blq
  def print_out(self,data=None,out=None):
    if(data is None):
      format_blq = self.format_blq()
    else:
      format_blq = self.format_blq(data)
    if(out is None):
      out=sys.stdout
    print>>out,"# q, B0, B2, B4, ... (since Odd orders are ZERO) "
    for q,blq in zip(self.q,format_blq):
      print>>out,q,
      for cc in blq:
        print>>out,cc,
      print>>out


def read_blq(filename, lmax=None):
  comments=['#','%']
  qs= flex.double()
  blq=flex.double()
  input=open(filename,'r')
  if(lmax is not None): n_l = int(lmax/2) + 1
  for line in input:
    keys = line.split()
    if(keys[0] in(comments)):
      continue
    qs.append(float(keys[0]))
    if( (lmax is None) or (len(keys) < n_l ) ):
      n_l = len(keys)
      lmax = n_l*2 - 1
    for kk in keys[1:n_l+1]:  # Because coefs start from col-2
      blq.append(float(kk) )
  input.close()
  #lmax = int(blq.size()/qs.size()) - 1  # because l starts with 0
  data = blq_data(qs,blq,lmax)
  return data


### Testing routines for the functions defined in this file ###
def test_read_blq():
  import os
  tmp_file='tmp.blq'
  output=open(tmp_file,'w')
  print>>output,'#'
  print>>output,'%'
  qs = flex.double(range(10))/20.0
  max_L = 10
  blq = flex.double(xrange(qs.size()*max_L))
  count=0
  for q in qs:
    print>>output, q,
    for ii in xrange(max_L):
      print>>output, blq[count],
      count=count+1
    print>>output
  output.close()
  data = read_blq(tmp_file)
  new_qs = data.q
  new_blq = data.blq
  assert new_qs.size() == qs.size()
  assert new_blq.size() == blq.size()
  for q0,q1 in zip(qs,new_qs):
    assert(abs(q0-q1)<1e-6)
  for b0,b1 in zip(blq,new_blq):
    assert(abs(q0-q1)<1e-6)
  os.remove(tmp_file)

if __name__=="__main__":
  test_read_blq()
