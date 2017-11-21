from scitbx.array_family import flex
import sys,os,math

class simple_saxs_data(object):
  def __init__(self, q,i,s,comments=[],concentration_mg_ml=None, data_type=None, data_id=None, cleanup=True):
    self.q = q
    self.i = i
    self.s = s
    self.comments = comments
    self.concentration=concentration_mg_ml
    self.data_type = data_type
    self.data_id = data_id
    self.clean_data()


  def multiply_add(self,scale,offset=0):
    self.i = self.i*scale+offset
    self.s = self.s*scale

  def deep_copy(self):
    nc = list( self.comments )
    new_curve = simple_saxs_data( self.q.deep_copy(),self.i.deep_copy(),self.s.deep_copy(),nc,self.concentration,self.data_type, self.data_id)
    return new_curve

  def add_curve(self, external_data, factor=-1):
    self.i = self.i+factor*external_data.i
    self.s = self.s*self.s+factor*factor*external_data.s*external_data.s
    self.s = flex.sqrt( self.s )
    for line in external_data.comments:
      self.comments.append( line )

  def clean_data(self,eps=1e-18):
    tmp_s = flex.double()
    tmp_i = flex.double()
    tmp_q = flex.double()
    for qq,ii,ss in zip(self.q, self.i, self.s):
      not_okai=False
      if abs(ii-ss)<=eps:
        if abs(ii)<=eps:
          not_okai = True
      if not not_okai:
        tmp_q.append( qq )
        tmp_i.append( ii )
        tmp_s.append( ss )
    self.i = tmp_i
    self.s = tmp_s
    self.q = tmp_q

  def show_summary(self, out=None,nn=50):
    if out is None:
      out = sys.stdout
    qq = self.q.size()
    if qq<=nn:
      nn = int(qq)/2
    nn = int(nn)
    print >> out, "Data id         : ", self.data_id
    print >> out, "Data type       : ", self.data_type
    print >> out, "Concentratation :  %s"%(str(self.concentration))
    print >> out, "q-range         :  %5.3e   %5.3e"%(flex.min(self.q),flex.max(self.q))
    print >> out, "N obs           :  %5.3e"%(self.q.size())
    print >> out, "<I/sigI> low    :  %5.3e"%( flex.mean( self.i[0:nn]/(self.s[0:nn]+1e-13) ) )
    print >> out, "<I/sigI> high   :  %5.3e"%( flex.mean( self.i[qq-nn:]/(self.s[qq-nn:]+1e-13) ) )
    print >> out, "    low q range :  %5.3e  %5.3e"%( self.q[0], self.q[nn] )
    print >> out, "   high q range :  %5.3e  %5.3e"%( self.q[qq-nn], self.q[qq-1] )
    print >> out
    for cmm in self.comments:
      print >> out, cmm
    print >> out

  def cut_data(self, q_min, q_max):
    selection_low  = flex.bool(self.q >= q_min)
    selection_high = flex.bool(self.q <= q_max)
    selection = selection_low & selection_high
    new_q = self.q.select( selection )
    new_i = self.i.select( selection )
    new_s = self.s.select( selection )
    new_object = simple_saxs_data(new_q,
                                  new_i,
				  new_s,
                                  self.comments,
				  self.concentration,
				  self.data_type, 
				  self.data_id)
    return new_object 


def background_subtract(sample,buffer,sv=7.3e-4,conc=None):
  assert (sample.concentration is not None) or (conc is not None)
  new_curve = sample.deep_copy()
  if sample.concentration is None:
    new_curve.concentration = conc

  new_curve.multiply_add( 1,0) #sv*new_curve.concentration,0 )
  tmp_curve = buffer.deep_copy()
  tmp_curve.multiply_add( (1.0-sv*new_curve.concentration),0 )
  new_curve.add_curve(tmp_curve,factor=1.0)
  return new_curve
