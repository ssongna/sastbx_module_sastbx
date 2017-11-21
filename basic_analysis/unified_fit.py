from scitbx.array_family import flex
import scitbx.math as sm
from scitbx import differential_evolution as de
import math,sys,os




class n_level_model(object):
  def __init__(self, data, n_levels=1, min_q=0.0, max_q=0.3, porod_power=4.0, io_guess=None, rg_guess=None):
    self.q = data.q
    self.i = data.i
    self.s = data.s

    self.q_ori = self.q.deep_copy()
    sel = self.q < max_q+0.5 
    self.q_ori = self.q_ori.select( sel )

    sel = self.q < max_q
    self.q = self.q.select( sel )
    self.i = self.i.select( sel )
    self.s = self.s.select( sel ) 
    
    sel = self.q > min_q
    self.q = self.q.select( sel )
    self.i = self.i.select( sel )
    self.s = self.s.select( sel )

    self.n_levels = n_levels
    self.max_q = max_q

    self.i_scale = io_guess/10.0
    self.i = self.i/self.i_scale
    self.s = self.s/self.i_scale 


    self.eps = 1e-8
    if flex.min( self.q )>self.eps:
      self.eps =1e-8 

    self.cnst = 1.0/math.sqrt(6.0)

    self.porod_power = porod_power
    self.rg_guess = rg_guess
    self.io_guess = io_guess

    self.x = None
    self.n = 1+5*self.n_levels
    self.domain = [ ]
    self.domain = [ (0,50), (5.0, 15.0), (rg_guess*0.5, rg_guess*1.5) ]
    for ii in xrange( 3, self.n) :
      self.domain.append( (0,10) )
    self.optimizer =  de.differential_evolution_optimizer(self,
                                                          monitor_cycle=300,
                                                          max_iter=1000000,
                                                          population_size=3*self.n,
                                                          n_cross=3,
                                                          show_progress_nth_cycle=10,
                                                          show_progress=False,
                                                          f=0.95,eps=1e-8)
    self.unifit_curve = self.curve( self.x, self.q_ori )
    #for qq,jj in zip(self.q_ori,  a):
    #  print qq,jj*self.i_scale


  def single_level_curve(self, gi, rgi, bi, rgcfi, pi, q=None):
    gi=abs(gi)
    rgi=abs(rgi)
    bi=abs(bi)
    rgcfi=abs(rgcfi)
    if q is None:
      q = self.q
    result = abs(gi)*flex.exp( -q*q*rgi*rgi/3.0 ) + bi*flex.exp(-q*q*rgcfi*rgcfi/3.0)*flex.pow(
                     flex.pow( sm.erf(q*rgi*self.cnst), 3.0 )/(q+self.eps), pi)
    return result 


  def get_param_vectors(self, vector):
    fb = abs(vector[0])
    single_vectors = []
    for ii in range( self.n_levels ):
      tmp_vector=list( vector[ii*5+1:ii*5+1+5] )
      single_vectors.append( tmp_vector )
    return fb, single_vectors

  def curve(self, vector, q=None):
    if q is None:
      q = self.q
    fb,vecs = self.get_param_vectors(vector)
    result = q*0
    for vec in vecs:
      if self.porod_power is not None:
        pp =  self.porod_power
      else:
        pp = vec[4]
      if self.fix_rg_io:
        vec[0]=self.io_guess
        vec[1]=self.rg_guess

      result += self.single_level_curve( vec[0], vec[1], vec[2], vec[3], pp, q )
    return result+fb
     
  def target(self,vector):
    fb,vecs = self.get_param_vectors(vector)
    tmp = self.curve(vector)
    tmp =tmp-self.i
    tmp = tmp/(self.s+self.eps)
    tmp = tmp*tmp
    tmp = flex.sum( tmp )
    return tmp/self.q.size()



  def print_status(self, min_s, mean_s, sol, txt):
    print "CYCLE:: ", txt
    print "MIN_S, MEAN_S:", min_s,mean_s
    print list(sol)
    io = sol[1]
    b1 = sol[3]
    print abs(io/b1), abs(b1/io) 
    #for qq, ii,cc in zip(self.q, self.i, ic):
    #  print qq,ii,cc
    print "-----------------------------------------"
 

  def show(self, out=None):
    if out == None:
      out = sys.stdout
    fb,vecs = self.get_param_vectors(self.x)
    print >> out, " Unified fit results, with max_q %5.4e   ( %3.2f * external Rg estimate)"%(self.max_q, self.max_q*self.rg_guess)
    print >> out, " Io   : %5.3e"%(abs(vecs[0][0])*self.i_scale)
    print >> out, " Rg   : %5.3e"%(abs(vecs[0][1]))
    print >> out, " B    : %5.3e"%(abs(vecs[0][2])*self.i_scale)
    print >> out, " Rgcf : %5.3e"%(abs(vecs[0][3]))
    io = abs(vecs[0][0])
    b1 = abs(vecs[0][2])
    print >> out, " B+I0 : %5.3e"%( b1+max(io,1e-12) )
    print >> out


 





