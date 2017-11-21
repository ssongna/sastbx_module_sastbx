from scitbx.array_family import flex
import scitbx.math as sm
from scitbx import differential_evolution as de
import math,sys,os


class porod_engine(object):
  def __init__(self, data, min_q, delta=0.15):
    self.data = data
    self.min_q = min_q
    self.max_q = self.min_q+delta

    self.q = data.q
    self.i = data.i
    self.s = data.s

    sel = (self.q > self.min_q) & (self.q < self.max_q)
    self.q = self.q.select( sel )
    self.i = self.i.select( sel )
    self.s = self.s.select( sel )

    self.n = 2
    self.x = flex.double()

    self.nn = self.q.size()

    self.domain = [ (0,1.0e12),(0,5e3) ]

    self.optimizer =  de.differential_evolution_optimizer(self,
                                                          monitor_cycle=300,
                                                          max_iter=1000000,
                                                          population_size=50,
                                                          n_cross=3,
                                                          show_progress_nth_cycle=1000,
                                                          show_progress=True,
                                                          f=0.75,eps=1e-8)

  def target(self, vector):
    a = abs( vector[0] )
    b = abs( vector[1] )
    tmp = a/(self.q*self.q*self.q*self.q)+b
    tmp = tmp-self.i
    tmp = tmp/self.s
    tmp = tmp*tmp
    tmp = flex.sum(tmp)
    return tmp/self.nn



  def print_status(self, min_s, mean_s, sol, txt):
    print "CYCLE:: ", txt
    print "MIN_S, MEAN_S:", min_s,mean_s
    print list(sol)
    io = abs(sol[0])
    b1 = abs(sol[1])
    print io, b1, self.min_q, self.max_q,
    if txt == "Final":
      print  min_s, "   GREP"
    else:
      print
    print "-----------------------------------------"
 



