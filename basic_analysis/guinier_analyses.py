import scitbx.math
from scitbx.array_family import flex
import math,os,sys
from scitbx import differential_evolution as de
from scitbx.math import gamma_incomplete, gamma_incomplete_complement
from scitbx.math import gamma_complete


class guinier_accumulator(object):
  def __init__(self):
    self.rg = flex.double()
    self.i0 = flex.double()
    self.scores = flex.double()

  def add_data(self, rg, i0, free_score):
    self.rg.append( rg )
    self.i0.append( i0 )
    self.scores.append( free_score  )

  def weighted_means(self):
    min_score = flex.min( self.scores )
    mw = flex.mean( self.scores )
    self.scores = self.scores-min_score
    wghts = flex.exp( -self.scores*0.50 )
    sw = 1e-12+flex.sum( wghts )
    mrg   = flex.sum(wghts*self.rg)/sw
    srg   = flex.sum(wghts*self.rg*self.rg)/sw
    mi0   = flex.sum(wghts*self.i0)/sw
    si0   = flex.sum(wghts*self.i0*self.i0)/sw
    si0   = math.sqrt(si0-mi0*mi0)
    srg   = math.sqrt(srg-mrg*mrg)
    return mrg,srg,mi0,si0,mw

  def show(self, min_q, max_q, out=None):
    if out is None:
      out = sys.stdout
    mrg, srg,mi0,si0,mw = self.weighted_means()
    print >> out, "    Guinier analyses "
    print >> out, "       min_q :  %6.4e      "%(min_q)
    print >> out
    print >> out, "    Results of Guinier analyses of small q-range chunks  "
    print >> out
    print >> out, "       Rg    :  %6.4e      sigma :  %6.4e     (ratio: %4.3e)"%(mrg,srg,mrg/(srg+1e-12))
    print >> out, "       I0    :  %6.4e      sigma :  %6.4e     (ratio: %4.3e)"%(mi0,si0,mi0/(si0+1e-12))
    print >> out, "       mean chi_squared :  %5.4e "%(mw)
    print >> out



class guinier_engine(object):
  def __init__(self, data, min_q=None, max_q=None):
    self.data = data
    self.min_q = min_q
    if self.min_q == None:
      self.min_q = self.data.q[0]
    self.max_q = max_q
    self.lims = self.create_limits()
    self.accumulator = guinier_accumulator()
    self.all_data = self.compute_all()
    self.w_rg, self.w_srg, self.w_mi0, self.w_si0, self.mw = self.accumulator.weighted_means()
    self.bad_sections = self.q_range_analyses( self.w_rg, self.w_mi0)
    self.rg, self.io, self.target, self.n = self.estimate_rg_on_basis_of_q_range_analyses( self.bad_sections, self.w_rg, self.w_mi0)

  def show(self, out=None):
    if out is None:
      out = sys.stdout
    mrg = self.w_rg
    srg = self.w_srg
    mi0 = self.w_mi0
    si0 = self.w_si0
    mw  = self.mw
    print >> out
    print >> out
    print >> out
    print >> out
    print >> out, "-----------------------------------------------------------------"
    print >> out, "    Guinier analyses "
    print >> out, "       min_q :  %6.4e      "%(self.min_q)
    print >> out
    print >> out, "    Results of Guinier analyses of small q-range chunks  "
    print >> out
    print >> out, "       Rg    :  %6.4e      sigma :  %6.4e     (ratio: %4.3e)"%(mrg,srg,mrg/(srg+1e-12))
    print >> out, "       I0    :  %6.4e      sigma :  %6.4e     (ratio: %4.3e)"%(mi0,si0,mi0/(si0+1e-12))
    print >> out, "       mean chi_squared :  %5.4e "%(mw)
    print >> out
    print >> out, "    'Bad' regions "
    for sec in self.bad_sections:
      print >> out, "           stop_q: %5.4e   start_q: %5.4e "%(sec[0], sec[1])
    print >> out
    print >> out, "    Using data between %6.4e and %6.4e and excluding bad regions"%(self.min_q,self.max_q)
    print >> out
    print >> out, "       Rg   :  %6.4e "%( self.rg  )
    print >> out, "       I0   :  %6.4e "%( self.io  )
    print >> out, "-----------------------------------------------------------------"
    print >> out
    print >> out
    print >> out
    print >> out


  def compute_all( self ):
    self.rg2s = []
    self.lnis = []
    self.free_scores = []
    self.stop_qs = []
    self.start_qs = []
    for lim in self.lims:
      rg2, lni, score, free_score, stop_q, start_q = self.compute_rg( lim[0], lim[1] )
      rg2, lni, score, free_score, start_q, stop_q = self.filter( rg2, lni, score, free_score, stop_q, start_q  )
      if rg2 is not None:
        self.accumulator.add_data( math.sqrt(rg2), math.exp(lni), free_score )
        self.rg2s.append( math.sqrt(rg2) )
        self.lnis.append( lni )
        self.free_scores.append( score )
        self.stop_qs.append( stop_q )
        self.start_qs.append( stop_q )

  def filter(self, rg2, lni, score, free_score, stop_q, start_q, rat_lim=1.3 ):
    if rg2 > 0:
        if (math.sqrt( rg2)*stop_q)> rat_lim:
          if free_score is not None:
            return rg2, lni, score, free_score, start_q, stop_q
    return None,None,None,None,None,None

  def chi_square( self, lni, rg2, q,i,s, mean=True):
    i_calc = flex.exp(lni-rg2*q*q/3.0)
    d = (i-i_calc)/s
    if mean:
      d = flex.mean( d*d )
    else:
      d = d*d
    return d

  def free_score(self, lni, rg2, start_q, stop_q):
    selection_low  = flex.bool(self.data.q<start_q)
    selection_high = flex.bool(self.data.q>stop_q)
    selection_rg   = flex.bool(self.data.q*math.sqrt(abs(rg2))<1.3)
    selection_above_min_q = flex.bool(self.data.q >= self.min_q )
    #for q,i,j,k in zip(self.data.q,selection_low,selection_high,selection_rg):
    #  print q,i,j,k, start_q, stop_q, math.sqrt(abs(rg2))*q
    tot_sel = ((selection_low | selection_high) & selection_above_min_q ) &  selection_rg
    tmp_q = self.data.q.select( tot_sel )
    tmp_i = self.data.i.select( tot_sel )
    tmp_s = self.data.s.select( tot_sel )
    score = None
    if tmp_s.size()>0:
      score = self.chi_square( lni, rg2, tmp_q, tmp_i, tmp_s )
    return score

  def estimate_rg_on_basis_of_q_range_analyses(self, bad_sections, rg_estimate,io_estimate, rat=1.3):
    master_selector = flex.bool( self.data.q > self.min_q ) # all true
    rg_based_selector = flex.bool(self.data.q*rg_estimate<rat)
    self.max_q = rat/rg_estimate
    master_selector  = master_selector & rg_based_selector
    for section in bad_sections:
      low_toss  = flex.bool( self.data.q >= section[0] )
      high_toss = flex.bool( self.data.q <= section[1] )
      combo_toss = ~(low_toss & high_toss)
      master_selector = master_selector&combo_toss
    tmp_q = self.data.q.select(master_selector)
    tmp_i = self.data.i.select(master_selector)
    tmp_s = self.data.s.select(master_selector)
    de_engine =  rg_fitter(tmp_q,tmp_i,tmp_s,rg_estimate*2,io_estimate*3 )
    rg,io = de_engine.return_solution()
    t = de_engine.return_final_target()
    return rg,io,t,tmp_q.size()

  def q_range_analyses(self, rg, io, rat_lim=1.5, window_size=3, level=10.0, sigma=False):
    selector = flex.bool(self.data.q*rg<rat_lim)
    tmp_q = self.data.q.select( selector )
    tmp_i = self.data.i.select( selector )
    tmp_s = self.data.s.select( selector )
    rg2 = rg*rg
    lni = math.log( io )
    cs = None
    if sigma:
      cs = self.chi_square( lni, rg2, tmp_q, tmp_i, tmp_s,False )
      cs = flex.sqrt( cs )
    else:
      cs = flex.exp(lni-rg2*tmp_q*tmp_q/3.0)
      ss = cs/100.0
      cs = flex.abs(tmp_i-cs)/ss

    not_okai_ranges = []
    previous_one_was_bad=False
    tmp_range = []
    for ii in xrange( window_size, cs.size() ):
      tmp_cs = flex.mean( cs[ii-window_size:ii] )
      if tmp_cs > level:
        if not previous_one_was_bad:
          tmp_range.append( tmp_q[ii] )
          tmp_range.append( tmp_q[ii] )
          previous_one_was_bad=True
        else:
          tmp_range[1] = tmp_q[ii]
      else:
        previous_one_was_bad=False
        if len(tmp_range)>0:
          not_okai_ranges.append( tmp_range )
        tmp_range=[]
    return not_okai_ranges

  def compute_rg_from_data(self,q,i):
    q_sq = q*q
    ln_i = flex.log( i )
    cc_obj = flex.linear_regression( q_sq, ln_i )
    rg2 = -cc_obj.slope()*3.0
    lni = cc_obj.y_intercept()
    return rg2, lni


  def compute_rg(self, start_i, stop_i):
    tmp_q = self.data.q[start_i:stop_i+1]
    tmp_i = self.data.i[start_i:stop_i+1]
    tmp_s = self.data.s[start_i:stop_i+1]

    selection = flex.bool( tmp_i > 0)
    tmp_q = tmp_q.select( selection )
    tmp_i = tmp_i.select( selection )
    tmp_s = tmp_s.select( selection )
    rg2, lni = self.compute_rg_from_data( tmp_q, tmp_i )
    score = self.chi_square( lni, rg2, tmp_q, tmp_i, tmp_s )
    fscore= self.free_score( lni, rg2, tmp_q[0], tmp_q[ len(tmp_q)-1 ] )
    return rg2, lni, score, fscore, tmp_q[0], tmp_q[ len(tmp_q)-1 ]

  def create_limits(self,min_rg=5,points=5):
    low = 0.0
    if self.min_q is not None:
      low= self.min_q
    else:
      low= self.data.q[0]

    high = 1.0/(min_rg*1.3)
    if self.max_q is None:
      self.max_q = high
    else:
      high = self.max_q
    if high > flex.max(  self.data.q ):
      high = flex.max(  self.data.q )


    lims = []
    done=False
    add = 0
    while not done:
      start_q_i = 0+add
      end_q_i = start_q_i+points
      if end_q_i < self.data.q.size():
        if self.data.q[end_q_i] > high:
          done=True
        else:
          lims.append( (start_q_i,end_q_i ) )
          add +=1
      else:
        done=True
    return lims

  def curve_from_fit(self, rg, lni0):
    i_calc = flex.exp( - rg*rg*self.data.q*self.data.q/3.0 +lni0)
    for q,ic,io,s in zip(self.data.q, i_calc, self.data.i, self.data.s):
      print q, ic, io, s


class multi_step_rg_engine(object):
  def __init__(self, data, out=None):
    self.out = out
    if self.out is None:
      self.out = sys.stdout

    self.min_q_trials = range(42)  # [0.0, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02]
    self.data = data
    self.user_q_min = data.q[0]

    self.gao = []
    self.rg_estimates = flex.double()
    self.io_estimates = flex.double()
    self.n_bad = flex.double()
    self.n_chunks = flex.double()
    self.target = flex.double()
    self.n_obs = flex.double()
    self.bad_sections = []
    self.low_cut = flex.double()
    self.min_q = []
    self.median_n = 20
    self.span_n = 5

    print >> out
    print >> out, "Performing Rg/Io analyses at various truncation levels ..."
    print >> out
    
    for min_q in self.min_q_trials:

      tmp_min_q = min_q/1000.0+self.user_q_min
      self.min_q.append( tmp_min_q )
      ge = guinier_engine(self.data,tmp_min_q)
      self.gao.append( ge )
      self.rg_estimates.append( ge.rg )
      self.io_estimates.append( ge.io )
      n_bad, n_chunks, suggested_low_res_cut = self.analyze_trunc_and_bad_sections( tmp_min_q, ge.bad_sections )
      t = ge.target
      n = ge.n
      self.n_bad.append(n_bad)
      self.n_chunks.append( n_chunks )
      self.low_cut.append( suggested_low_res_cut )
      self.target.append( ge.target )
      self.bad_sections.append( ge.bad_sections )
    ge_order = flex.sort_permutation( self.rg_estimates)

    this_one                 = ge_order[self.median_n]
    self.median_rg           = self.rg_estimates[ this_one ]
    self.median_io           = self.io_estimates[ this_one ]
    self.median_min_q        = self.min_q[ this_one ]
    self.median_bad_sections = self.bad_sections[ this_one ]
    self.median_low_res_cut  = self.low_cut[ this_one ]

    suggested_truncation = 0.0
    if len(self.median_bad_sections)>0:
      suggested_truncation= self.median_bad_sections[0][1]
    cons_ind=0
    for q in self.min_q:
      cons_ind+=1
      if suggested_truncation>q:
        break
    #that_one = cons_ind
    #self.cnss_rg           = self.rg_estimates[ that_one ]
    #self.cnss_io           = self.io_estimates[ that_one ]
    #self.cnss_min_q        = self.min_q[ that_one ]
    #self.cnss_bad_sections = self.bad_sections[ that_one ]



    self.build_report()

  def build_report(self,out):
    # here we build a brief report
    result = """"""
    result += " ==== Estimates of Rg and Io ====\n"
    result += "\n"
    result += " Analyses of Rg estimates at various data truncation levels suggests \n"
    result += " Rg : %8.2f\n"%(self.median_rg)
    result += " I0 : %8.2e\n"%(self.median_io)
    if len(self.median_bad_sections) > 0:
      result += " Worrysome sections (posibly due to aggregation / structure / bad data): \n"
      for bs in self.median_bad_sections:
        result += "   %4.3e  -- %4.3e  \n"%(bs[0],bs[1])
      result += "\n"
      result += "\n"
      result += " Low q cut suggested for further automated analyses \n"
      result += "   %6.4e \n"%( self.median_low_res_cut)

    else:
      result += " No evidence for aggregation / structure / bad data found \n"
    #print result
    f = open(out,"a")
    f.write(result)
    f.close()
    return result


  def mod_scores( self, x, k ):
    if x> 10: x=10.0
    gam = gamma_complete(k/2.0)
    part1 = 2**(k/2.0)
    part2 = x**((k/2.0)-1)
    part3 = math.exp(-x/2.0)
    tot1  = part2*part3/(part1*gam)

    bart1 = gamma_incomplete_complement(k/2.0, x/2.0)
    tot2  = bart1/gam

    return math.log(tot1), math.log(tot2)

  def analyze_trunc_and_bad_sections(self, min_q, bad_sections, min_step=3, max_low_res_cut_value=5.0e-2 ):
    step_size = self.data.q[1]-self.data.q[0]
    min_q_fom_data = self.data.q[0]
    tot_delta = 0
    tot_large_chunks = 0
    suggested_low_res_cut_value = self.data.q[0]
    for section in bad_sections:
      delta = section[1]-section[0]
      n_delta = int(delta/step_size)+1
      tot_delta += n_delta
      if n_delta >= min_step:
        tot_large_chunks += 1
        if section[1]<max_low_res_cut_value:
          suggested_low_res_cut_value = section[1]

    return tot_delta, tot_large_chunks, suggested_low_res_cut_value









class rg_fitter(object):
  def __init__(self,q,i,s,up_rg, up_io):
    self.q = q
    self.i = i
    self.s = s

    self.x = None
    self.n = 2
    self.domain = []
    self.domain.append( (0,up_rg) )
    self.domain.append( (0,up_io ) )
    self.optimizer =  de.differential_evolution_optimizer(self,
                                                          monitor_cycle=100,
                                                          max_iter=1000000,
                                                          population_size=20,
                                                          n_cross=3,
                                                          show_progress_nth_cycle=500,
                                                          show_progress=False,
                                                          f=0.95,eps=1e-18)

  def return_solution(self):
    return abs(self.x[0]), abs(self.x[1])

  def return_final_target(self):
    return self.target( self.x )

  def i_calc(self,rg,io,q):
    result = io*flex.exp(-rg*rg*q*q/3.0)
    return result

  def get_params_from_vector(self,vector):
    rgs  = abs( vector[0] )
    ios  = abs( vector[1] )
    return rgs, ios

  def target(self,vector):
    rgs,ios = self.get_params_from_vector(vector)
    i_tot = self.i_calc(rgs,ios,self.q)
    score = flex.sum( flex.pow(((i_tot-self.i) /(self.s+1e-12)),2) )
    return score

  def print_status(self, min_s, mean_s, sol, txt):
    print "CYCLE:: ", txt
    print "MIN_S, MEAN_S:", min_s,mean_s
    rgs,ios = self.get_params_from_vector(sol)
    print rgs, ios
    ic =  self.i_calc(rgs,ios,self.q)
    #for qq, ii,cc in zip(self.q, self.i, ic):
    #  print qq,ii,cc
    print "-----------------------------------------"
