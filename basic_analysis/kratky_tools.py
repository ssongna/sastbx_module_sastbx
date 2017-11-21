from scitbx.array_family import flex
from sastbx.data_reduction import curves, saxs_read_write
from scitbx.math import chebyshev_polynome
from scitbx.math import scale_curves
import sys,os,math

class std_kratky(object):
  def __init__(self):
    self.lims = (0,3.999)
    self.n = 25
    self.cheb_coefs_mean = flex.double( [-5.7657057160113165, -1.4401697329647773, 0.27024527974695978, 0.0078985123178169463, -0.09631162819805221, 0.10735910911235932, -0.2233791415486541, 0.19498973606201189, -0.26350484909274735, 0.27066922616706718, -0.25090257310775238, 0.22483113451869022, -0.16726514028401993, 0.1027196386216821, -0.091673985006882924, 0.052259467525186447, -0.016970937165037378, 0.031137149395069535, 0.007431391017472774, 0.014012399241218143, 0.010042763559628737, -0.015341614910420142, 0.0068872073662354753, 0.012869250293871155, 0.008959831186295333] )
    self.cheb_coefs_low=flex.double( [-9.2669907863291705, -1.9456068148305121, 0.86779111670403164, 0.14288741649825987, 0.052348627094772476, -0.075972113378021927, -0.24684531148041183, 0.26186454415838456, -0.11901905265242005, 0.51976185959989885, -0.45937166668405666, 0.31225783754328595, -0.16416500741911164, 0.32535788332409726, -0.17264130742234479, 0.049992516172585558, -0.11856771807344196, 0.044221277456559308, 0.18118326484913538, 0.0550906704796576, 0.047995789360044186, -0.083644664646019343, 0.029230535847010842, -0.040766646640578184, -0.070401435043602661] )
    self.cheb_coefs_high=flex.double( [-2.4313071167283367, -1.0123006175290996, -0.15088667146206206, -0.13439270400442058, 0.0096372912886536471, 0.063670920315285234, -0.21889541306029214, -0.0062901907374255505, -0.17570764254872775, 0.0052730332601467964, 0.010615988987238418, 0.14617014871482245, 0.021082968158927021, -0.0046556220561053659, -0.02634693069601177, 0.026402323365458613, -0.051317703298683792, 0.00089074796092045432, 0.011927144561737445, 0.040880441622645626, 0.028415795670163915, -0.01415133080863344, -0.02323296702274262, -5.8561934669618548e-05, 0.051520439604993466] )

    self.loc = 0.10
    self.scale = 1.0
    self.m   = chebyshev_polynome(self.n, self.lims[0]*0.99999, self.lims[1]/0.99999, self.cheb_coefs_mean )
    self.low = chebyshev_polynome(self.n, self.lims[0]*0.99999, self.lims[1]/0.99999, self.cheb_coefs_low  )
    self.high= chebyshev_polynome(self.n, self.lims[0]*0.99999, self.lims[1]/0.99999, self.cheb_coefs_high )

  def std_kratky_plot(self, q, xscale, yscale, margin=0.95):
    x   = q*self.loc/xscale
    k_m = flex.exp( self.m.f(x)    ) * yscale / self.scale
    k_l = flex.exp( self.low.f(x)  ) * yscale / self.scale*margin
    k_h = flex.exp( self.high.f(x) ) * yscale / self.scale/margin
    return k_m, k_l, k_h

class kratky_scaler(object):
  def __init__(self, max_w = 4.0, n=4000,loc=0.1, dob=True):
    self.max_w    = max_w
    self.interpol = scale_curves.curve_interpolator( 0, max_w, n )
    self.w        = self.interpol.target_x
    self.w_step   = self.interpol.delta
    self.n        = n
    self.loc      = loc
    self.std_kp   = std_kratky()

  def kratky_plots(self,data):
    q_in  = data.q
    i_in  = data.i
    k     = i_in*q_in*q_in
    k     = self.smooth_curve(k,2)
    loc,scale = self.normalize(q_in,k)
    k     = i_in*q_in*q_in
    mk,lk,hk =  self.std_kp.std_kratky_plot(q_in, loc, scale)
    return q_in,k,mk,lk,hk,loc

  def smooth_curve(self,curve,n=5):
    new_curve = curve*0
    for ii in range( curve.size()-n ):
      tmp = flex.mean(curve[ii:ii+n])
      new_curve[ii] = tmp
    for ii in range(int(curve.size()-n),curve.size()):
      new_curve[ii] = new_curve[ii-1]
    return new_curve

  def normalize(self, q, curve):
    last_height = flex.mean( curve[ curve.size()-10:curve.size()-1 ] )
    funct = q/flex.max(q)*last_height
    d = curve - funct
    max_index = flex.max_index( d )
    loc = q[max_index]
    scale = curve[max_index]
    return loc, scale

  def fraction_of_high_violations(self,k,h):
     frac = (k > h).count(True) /(1.0*k.size())
     return frac

  def fraction_above_mean_after_peak(self,peak,q,k,mk,peak_factor=1.25,kratky_factor=1.0):
    ipeak=0
    for ii in range(k.size()):
      if peak*peak_factor > q[ii]:
        ipeak=ii
        break


    dd = k-mk*kratky_factor
    dd = dd[ipeak:]
    n = dd.size()*1.0
    below = (dd<0).count(False)
    return below/n

  def analyze_and_report(self,data,out=None,critical_frac=0.1):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "========  Kratky Analyses ========"
    print >> out

    q,k,mk,lk,hk,loc = self.kratky_plots(data)
    k_fact = 2.00
    p_fact = 2.00
    k_frac = self.fraction_above_mean_after_peak(loc,q,k,mk,p_fact,k_fact)
    print >> out, "The first maximum is found at (roughly) q=%5.2e"%(loc)
    print >> out, "The fraction of points that is systematically higher then "
    print >> out, "%4.2f the standard Kratky plot, is %5.2f for q>%5.2e"%(k_fact,k_frac,p_fact*loc)
    print >> out
    if k_frac>0.15:
      print >> out, "  This could incidate that the background subtraction was non-optimal "
      print >> out, "  or that the sample is partially unfolded."
      print >> out, "  Extremely serious sample issues will show up in the Kratky Outlier test below."

    frac = self.fraction_of_high_violations(k,hk)
    print >> out
    print >> out, "The fraction of Kratky outliers (outside th 99% confidenice interval) is: "
    print >> out, "%5.2f %% "%(frac*100.0)
    print >> out
    if frac > critical_frac:
      print >> out, "This could indicate that the protein is partially unfolded"
    else:
      print >> out




def kratky_analyses(data):
  ko = kratky_scaler()
  ko.analyze_and_report(data)


def run(file_name):
  ks = kratky_scaler()
  data = saxs_read_write.read_standard_ascii_qis(file_name)
  kratky_analyses(data)


if __name__ == "__main__":
  run(sys.argv[1])
