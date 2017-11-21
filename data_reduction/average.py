from scitbx.array_family import flex
from scitbx.math import scale_curves
import experimental_info, curves, saxs_read_write
import sys, os, math

class averaging_engine(object):
  def __init__(self, specimen, params):
    self.specimen = specimen
    self.params = params
    # these variables will be constructed from the specimen given
    self.q_range = None
    self.mean_i = None
    self.var_i = None

    self.report = """"""

    self.mean_trans = None
    self.std_trans = None
    self.set_transmissions()

  def set_transmissions(self):
    m = 0
    v = 0 
    n = 0
    for intdat in self.specimen.integrated_data:
      print intdat.transmission_factor
      m += intdat.transmission_factor
      v += intdat.transmission_factor*intdat.transmission_factor
      n += 1.0
    m = m/n
    v = v/n-m*m
    self.mean_trans = m 
    self.std_trans = math.sqrt(v)

  def scale_to_reference(self, ref_id,chuck_size=50):
    #here we scale the data
    total_scales = []
    header_scales = []
    refined_scales = []
    n_data_sets = len(self.specimen.integrated_data)
        
    if (self.params.scaling.scale == "from_header") or (self.params.scaling.scale == "linear"):
      header_scales = self.normalize_by_i1()
      self.build_report(scales=header_scales,txt="    Scales I0 from image header     ")
      
    if self.params.scaling.scale == "linear":
      m = []
      v = []
      for ii in xrange(n_data_sets): 
        n = self.specimen.integrated_data[ii].mean.size()
        while n < 3*chuck_size:
          chuck_size = chuck_size-1
        tmp_m = self.specimen.integrated_data[ii].mean[chuck_size:n-chuck_size]
        tmp_v = self.specimen.integrated_data[ii].var[chuck_size:n-chuck_size]
        m.append(tmp_m)
        v.append(tmp_v)
      refined_scales, refined_offsets = scale_curves.scale_it_pairwise( m,v,ref_id,factor=20,show_progress=False,add=False)
      self.scale_with_supplied_scale(refined_scales,refined_offsets)
      self.build_report(refined_scales,refined_offsets,"    Refined scales    ")

    if self.params.scaling.scale == "None":
      total_scales  = [1.0]*n_data_sets
      total_offsets = total_offsets[0.0]*n_data_sets
      self.scale_with_supplied_scale(total_scales, ref_id=ref_id)
      self.build_report(total_scales,total_offsets,"    No scaling   ")

 
  def normalize_by_i0(self):
    scales = []
    for ii in xrange( len(self.specimen.integrated_data) ):
      scale_factor = 1.0/self.specimen.integrated_data[ii].i_tot
      self.specimen.integrated_data[ii].multiply( scale_factor )
      scales.append( scale_factor )
    return scales

  def normalize_by_i1(self):
    scales = []
    for ii in xrange( len(self.specimen.integrated_data) ):
      scale_factor = 1.0/self.specimen.integrated_data[ii].i_after
      self.specimen.integrated_data[ii].multiply( scale_factor )
      scales.append( scale_factor )
    return scales

  def scale_with_known_scale(self):
    scales = []
    for ii in xrange( len(self.specimen.integrated_data) ):
      scale_factor = 1.0 / self.specimen.integrated_data[ii].i_after
      self.specimen.integrated_data[ii].multiply( scale_factor )
      scales.append( scale_factor )
    return scales

  def scale_with_supplied_scale(self,scales,offsets):
    for ii in xrange( len(self.specimen.integrated_data) ):
      scale = scales[ii]
      offset = offsets[ii]
      self.specimen.integrated_data[ii].multiply( scale )
      self.specimen.integrated_data[ii].add( offset )

  def get_average_data(self):
    n = len(self.specimen.integrated_data)
    q = self.specimen.integrated_data[0].q_range.q_mean
    i_tot, v_tot = self.specimen.integrated_data[0].mean_var()
    for ii in xrange( 1,n ):
       i,v = self.specimen.integrated_data[ii].mean_var()
       print ii, flex.mean(i)
       i_tot += i
       v_tot += v
    i_tot = i_tot/n
    v_tot = v_tot/n
    return q, i_tot, flex.sqrt(v_tot)

  def build_report(self,scales,offsets=None,txt=""):
    self.report+="---- Results of scaling data ----\n"
    self.report+=txt
    self.report+="  \n"
    if offsets is None:
      offsets = [0]*len(scales)
    for image_name,scale,offset,id in zip( self.specimen.image_names(), scales,offsets, xrange(len(scales)) ):
      self.report+="  id: %5i    image: %s  -- scale, offset -->  %8.3e     %8.3e"%(id,image_name,scale,offset)
      self.report+="\n"
    self.report += "\n\n"
    table = """R value analyses\n"""
    table+= "\n"
    table+= "R = 100%*sum( abs(I_i-I_j) )/sum( abs(I_i+I_j)/2 )"
    table+= "\n"
    table+= "\n"
    table += "       "

    for ii in xrange(len(self.specimen.integrated_data)):
      table += "%6i"%ii
    table +="\n"

    for ii in xrange(len(self.specimen.integrated_data)):
      table += "%6i  "%ii
      for jj in range(0,ii):
        table += "      "
      for jj in xrange( ii,len(self.specimen.integrated_data) ):
        rval = self.r_value(ii,jj)
        table += " %3.2f "%(rval)
      table +="\n"
    self.report += table
    self.report +="\n\n\n\n"



  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out, self.report


  def r_value(self,ii,jj):
    i_ii,s_ii = self.specimen.integrated_data[ii].mean_var()
    i_jj,s_jj = self.specimen.integrated_data[jj].mean_var()
    top = flex.sum( flex.abs(i_ii-i_jj) )
    bottom = flex.sum( flex.abs(i_ii+i_jj)/2.0 )
    return 100.0*top/max(bottom,1e-13)

  def get_average_simple_saxs_data(self):
    q,i,s = self.get_average_data()
    ssd = curves.simple_saxs_data( q,i,s, 
                                   concentration_mg_ml=self.specimen.params.concentration,
                                   data_type=self.specimen.params.type, 
                                   data_id=self.specimen.params.id)

    return ssd

  def write_average_intensity(self):
    ssd = self.get_average_simple_saxs_data()
    file_name = self.specimen.params.output
    saxs_read_write.write_standard_ascii_qis( ssd, file_name )

"""
class kapton_scaler(object):
  def __init__(self, sample, background, start_q=0.33, stop_q=0.43):
    #  Assumption: at high q-values (arouind the kapton peak to be specific),
    #  we assume that 
    #    background*scale + offset = sample
    #  minimizing both scale and offset in a least squares fashion, will provide us with additional scale and offset 
    #  
    self.sample = sample
    self.background = background
    self.start_q = start_q
    self.stop_q = stop_q
    self.tmp_q = self.sample
"""


"""

We have two main scaling methods:

- on the basis of ion chamber readings

- a least squares protocol

** Method 1
    I_a(Q)' = I_a(Q)*(I_tot_ref/I_tot_a)
    with I_tot_x being the total received intensity on sample x

** Method 2
    A least squares method: this method uses a pairwise scaling (fast) to a reference dataset

"""
