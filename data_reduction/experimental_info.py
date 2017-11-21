from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import pdb, ranges, detectors
from iotbx.option_parser import option_parser
import libtbx.phil.command_line
from cStringIO import StringIO
from libtbx.utils import null_out
import sys, os


master_params = iotbx.phil.parse("""\
experiment{
  date=None
  .type=str
  .help="Date of the experiment"
  beamline=None
  .type=None
  .help="End station on which data was collected"
  id=None
  .type=str
  .help="Experiment id"
  concise_output = True
  .type=bool
  .help="When False, individual curves for each image are written out in working directory"
  data_set
  .multiple=True
  .help="Data set defintion"
  {
    specimen
    .multiple=True
    .help="Specimen definition"
    {
      type = sample buffer window
      .type=choice
      .help="Determines if this sample is from sample or background"
      id=None
      .type=str
      .help="Sample identifier"
      details{
        path_length = 1
        .type=float
        .help="Sample thickness"
        concentration = None
        .type = float
        .help="Concentration in mg/ml"
      }
      data
      .multiple=True
      .help="Data files definition"
      {
        base = None
        .help="Image base name"
        .type=str
        range=None
        .type=str
        .help="Range of images"
        file_name_settings
        {
          wildcard = "#"
          .help="Serial id wildcard"
          .type=str
          has_leading_zeros=True
          .type=bool
          .help="Does the serial id have leading zeros?"
        }
      }
      output
      .help="Defines output file names"
      {
        file_name=None
        .help="Output file name for integrated data from this specimen. If None, filename will be constructued from input parameters."
      }
   }
    hutch_setup{
      wavelength=None
      .type=float
      .help="Wavelength"
      detector{
        type=ADSC *EDF
        .help="Detector type"
        .type=choice
        nx=None
        .type=int
        .help="Number of pixels along X"
        ny=None
        .type=int
        .help="Number of pixels along Y"
        scale_x = None
        .type=float
        .help="Pixel scale along X"
        scale_y = None
        .type = float
        .help = "Pixel scale along Y"
        distance=None
        .type=float
        .help="Sample detector distance"
        beam_origin_x=None
        .type=float
        .help="Direct beam coordinates X (in pixels from detector origin)"
        beam_origin_y=None
        .type=float
        .help="Direct beam coordinates Y (in pixels from detector origin)"
        angle_alpha = 0
        .type=float
        .help="Detector misset angle alpha"
        angle_beta = 0
        .type=float
        .help="Detector misset angle beta"
        angle_gamma = 0
        .type=float
        .help="Detector misset angle gamma"
        aux_info{
          location = *image_header additional_file
          .type = choice
          .help = "Where to find much needed information?"
          file_name = None
          .type=path
          .help="The file with aux info if not available in header"
          file_type = *image_name_i0_i1
          .type=choice
          .help="Format of file"
          prefix=''
          .type=str
          .help="prefix to add to file name in aux file to make things work"
        }
        exclude
        .help="Exclusion areas on the detector"
        {
          circle
          .help="A circular exlusion area"
          {
            origin = None
            .type= str
            .help="coordinates of exclusion circle"
            radius=None
            .type=float
            .help="Radius of exclusion circle"
          }
          polygon
          .multiple=True
          .help="Polygon exclusion area"
          {
            corner = None
            .type = str
            .help="Corner of exclusion zone. We need at least 3 of these coordinates."
            .multiple = True
          }
          from_image
          .help="Auto detect dead pixels from an input image"
          {
            file_name=None
            .type=path
            .help="file name of image for auto detection of dead pixel mask"
          }
        }
        shadow
        .help="Shadow area for offset and approximate PSF corrections"
        {
          circle
          .help="A circular shadow area"
          {
            origin = None
            .type= str
            .help="Coordinates of shadow circle"
            radius=None
            .type=float
            .help="Radius of shadow circle"
          }
          polygon
          .multiple=True
          .help="Polygon shadow area"
          {
            corner = None
            .type = str
            .help="Corner of shadow zone. We need at least 3 of these coordinates."
            .multiple = True
          }
        }
        header_keys{
            exposure_time = 'TIME'
            .type = str
            .help = 'Header keyword describing exposure time.'

            I0 = I0_1,I0_2
            .type = str
            .help = 'Header keyword(s) describing intensity before sample. Used for transmission correction.'
            I1 = I1_1,I1_2
            .type = str
            .help = 'Header keyword(s) describing intensity after sample. Used for transmission correction.'

            distance = 'DISTANCE'
            .type=str
            .help='Header keyword describing detector to sample distance.'

            wavelength = 'WAVELENGTH'
            .type = str
            .help = 'Header keyword describing wavelength in Angstrom'

            beam_origin_x = 'BEAM_CENTER_X'
            .type = str
            .help = 'Header keyword describing beam center x'
            beam_origin_y = 'BEAM_CENTER_Y'
            .type = str
            .help = 'Header keyword describing beam center y'

            overload = 'CCD_IMAGE_SATURATION'
            .type = str
            .help = 'Header keyword describing overload value of pixels'

            nx = 'SIZE1'
            .type=str
            .help = 'Header keyword describing Number of pixels along x'
            ny = 'SIZE2'
            .type=str
            .help='Header keyword describing Number of pixels along y'

            scale_x = 'PIXEL_SIZE'
            .type = str
            .help = 'Header keyword describing pixel scale along x'
            scale_y = 'PIXEL_SIZE'
            .type = str
            .help = 'Header keyword describing pixel scale along y'


          }

      }
    }
 }
}
processing_options{
  integration{
    q_range{
      min_q = 0
      .type=float
      .help="Lowest q in data procesing."
      max_q = None
      .type=float
      .help="Maximum q in data processing. Leave unspecified for full detector."
      q_step=0.0005
      .type=float
      .help="step size in binning"
    }
  }
  scaling{
    averaging{
      scaling{
        scale =  *from_header linear None
        .type=choice
      }
    }
  }
}


""")


class polygon(object):
  def __init__(self,params):
    self.corners_x = []
    self.corners_y = []
    for xy in params.corner:
      x,y = self.parse( xy )
      self.corners_x.append( x )
      self.corners_y.append( y )

  def parse(self, txt):
    keys = txt.split(" ")
    result = []
    for key in keys:
      tmp=None
      try:
        tmp = int( key )
      except: pass
      if tmp is not None:
        result.append( tmp )
    if len(result)==2:
      return result[0], result[1]
    else:
      return None,None


  def show(self,out=None):
    if out is None:
      out = sys.stdout
    for x,y in zip( self.corners_x, self.corners_y ):
      print >> out, "   ---->   ", x, y

class circle(object):
  def __init__(self,params):
    self.x, self.y = self.parse( params.origin)
    self.radius = None
    if  params.radius  is not None:
      self.radius  = int( params.radius )

  def parse(self, txt):
    if txt is not None:
      keys = txt.split(" ")
      result = []
      for key in keys:
        tmp=None
        try:
          tmp = int( key )
        except: pass
        if tmp is not None:
          result.append( tmp )
      if len(result)==2:
        return result[0], result[1]
      else:
        return None,None
    else:
      return None,None

class specimen_specs(object):
  def __init__(self, params):
    self.id = params.id
    self.type = params.type
    self.path_length = params.details.path_length
    self.concentration = params.details.concentration
    self.file_names = []
    for data_set in params.data:
      tmp = self.get_file_names(data_set.base, data_set.file_name_settings.wildcard,  data_set.file_name_settings.has_leading_zeros, data_set.range)
      for file in tmp:
        self.file_names.append( file )
    self.output = params.output.file_name
    if self.output is None:
      self.output = self.id+".qis"
    if self.output is None:
      self.output = data_set.base+".qis"

  def get_file_names(self, base_name, wildcard, has_leading_zeros, user_range_txt):
    file_name_engine = ranges.serial_file_name_handler(base_name, wildcard, has_leading_zeros)
    id_list = ranges.range_to_list( ranges.range_parser(user_range_txt) )
    result = file_name_engine.names_from_range_list( id_list )
    return result


  def show(self, out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "Specimen ID            :  %s "%str(self.id)
    print >> out, "Specimen type          :  %s "%str(self.type)
    print >> out, "Specimen concentration :  %s "%str(self.concentration)
    print >> out
    for file in self.file_names:
      print >> out, "             ", file

class hutch_specs(object):
  def  __init__(self, params):
    self.params = params

    self.detector_type = self.params.detector.type

    self.nx = None
    self.ny = None

    self.distance = None
    self.wavelength = None

    self.beam_origin_x = None
    self.beam_origin_y = None
    self.scale_x = None
    self.scale_y = None

    self.angle_alpha = 0
    self.angle_beta = 0
    self.angle_gamma = 0

    self.exclude_polygons = []
    self.exclude_circle = None
    self.exclude_from_image = None

    self.shadow_polygons = []
    self.shadow_circle = None

    self.aux_info_object = None



  def set_values(self,sample_image):
    # read the header of a sample image
    img_factory = detectors.TrySingleImageType(sample_image,self.detector_type)
    img_factory.readHeader()
    header = img_factory.parameters
    # fill in values we need on the basis of a sample image
    self.wavelength = self.set_keys( header, self.params.detector.header_keys.wavelength, self.params.wavelength )
    self.distance = self.set_keys( header, self.params.detector.header_keys.distance, self.params.detector.distance )
    self.beam_origin_x = self.set_keys( header, self.params.detector.header_keys.beam_origin_x, self.params.detector.beam_origin_x )
    self.beam_origin_y = self.set_keys( header, self.params.detector.header_keys.beam_origin_y, self.params.detector.beam_origin_y )
    self.nx = self.set_keys( header, self.params.detector.header_keys.nx, self.params.detector.nx )
    self.ny = self.set_keys( header, self.params.detector.header_keys.ny, self.params.detector.ny )
    self.scale_x = self.set_keys( header, self.params.detector.header_keys.scale_x, self.params.detector.scale_x)
    self.scale_y = self.set_keys( header, self.params.detector.header_keys.scale_y, self.params.detector.scale_y)

    if self.params.detector.aux_info.location == "additional_file":
      from aux_info_external_file_parsing import aux_info_simple
      self.aux_info_object = aux_info_simple( self.params.detector.aux_info.file_name, prefix=self.params.detector.aux_info.prefix )
         

    # exclusion
    for poly in self.params.detector.exclude.polygon:
      if len(poly.corner)>2:
        tmp = polygon( poly )
        self.exclude_polygons.append( tmp )
    if self.params.detector.exclude.circle.radius is not None:
      self.exclude_circle = circle( self.params.detector.exclude.circle )

    self.exclude_from_image = self.params.detector.exclude.from_image.file_name

    # shadow
    for poly in self.params.detector.shadow.polygon:
      tmp = polygon(poly)
      self.shadow_polygons.append( tmp )
    if self.params.detector.shadow.circle.radius is not None:
      self.shadow_circle = circle( self.params.detector.shadow.circle )

  


  def set_keys(self, header, header_key, user_input):
    target=None
    if user_input is None:
      try:
        target = header[ header_key ]
      except: pass
    else:
      target = user_input
    if target is None:
      raise Sorry("Failed to set essential key %s. Check input or image header"%(header_key) )

    return target

  def construct_polygon_parameters(self, polygon_parameters):
    corners_x = []
    corners_y = []
    for xy in polygon_parameters.corner:
      print xy

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "---- Hutch setup  ----"
    print >> out
    print >> out, "   wavelength    : ", self.wavelength
    print >> out, "   distance      : ", self.distance
    print >> out, "   xbeam, ybeam  : ", self.beam_origin_x , self.beam_origin_y
    print >> out, "   nx,ny         : ", self.nx, self.ny
    print >> out
