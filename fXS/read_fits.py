import sys, os
from scitbx.array_family import flex
from libtbx.utils import Sorry
has_numpy=False
try:
  import numpy
  has_numpy=True
except: pass
if not has_numpy:
  raise Sorry("Numpy is not installed")



def read_and_return_fits_file(filename):
  hdulist = pyfits.open( filename )
  hdulist.readall()
  hdulist.info()
  data = hdulist['IMAGES'].data
  nxny = data.shape
  grid = flex.grid( nxny )
  ndata = []
  xc=0
  yc=0
  for row in data:
    nrow=[]
    yc=0
    for datum in row:
      ndata.append( float(datum) )
    xc+=1
  flex_data = flex.double( ndata )
  return flex_data, nxny


class polygon(object):
  def __init__(self,zone):
    self.corners_x = []
    self.corners_y = []
    for xy in zone.corners:
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
  def __init__(self,circle):
    self.x, self.y = self.parse(circle.origin)
    self.radius = None
    if  params.radius  is not None:
      self.radius  = int( circle.radius )

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


class image_manager(object):
   def __init__(self, shadow, exclude):
     self.shadow = shadow
     self.exclude = exclude
