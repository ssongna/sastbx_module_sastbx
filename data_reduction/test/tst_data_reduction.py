from sastbx import data_reduction
from iotbx import detectors
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import math
import pickle

def tst_mask_definitions():
  corners = ( ( 5, 5),
              ( 5,10),
              (10,10),
              (10, 5) )
  cx = flex.int([5,5,10,10] )
  cy = flex.int([5,10,10,5])
  N=15
  mask_gen = data_reduction.image_mask(N,N)
  mask_gen.add_polygon(cx,cy,0,0)
  mask = mask_gen.mask()
  count=0
  desired_mask="""0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
"""

  result = """"""
  for ii in xrange(N):
    for jj in xrange(N):
      result+=str(mask[count])+" "
      count +=1
    result=result[0:len(result)-1]
    result+="\n"
  assert desired_mask == result

  """
  y--->
  0 1 2 3 4 5 6 7 8 9
x - - - - - - - - - - 0
| - - - - * - - - - - 1
| - - - # # # - - - - 2
| - - # # # # # - - - 3
V - - # # # # # - - - 4
  - * # # # # # * - - 5
  - - # # # # # - - - 6
  - - - # # # - - - - 7
  - - - - * - - - - - 8
  - - - - - - - - - - 9

  (1,4)
  (5,1)
  (8,4)
  (5,7)

  """



  desired_mask =  """0 0 0 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 1 1 1 0 0 0 0
0 0 0 1 1 1 1 0 0 0
0 0 1 1 1 1 1 0 0 0
0 1 1 1 1 1 1 1 0 0
0 0 1 1 1 1 1 0 0 0
0 0 0 1 1 1 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
"""


  cx = flex.int([1,5,8,5])
  cy = flex.int([4,1,4,7])
  N=10
  mask_gen = data_reduction.image_mask(N,N)
  mask_gen.add_polygon(cx,cy,0,0)
  mask = mask_gen.mask()
  count=0
  result = """"""
  for ii in xrange(N):
    for jj in xrange(N):
      result+=str(mask[count])+" "
      count +=1
    result=result[0:len(result)-1]
    result+="\n"
  assert result == desired_mask


  desired_mask = """0 0 0 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 1 1 0 0 0 0 0
0 0 0 1 1 0 0 0 0 0
0 0 1 1 1 0 0 0 0 0
0 1 1 1 1 1 1 1 0 0
0 0 0 0 1 1 1 0 0 0
0 0 0 0 1 1 0 0 0 0
0 0 0 0 1 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
"""

  cx = flex.int([1,5,5,8])
  cy = flex.int([4,1,7,4])
  N=10
  mask_gen = data_reduction.image_mask(N,N)
  mask_gen.add_polygon(cx,cy,0,0)
  mask = mask_gen.mask()
  count=0
  result = """"""
  for ii in xrange(N):
    for jj in xrange(N):
      result+=str(mask[count])+" "
      count +=1
    result=result[0:len(result)-1]
    result+="\n"
  assert desired_mask == result


  desired_mask="""0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0
0 0 0 1 1 1 1 1 0 0
0 0 0 1 1 1 1 1 0 0
0 0 1 1 1 1 1 1 1 0
0 0 0 1 1 1 1 1 0 0
0 0 0 1 1 1 1 1 0 0
0 0 0 0 0 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0
"""
  mask_gen = data_reduction.image_mask(N,N)
  mask_gen.add_circle( 5,5,3 )
  mask = mask_gen.mask()
  count=0
  result = """"""
  for ii in xrange(N):
    for jj in xrange(N):
      result+=str(mask[count])+" "
      count +=1
    result=result[0:len(result)-1]
    result+="\n"
  assert desired_mask == result

  desired_mask = """1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 0 1 1 1 1
1 1 1 0 0 0 0 0 1 1
1 1 1 0 0 0 0 0 1 1
1 1 0 0 0 0 0 0 0 1
1 1 1 0 0 0 0 0 1 1
1 1 1 0 0 0 0 0 1 1
1 1 1 1 1 0 1 1 1 1
1 1 1 1 1 1 1 1 1 1
"""

  mask_gen.mask_invert()
  mask = mask_gen.mask()
  count = 0
  result = """"""
  for ii in xrange(N):
    for jj in xrange(N):
      result+=str(mask[count])+" "
      count +=1
    result=result[0:len(result)-1]
    result+="\n"
  assert desired_mask == result


def tst_data_processing():
  cx = flex.int([1,2,0])
  cy = flex.int([1,2,0])
  N=3
  mask = data_reduction.image_mask(N,N)
  mask.add_polygon(cx,cy,0,0)
  detector = data_reduction.perpendicular_detector(N,N,1,1,0,0,100.0,1.0)
  qr = detector.q_range()
  assert approx_equal(qr[0],0,1e-3)
  assert approx_equal(qr[1],0.177,1e-3)
  detector.compute_bin_indices(0.05, 0.0, 0.20);
  integrator = data_reduction.radial_integrator( detector, mask.mask()*0, None )
  data = []
  for xx in xrange(N):
    for yy in xrange(N):
      d = int( math.sqrt(xx*xx+yy*yy)*100 )
      data.append( d +10)
  data = flex.int( data )

  integrator.integrate( data, False )
  mi = integrator.mean_intensity()
  vi = integrator.variance_intensity()
  target_mi = [9.9999999999990017, 123.66666666666255, 221.49999999999443, 291.99999999997084, 0.0, 0.0 ]
  target_vi = [9.9760200100718066e-12, 373.55555555605133, 132.25000000122964, 8.5128704085946083e-09, 0.0, 0.0]
  for o,c in zip(mi,target_mi):
    assert approx_equal(o,c,1e-5)
  for o,c in zip(vi,target_vi):
    assert approx_equal(o,c,1e-5)

  print list(mi)

  integrator.integrate( data, False )
  mi = integrator.mean_intensity()
  vi = integrator.variance_intensity()
  target_mi = [9.9999999999990017, 123.66666666666255, 221.49999999999443, 291.99999999997084, 0.0, 0.0 ]
  target_vi = [9.9760200100718066e-12, 373.55555555605133, 132.25000000122964, 8.5128704085946083e-09, 0.0, 0.0]
  for o,c in zip(mi,target_mi):
    assert approx_equal(o,c,1e-5)
  for o,c in zip(vi,target_vi):
    assert approx_equal(o,c,1e-5)

  print list(mi)

  integrator.integrate( data, False )
  mi = integrator.mean_intensity()
  vi = integrator.variance_intensity()
  target_mi = [9.9999999999990017, 123.66666666666255, 221.49999999999443, 291.99999999997084, 0.0, 0.0 ]
  target_vi = [9.9760200100718066e-12, 373.55555555605133, 132.25000000122964, 8.5128704085946083e-09, 0.0, 0.0]
  for o,c in zip(mi,target_mi):
    assert approx_equal(o,c,1e-5)
  for o,c in zip(vi,target_vi):
    assert approx_equal(o,c,1e-5)

  print list(mi)

  integrator.integrate( data, False )
  mi = integrator.mean_intensity()
  vi = integrator.variance_intensity()
  target_mi = [9.9999999999990017, 123.66666666666255, 221.49999999999443, 291.99999999997084, 0.0, 0.0 ]
  target_vi = [9.9760200100718066e-12, 373.55555555605133, 132.25000000122964, 8.5128704085946083e-09, 0.0, 0.0]
  for o,c in zip(mi,target_mi):
    assert approx_equal(o,c,1e-5)
  for o,c in zip(vi,target_vi):
    assert approx_equal(o,c,1e-5)

  print list(mi)



def tst_pickling():
  # first check detector pickling
  N=20
  detector = data_reduction.perpendicular_detector(N,N,1,1,0,0,100.0,1.0)
  #detector.compute_bin_indices(0.05, 0.0, 0.20)
  s = pickle.dumps( detector )
  l = pickle.loads( s )
  sl =pickle.dumps( l )
  assert sl == s

  ipp = open('dp.pckl','w')
  dp = pickle.dump(detector,ipp)
  ipp.close()
  ipp = open('dp.pckl','r')
  nd = pickle.load( ipp )
  ipp.close()
  for a,b in zip(detector.q(), nd.q()):
     assert approx_equal(a,b,eps=1e-8)


if __name__ =="__main__":
  tst_mask_definitions()
  tst_data_processing()
  tst_pickling()
  print "OK"
