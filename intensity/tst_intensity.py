from sastbx import intensity
from scitbx.array_family import flex
def tst_sphere():
  bli = intensity.block_integrator(50.0, 0.1, 1.0, 0.001)
  q = flex.double( range(1,300) )/1000.0
  bli.setup_arrays( q )
  r = flex.double( range(500) )/500.0
  r = r
  pr = 6.0*r*r*(2-3.0*r+r*r*r)
  r = r*50

  #for rr, pp in zip(r,pr):
  #  print rr, pp

  iii = bli.get_intensity( pr )
  iii = iii/25.0

  jjj = ( flex.sin( 25.0 *q) - 25.0*q*flex.cos(q*25.0) )/( 1e-12+q*q*q*25*25*25 )
  jjj = jjj*jjj
  jjj = jjj/jjj[0]
  for q,i,j in zip(q,iii,jjj):
    assert( abs(i-j)<1e-2)






tst_sphere()
