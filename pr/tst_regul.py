from scitbx.array_family import flex
import sastbx.pr
from libtbx.test_utils import approx_equal
import pr_tools
import math




def sphere_data(q, d_max):
  r=d_max/2.0
  a = r*r*r*( flex.sin(q*r)-q*r*flex.cos(q*r) )/( r*r*r*q*q*q )
  return a*a




def sphere_pofr(r_in, d_max):
  r = r_in/d_max
  p = 6.0*r*r*(2-3.0*r+r*r*r)
  return p




def tst_integration(d_max=70.0):
  q1 = 0.35*flex.double( range(1,100) )/100.0
  q2 = 0.25*flex.double( range(1,73) )/100.0

  i1 = sphere_data(q1, d_max)
  i2 = sphere_data(q2, d_max)
  i1 = 1000.0*i1/i1[0]
  i2 = 1000.0*i2/i2[0]

  s1 = i1*0.001
  s2 = i2*0.001

  d1 = q1[1]-q1[0]
  d2 = q1[1]-q1[0]
  regul_obj = sastbx.pr.regul_basic(d_max,70, max(flex.max(q1), flex.max(q2) )+0.05, min(d1,d2)/2.0)


  regul_obj.load_data(q1,i1,s1)
  regul_obj.load_data(q2,i2,s2)

  r = regul_obj.r()
  r = r
  p = sphere_pofr(r,d_max)

  regul_obj.load_scales( flex.double([1.0,1.0]) )
  regul_obj.load_p( p )

  ic1 = regul_obj.get_intensity(0)
  ic2 = regul_obj.get_intensity(1)

  ic1 = 1000.0*ic1/ic1[0]
  ic2 = 1000.0*ic2/ic2[0]

  for qq,ii,jj in zip(q1,i1,ic1):
    li = math.log(ii)
    lj = math.log(jj)
    d = abs(li-lj)/abs(li)
    assert approx_equal(d,0,eps=1e0)

  for qq,ii,jj in zip(q2,i2,ic2):
    li = math.log(ii)
    lj = math.log(jj)
    d = abs(li-lj)/abs(li)
    assert approx_equal(d,0,eps=1e0)





def test_gradients_to_data(d_max=70):
  q1 = 0.35*flex.double( range(1,100) )/100.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01+1.0
  regul_obj = sastbx.pr.regul_basic(d_max,int(d_max), flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)

  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)

  # checking generation of delta arrays
  ic = regul_obj.get_intensity( 0 )
  delta = regul_obj.get_delta(0)
  for ii,jj,dd,ss in zip( i1,ic,delta,s1):
    assert approx_equal( (ii-jj)/(ss*ss), dd, eps=1e-4 )

  scale = i1[0]/ic[0]
  regul_obj.load_scales( flex.double([scale]) )
  delta = regul_obj.get_delta(0)
  ic = regul_obj.get_intensity( 0 )
  for ii,jj,dd,ss in zip( i1,ic,delta,s1):
    assert approx_equal( (ii-jj)/(ss*ss), dd, eps=1e-4 )

  scale = flex.sum(i1)/flex.sum(ic)
  regul_obj.load_scales( flex.double([scale]) )
  delta = regul_obj.get_delta(0)
  ic = regul_obj.get_intensity( 0 )
  for ii,jj,dd,ss in zip( i1,ic,delta,s1):
    assert approx_equal( (ii-jj)/(ss*ss), dd, eps=1e-4 )
  #

  scale = 3000
  regul_obj.load_scales( flex.double([scale]) )
  delta = regul_obj.get_delta(0)
  ic = regul_obj.get_intensity( 0 )
  for ii,jj,dd,ss in zip( i1,ic,delta,s1):
    assert approx_equal( (ii-jj)/(ss*ss), dd, eps=1e-4 )
  #


  # test first derivatives
  h = 1e-8
  scale = 1.0*i1[0]/ic[0]
  regul_obj.load_p( p_theory )
  regul_obj.load_scales( flex.double([scale]) )
  dp = regul_obj.dt_dp(0)
  ic = regul_obj.get_intensity( 0 )
  t = flex.sum( (i1-ic)*(i1-ic)/(s1*s1) )
  for ii in range(len(p_theory)):
    new_p = p_theory.deep_copy()
    new_p[ii] = p_theory[ii]+h
    regul_obj.load_p( new_p )
    new_ic = regul_obj.get_intensity( 0 )
    th = regul_obj.t_int(0)
    dd = (th-t)/h
    assert approx_equal( 1, dd/dp[ii],eps=1e-1 )


  scale = 10.0*i1[0]/ic[0]
  regul_obj.load_p( p_theory )
  regul_obj.load_scales( flex.double([scale]) )
  dp = regul_obj.dt_dp(0)
  ic = regul_obj.get_intensity( 0 )
  t = flex.sum( (i1-ic)*(i1-ic)/(s1*s1) )
  for ii in range(len(p_theory)):
    new_p = p_theory.deep_copy()
    new_p[ii] = p_theory[ii]+h
    regul_obj.load_p( new_p )
    new_ic = regul_obj.get_intensity( 0 )

    th = regul_obj.t_int(0)
    dd = (th-t)/h
    assert approx_equal( 1.0,dp[ii]/dd,eps=1e-3 )


  scale = 100.0*i1[0]/ic[0]
  regul_obj.load_p( p_theory )
  regul_obj.load_scales( flex.double([scale]) )
  dp = regul_obj.dt_dp(0)
  ic = regul_obj.get_intensity( 0 )
  t = flex.sum( (i1-ic)*(i1-ic)/(s1*s1) )
  for ii in range(len(p_theory)):
    new_p = p_theory.deep_copy()
    new_p[ii] = p_theory[ii]+h
    regul_obj.load_p( new_p )
    new_ic = regul_obj.get_intensity( 0 )

    th = regul_obj.t_int(0)
    dd = (th-t)/h
    assert approx_equal( 1.0,dd/dp[ii],eps=1e-1 )


def test_hessian_to_data(d_max=10,scale=0.1):
  q1 = 0.35*flex.double( range(1,100) )/100.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01+1.0
  regul_obj = sastbx.pr.regul_basic(d_max,int(d_max), flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)
  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)
  regul_obj.load_scales( flex.double([scale]) )

  hessian = regul_obj.hessian_t(0)
  count=0


  dt_dp = regul_obj.dt_dp(0)
  h = 0.0001
  fd_hes= []
  for ii in range(r.size()):
    for jj in range(r.size()):
      new_p = p_theory.deep_copy()
      new_p[ii] = p_theory[ii]+h
      regul_obj.load_p( new_p )
      new_dt_dp = regul_obj.dt_dp(0)
      delta_h = (new_dt_dp - dt_dp)/h
      index = ii+jj*(r.size())
      if ii == jj:
        fd_hes.append( delta_h )


  count=0
  for ii in range(r.size()):
    for jj in range(r.size()):
      assert approx_equal( fd_hes[jj][ii], hessian[count], 1e-3 )
      count += 1


def tst_scales(d_max=80,scale=1):
  q1 = 0.50*flex.double( range(1,200) )/200.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01 + 1.0
  regul_obj = sastbx.pr.regul_basic(d_max,int(d_max), flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)
  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)
  regul_obj.load_scales( flex.double([scale]) )

  inten = regul_obj.get_intensity(0)

  tt = (i1-regul_obj.get_intensity(0))
  tt = tt*tt/(s1*s1)
  tt = flex.sum(tt)
  t=regul_obj.t_int(0)
  dt_ds = regul_obj.dt_ds(0)
  assert approx_equal(tt, t, 1e-3 )

  h=1e-7
  new_scale = scale+h
  regul_obj.load_scales( flex.double([new_scale]) )
  th = regul_obj.t_int(0)
  fd_dt_ds = (th-t)/h
  assert approx_equal( dt_ds,fd_dt_ds, eps=1e-3 )
  ddd = regul_obj.dt_ds(0)
  assert approx_equal( (ddd-dt_ds)/h, regul_obj.ddt_dds(0), eps=1e-2)


  scale = scale
  for ii in range(3):
    g = regul_obj.dt_ds(0)
    G = regul_obj.ddt_dds(0)
    d = -g/G
    scale = scale+d
    regul_obj.load_scales( flex.double([scale]) )
  assert approx_equal( regul_obj.t_int(0), 0, eps=1 )


def tst_regul(d_max=30,scale=1.0):
  q1 = 0.5*flex.double( range(10,200) )/200.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01
  regul_obj = sastbx.pr.regul_basic(d_max,30, flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)*0+1.0
  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)
  assert approx_equal( regul_obj.t_regul(), 1.0, eps=1e-5)



def tst_dregul(d_max=10,scale=1.0):
  q1 = 0.5*flex.double( range(10,200) )/200.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01
  regul_obj = sastbx.pr.regul_basic(d_max, int(d_max), flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)
  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)

  dr_dp = regul_obj.dregul_dp()
  t = regul_obj.t_regul()
  h=1e-6
  for ii in range(int(d_max)):
    new_p    = p_theory.deep_copy()
    new_p[ii] = new_p[ii] + h
    regul_obj.load_p( new_p )
    th = regul_obj.t_regul()
    assert approx_equal( (th-t)/h / dr_dp[ii], 1.0, eps=1e-3)

  regul_obj = sastbx.pr.regul_basic(7, 7, flex.max(q1)+0.05, 0.01)
  hess = regul_obj.hessian_regul()
  count=0



  regul_obj.load_p( p_theory )
  dr_dp = regul_obj.dregul_dp()
  h = 0.0001
  fd_hes= []

  for ii in range(r.size()):
    for jj in range(r.size()):
      new_p = p_theory.deep_copy()
      new_p[ii] = p_theory[ii]+h
      regul_obj.load_p( new_p )
      new_dr_dp = regul_obj.dregul_dp()
      delta_h = (new_dr_dp - dr_dp)/h
      if ii == jj:
        fd_hes.append( delta_h )



  count=0


  for ii in range(r.size()):
    for jj in range(r.size()):
      assert approx_equal( fd_hes[ii][jj], hess[count], eps=1e-3)
      count += 1



def tst_refi(d_max=45,scale=1.0):
  q1 = 0.5*flex.double( range(10,200) )/200.0
  i1 = sphere_data(q1, d_max)
  i1 = 1000.0*i1/i1[0]
  s1 = i1*0.01
  regul_obj = sastbx.pr.regul_basic(d_max,20, flex.max(q1)+0.05, 0.01)
  r = regul_obj.r()
  p_theory = sphere_pofr(r,d_max)
  regul_obj.load_p( p_theory )
  regul_obj.load_data(q1,i1,s1)
  i1 =  regul_obj.get_intensity(0)
  regul_obj = sastbx.pr.regul_basic(d_max,20, flex.max(q1)+0.05, 0.01)
  regul_obj.load_data(q1,i1,s1)
  p_theory = sphere_pofr(r,d_max)*0+1
  regul_obj.load_p( p_theory )
  regul_obj.load_scales( flex.double([scale]) )
  jc = regul_obj.get_intensity(0)

  # get the scales right please
  for ii in range(3):
    g = regul_obj.dt_ds(0)
    G = regul_obj.ddt_dds(0)
    d = -g/G
    scale = scale+d
    regul_obj.load_scales( flex.double([scale]) )
  print regul_obj.t_int(0), scale
  print list(p_theory)

  grad = regul_obj.dt_dp(0)
  hess = regul_obj.hessian_t(0)


  import scitbx.matrix

  scalar=10.0

  for n_rounds in range(1000):
    for ii in range(2):
      print ii, regul_obj.t_int(0)
      grad = regul_obj.dt_dp(0)
      hess = regul_obj.hessian_t(0)
      g = scitbx.matrix.row( grad )
      msq = scitbx.matrix.sqr( hess )
      hinv = scitbx.matrix.inverse_via_lu( msq )
      d = g*hinv
      d = flex.double(d)
      d = d/(flex.max(flex.abs(d)))
      new_p = p_theory-d/scalar
      regul_obj.load_p( new_p )
      p_theory = new_p.deep_copy()

    for ii in range(2):
      g = regul_obj.dt_ds(0)
      G = regul_obj.ddt_dds(0)
      d = -g/G
      scale = scale+d
      regul_obj.load_scales( flex.double([scale]) )
      print regul_obj.t_int(0), scale, "<---", n_rounds

  print list(p_theory)
  ic = regul_obj.get_intensity(0)
  for qq,ii,jj,kk in zip(q1,i1,ic,jc):
    print qq,ii,jj,kk





if __name__ == "__main__":
  tst_integration()
  test_gradients_to_data()
  test_hessian_to_data()
  test_hessian_to_data(d_max=20,scale=10)
  tst_scales(d_max=60, scale=0.0001)
  #tst_refi(d_max=60, scale=1.0)
  tst_regul()
  tst_dregul(7)

  print "OK"
