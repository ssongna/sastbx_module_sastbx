import os, sys
from scitbx.array_family import flex
import math, time
import random
from sastbx.data_reduction import saxs_read_write
from sastbx.rigid_body import rigidbody as rb
from sastbx.rigid_body import rb_engine as rbe


def polar2xyz( theta, phi, r=1):
  sint = math.sin( theta )
  cost = math.cos( theta )
  sinp = math.sin( phi )
  cosp = math.cos( phi )
  x = cosp * sint
  y = sinp * sint
  z = cost
  if(r==1): 
    return [ x, y , z ]
  else:
    return [ r*x, r*y, r*z ]

def topYandX(x_array,y_array, top_n):
  assert (len(x_array) == len(y_array))
  Y_top = []
  X_top = []
  X_top_indices = []
  
  for ii in range( len( x_array ) ):
    inserted = False
    for jj in range( len(X_top) ):
      if( x_array[ii] <= X_top[jj] ):
        X_top.insert(jj, x_array[ii] )
	X_top_indices.insert(jj, ii)
        if( len(X_top) > top_n ):
      	  X_top.pop()
	inserted = True
        break
    if((not inserted) and (len(X_top) <top_n) ):
      X_top.append(x_array[ii])
      X_top_indices.append(ii)

  for ii in range( top_n ):
    Y_top.append( y_array[ X_top_indices[ii] ] )
  return X_top, Y_top


class grid(object):
  def __init__(self, grids, r_range, r_step=1.0, delta_r=3.0):
    self.grids = grids
    self.r_range = r_range
    self.r_step = r_step
    self.delta_r = delta_r
    self.xyz_grid = []
    self.build_grid()

  def build_grid(self):
    for g in self.grids:
      base_xyz = flex.double( polar2xyz( g[0], g[1] ) )
      r = g[2]
      r_max = r + self.r_range + self.delta_r
      while(r < r_max):
        self.xyz_grid.append( list( base_xyz * r ) )
        r += self.r_step





class grid_sample(object):
  def __init__(self,rbe_obj,expt_data, top_n=200):
    self.data=expt_data
    self.pr = expt_data.i
    self.rbe=rbe_obj
    self.nbody=rbe_obj.nbody
    self.top_n = top_n

    self.main_grid = self.rbe.rbs[0].get_grid()
    self.r_range = self.rbe.rbs[1].range_f_w()
    self.delta_r = 3.0

    self.xyz_grid=grid(self.main_grid, self.r_range, self.delta_r).xyz_grid

    self.scores = []
    self.sampling()
    self.top_scores, self.top_solutions = topYandX(self.scores, self.xyz_grid, self.top_n)


  def sampling(self):
    for xyz in self.xyz_grid:
        score = self.target(xyz)
	self.scores.append( score )

  def target(self, vector):
    self.rbe.rbs[1].translate_after_rotation(vector)
    calc_pr = self.rbe.get_norm_pr()
    score = flex.sum( flex.abs(calc_pr-self.pr)*flex.exp(-self.pr) )
    return score



#=================================================================================================#
#
#=================================================================================================#

def test(args):
  rbs=[]
  file=args[0]
  expt_data=saxs_read_write.read_standard_ascii_qis( file )
  dmax=expt_data.q[-1]+1
  pdb=None
  max_num_fibonacci = 17
  count = 0
  group_size = 2 # this determines the size of each atom-group: the larger the size, the faster the calculation
  main_body = False
  for arg in args[1:]:
    pdb= rbe.PDB(arg)
    pdb.CA_indx = flex.int( range(0, pdb.xyz.size(), group_size ) )
    if(count == 0): 
      mainbody = True
    rbs.append(rb(pdb.xyz, pdb.CA_indx, dmax, max_num_fibonacci, main_body))
    count += 1

  target_xyz = rbs[1].get_crd().deep_copy();


  rb_eng = rbe.rb_engine(rbs,int(dmax) )

  top_n =100
  sample = grid_sample( rb_eng, expt_data, top_n=top_n)

  for ii in range(top_n):
    rbs[1].translate_after_rotation( sample.top_solutions[ii] )
    new_xyz = rbs[1].get_crd()
    outname = "refined"+str(ii)+".pdb"
    pdb.writePDB(new_xyz, outname)
    RMSD = target_xyz.rms_difference( new_xyz )
    for x in sample.top_solutions[ii]:
      print x,
    print RMSD


def test_topYandX():
  x_array = [2, 3, 1, 5, 8, 10, -1]
  y_array = [1, 2, 3, 4, 5, 6, 7]

  top_n = 3
  top_x, top_y = topYandX( x_array, y_array, top_n )
  print x_array
  print y_array
  print top_x
  print top_y

#==============================================================================================#
#
#==============================================================================================#
if __name__ == "__main__":
  t1 = time.time()
  test(sys.argv[1:])
#  test_topYandX()

  t2 = time.time()
  print "sampling time", t2-t1

