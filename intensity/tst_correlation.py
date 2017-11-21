from sastbx.intensity import she
from scitbx.array_family import flex
import sys
from stdlib import math as smath

def unpack( correlation, n_phi, n_q ):
  new_co = []
  index = 0
  for pp in range(n_phi):
    tmp_co = []
    for qq in range( n_q ):
      tmp_co_q = []
      for qqq in range(qq+1):
	tmp_co_q.append( correlation[index] )
	index = index + 1
      tmp_co.append( tmp_co_q )
    new_co.append( tmp_co )
  return new_co

def print_coef(coefs):
  for ll, bb in zip( range(expansion_coef.size()), expansion_coef):
    print ll, bb

def print_correlation( out_correlation, correlation, q_array, N_phi):
  print >>out_correlation, "# phi, q, correlation"
  two_pi = smath.pi*2.0
  ii = 0
  for qq in q_array:
    for pp in range(N_phi):
      phi = pp*two_pi/N_phi
      print >>out_correlation, qq*smath.cos(phi),qq*smath.sin(phi), correlation[pp][ii][ii]
    ii=ii+1


def run(pdb_file,q_array,output_prefix):
  N_phi = 51
  phi_array = flex.double( range(N_phi) ) / 25.0 * 3.1416
  she_obj = she.she(pdb_file, q_array=q_array, max_L=15)


  she_obj.engine.calc_spatial_correlation(phi_array)
  correlation = she_obj.engine.get_spatial_correlation()
  correlation = unpack( correlation, N_phi, q_array.size() )
  
  out_intensity = open(output_prefix+".int", 'w')
  out_correlation = open(output_prefix+".cor",'w')
  out_coef = open(output_prefix+".coe", 'w')

  print >>out_coef, "#q", 
  for ii in range(7):
    print >>out_coef, ii*2,
  print >>out_coef, "..."

  for qq in range(q_array.size() ):
    expansion_coef = she_obj.engine.get_expansion_coef(qq)
    print >>out_coef,qq*0.01, 
    for cc in expansion_coef:
      print >>out_coef, cc,
    print >>out_coef
  intensity = she_obj.engine.I()
  intensity_vac = she_obj.engine.get_IA()
  print >>out_intensity,'# q, intensity, intensity_in_vac'
  for qq,ii,ia in zip(q_array, intensity, intensity_vac):
    print >>out_intensity, qq,ii,ia

  #  print '&', qq
  #  print_coef( expansion_coef )
    
  print_correlation( out_correlation, correlation, q_array, N_phi)

  out_intensity.close()
  out_correlation.close()
  out_coef.close()



if __name__ == "__main__":
  q_array=flex.double(range(50))/100.0
  args = sys.argv[1:]
  if(len(args) != 2):
    print "Usage: command pdb_file prefix"
    exit()
  run(args[0], q_array, args[1])
