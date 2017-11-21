from scitbx import math
from scitbx.array_family import flex
from stdlib import math as smath

def good2n( nmax, coefs_list, ref_coefs, threshold=0.90,outfile='' ):

  ## cc=0.90 is equivalent to 5% mutation in real space at nmax<=10
  max_indx=math.nlm_array(nmax).nlm().size()
  for nn in range(nmax,1,-1):
    min_indx = math.nlm_array(nn-1).nlm().size()
    #coef_0 = ref_coefs[min_indx:max_indx]
    coef_0 = ref_coefs[0:max_indx]
    mean_0 = abs( ref_coefs[0] )
    sigma_0 = flex.sum( flex.norm( coef_0 ) ) - mean_0**2
    sigma_0 = smath.sqrt( sigma_0 )
    cc_array = flex.double()
    #out = open(outfile,"w")
    for coef in coefs_list:
      #coef_1 = coef[min_indx:max_indx]
      coef_1 = coef[0:max_indx]
      mean_1 = abs( coef[0] )
      sigma_1 = flex.sum( flex.norm( coef_1 ) ) -mean_1**2
      sigma_1 = smath.sqrt( sigma_1 )
      cov_01 = abs (flex.sum( coef_0*flex.conj(coef_1) ) )
      cov_01 = cov_01 - mean_0*mean_1
      this_cc = cov_01/sigma_1/sigma_0
      cc_array.append( this_cc )
      out = open(outfile,"a")
      print >>out,this_cc
      out.close()
      print this_cc
    mean_cc = flex.mean( cc_array )
    out = open(outfile,"a")
    print >>out, "level n: ", nn, mean_cc
    out.close()
    print "level n: ", nn, mean_cc
    if( mean_cc >= threshold ): return nn
    max_indx = min_indx
  return nn


def tst():
  nmax=20
  max_indx = math.nlm_array(nmax).nlm().size()
  a=flex.complex_double( flex.random_double(max_indx), flex.random_double(max_indx) )
  b=flex.complex_double( flex.random_double(max_indx), flex.random_double(max_indx) )
  c_list=[a]

  good_n = good2n( nmax, c_list, b, threshold=0.8 )
  print good_n

if __name__ == "__main__":
  tst()
