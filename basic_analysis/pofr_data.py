from scitbx.array_family import flex
import sys

class pofr(object):
  def __init__(self, filename=None, r=None, pr=None, error=None, from_gnom=True):
    self.cdf = flex.double()
    if(filename is None):
      self.r = r.deep_copy()
      self.pr = pr.deep_copy()
      self.dr = self.r[1] - self.r[0]
      self.pr = self.pr / ( self.dr * flex.sum( self.pr ) )
      if(error is None):
        self.error = self.pr * 0
      else:
        self.error = error.deep_copy()
    else:
      self.r = flex.double()
      self.pr = flex.double()
      self.error = flex.double()
      if(from_gnom):
        self.read_pr_from_gnom( filename )
      else:
	self.read_pr_from_std( filename )

  def read_pr_from_std(self, filename):
    file = open( filename, 'r' )
    for line in file:
      keys = line.split('\n')[0].split()
      if(len(keys) <= 3 and len(keys)>1):
        try:
          self.r.append( float(keys[0]) )
          self.pr.append( float(keys[1]) )
        except ValueError:
          break
    file.close()
    sum = flex.sum( self.pr ) * ( self.r[1] - self.r[0] )
    self.pr = self.pr / sum
 
  def read_pr_from_gnom(self,filename):
    file=open(filename, 'r')
    pr_start = False
    for line in file:
      keys = line.split('\n')[0].split()
      if(len(keys) == 3):
        if(not pr_start):
          if(keys[0] == 'R' and (keys[1] == 'P(R)') and (keys[2] == 'ERROR') ):
            pr_start=True
        else:
          self.r.append( float(keys[0]) )
          self.pr.append( float(keys[1]) )
          self.error.append( float(keys[2]) )
    file.close()
    sum = flex.sum( self.pr ) * ( self.r[1] - self.r[0] )
    self.pr = self.pr / sum
    self.error = self.error / sum

  def print_data(self, data):
    print '&'
    for ri,pri in zip(self.r,data):
      print ri, pri

  def linear_interpolation(self,new_r):
    if(self.r[-1]<new_r[-1]):
      print "WARNING: xx Interpolation out of bound xx "
      self.r.append(new_r[-1])
      self.pr.append(0)
    new_pr = flex.linear_interpolation(self.r, self.pr, new_r)
    return new_pr/(flex.sum(new_pr) * (new_r[1] -new_r[0]) )

  def pr2cdf(self):
    new_cdf = flex.double()
    dr = self.r[1]-self.r[0]
    sum = self.pr[0]
    for pri in self.pr[1:]:
      new_cdf.append( sum )
      sum += pri
    new_cdf.append( sum )
    self.cdf = new_cdf*dr
    return self.cdf
  
  def cdf2pr(self):
    #self.pr = flex.double()
    #self.pr.append( self.cdf[0] )
    #for ii in range(self.cdf.size()-1):
    #  self.pr.append( self.cdf[ii+1] - self.cdf[ii] )
    self.pr = self.cdf[1:] - self.cdf[:-1]
    self.pr.insert(0,0)
    #self.pr.insert(0,self.cdf[0])
    dr = self.r[1]-self.r[0]
    self.pr = self.pr /( flex.sum( self.pr )*dr )
    return self.pr

def test(filename):
  data = pofr(filename)
  data.print_data(data.pr)
  new_r = flex.double( range( int(max(data.r)+0.5)*2)) /2.0
  new_pr = data.linear_interpolation(new_r)
  new_data=pofr(r=new_r, pr=new_pr)
  print '&'
  new_data.print_data(new_data.pr)
  print '&'
  new_data.pr2cdf()
  new_data.print_data( new_data.cdf )


if __name__ == "__main__":
  test(sys.argv[1])
