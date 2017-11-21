from scitbx.array_family import flex
import sys,os
from stdlib import math as smath
from libtbx.utils import n_dim_index_from_one_dim 

def enlarge(data,factor,sigma=0.1,full=False):
  n,n = data.focus()
  m = int(n*factor)
  x    = flex.double()
  y    = flex.double()
  vals = flex.double()
  sigma=sigma*factor 
  new_data = flex.double( flex.grid(m,m), -9 )
  visited = flex.bool( flex.grid(m,m), False )
  oo = n/2.0
  for ii in range(n):
    for jj in range(n):
      dd = smath.sqrt( (ii-oo)**2.0 + (jj-oo)**2.0 )
      if dd <= oo:
        nx = ii*factor
        ny = jj*factor
        x.append( nx )
        y.append( ny )
        vals.append( data[ (ii,jj) ] )
        new_data[ (int(nx), int(ny)) ] = data[ (ii,jj) ]
        if not full:
          visited[ (int(nx), int(ny)) ] = True

  
  not_visited = ~visited
  not_visited = not_visited.iselection()
  # now we need to loop over all non-visited pixel values
  for pixel in not_visited:
        nv = -9
        index =  n_dim_index_from_one_dim(pixel, [m,m] )
        nvx  = index[1]
        nvy  = index[0]
        dx   = x-nvx
        dy   = y-nvy
        dd   = (dx*dx+dy*dy)
        ss   = flex.exp( -dd/(sigma*sigma) )
        nv   = flex.sum(ss*vals)/(1e-12+flex.sum(ss))
        new_data[ (nvx,nvy) ] = nv
        visited[  (nvx,nvy) ] = True
        #print nvx, nvy, nv
  not_visited = ~visited
  not_visited = not_visited.iselection()
  #print not_visited.size()
  return new_data
  oo=m/2.0
  for ii in range(m):
    for jj in range(m):
      dd = smath.sqrt( (ii-oo)**2.0 + (jj-oo)**2.0 )
      new_data[ (ii,jj) ] = new_data[ (ii,jj) ]/(1+smath.sqrt(dd))
  return new_data
      

def rod(np, length, width, val=1.0):
  result = flex.double( flex.grid(np,np), 0 )
  assert length < np
  assert width < length
  startx = (np-length)/2
  starty = (np-width)/2 
  for ii in range(length):
    for jj in range(width):
       result[ (ii+startx,jj+starty) ] = val
  return result

def ball(np, radius,val=1.0):
  result = flex.double( flex.grid(np,np), 0 )
  assert radius < np
  o1 = np/2.0-radius/1.5
  o2 = np/2.0+radius/1.5
  for ii in range(np):
    for jj in range(np):
      dd = (ii-o1)**2.0 + (jj-o1)**2.0
      if dd <= radius**2.0:
        result[ (ii,jj) ] = val
      dd = (ii-o2)**2.0 + (jj-o2)**2.0
      if dd <= radius**2.0:
        result[ (ii,jj) ] = val

  return result

def single_ball(np, radius,val=1.0):
  result = flex.double( flex.grid(np,np), 0 )
  assert radius < np
  o1 = np/2.0
  for ii in range(np):
    for jj in range(np):
      dd = (ii-o1)**2.0 + (jj-o1)**2.0
      if dd <= radius**2.0:
        result[ (ii,jj) ] = val

  return result

