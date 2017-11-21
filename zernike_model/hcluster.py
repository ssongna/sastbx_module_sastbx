from scitbx.array_family import flex


class node(object):
  def __init__(self, node0, node1=None, dist=0):
    if(node1 is None):
      self.elements = [node0]
      self.n_ele = 1
      self.leaf_eles = [ node0 ]
      self.height=1
    else:
      self.elements = [node0, node1]
      self.n_ele = node0.n_ele + node1.n_ele
      self.leaf_eles = node0.leaf_eles + node1.leaf_eles
      self.height = max( node0.height, node1.height ) + 1
    self.dist = dist

  def get_elements(self):
    return self.leaf_eles

  def print_node(self, height, out):
    space = "______"
    string = ""
    if( self.height == 1 ):
      print str(self.leaf_eles[0]+1),
    else:
      print "(",
      for nn in self.elements:
        nn.print_node(height-1,out),
      print ")",
   
   


  def print_dot(self, my_id, parent_id, output):
    my_id +=1
    if( self.height == 1):
      print >>output, "  subgraph %4d {\n  %4d [shape=box]\n  }"%(my_id, self.leaf_eles[0]+1) # make the node a box
      print >>output, "  %4d -> %4d [dir = none]"%(parent_id, self.leaf_eles[0]+1)
    else:
      print >>output, "  subgraph %4d {\n  %4d [shape=point]\n  }"%(my_id, my_id)  # make the node a point
      print >>output, "  %4d -> %4d [dir = none]"%(parent_id, my_id)
      parent_id=my_id
      for nn in self.elements:
        my_id +=1
        nn.print_dot(my_id, parent_id, output)
    return my_id



class hcluster(object):
  def __init__(self,pair_dist,cutoff,outfile):
    self.cutoff = cutoff
    self.pdist = pair_dist
    self.current_dmat=[]
    self.outfile = outfile
    for vector in self.pdist:
      self.current_dmat.append( vector.deep_copy() )

    self.nodes=[]
    self.n = pair_dist[0].size()
    for ii in range(self.n):
      self.nodes.append( node(ii) )

    self.merge(self.cutoff)

  def merge(self, cutoff):
    need_merge=True
    while( need_merge ):
      need_merge = self.merge_closest()

  def merge_closest(self):
    closest = 0
    n = len(self.nodes)
    for ii in range(n):
      for jj in range(ii):
        this_value = self.current_dmat[ii][jj]
        if( this_value > closest ):
          closest=this_value
          merge_pair = (ii,jj)

    if (closest > self.cutoff ):
      ii = merge_pair[0]
      jj = merge_pair[1]
      new_node = node( self.nodes[ii], self.nodes[jj], closest )
      self.nodes.pop(ii)
      self.nodes.pop(jj)
      self.nodes.append( new_node )
      self.current_dmat.pop( ii )
      self.current_dmat.pop( jj )
      new_dmat = []
      for dist_vec in self.current_dmat:
        new_vec = flex.double()
        for kk in range( dist_vec.size() ):
          if( kk == ii or kk == jj ): continue
          new_vec.append( dist_vec[kk] )
        new_dmat.append( new_vec )

      self.current_dmat = new_dmat
      new_vec = flex.double()
      for nn in self.nodes[0:-1]:
        dist = self.calc_dist(nn,new_node)
        new_vec.append( dist )
      self.current_dmat.append( new_vec )
   #   print "-----------eleminated", ii,jj
   #   for vector in self.current_dmat:
   #     print list(vector)
      return True
    else:
      return False


  def calc_dist(self, node1, node2 ):
    dist = 0
    for ii in node1.leaf_eles:
      for jj in node2.leaf_eles:
        dist = max( dist, self.pdist[ii][jj] )
    return dist


  def print_hclust(self,outfile):
    with open(outfile,"a") as out:
      print >>out, "%4d elements, %4d clusters, @cutoff=%f"%(self.n, len(self.nodes), self.cutoff)
    print "%4d elements, %4d clusters, @cutoff=%f"%(self.n, len(self.nodes), self.cutoff)
    height = 1
    for nn in self.nodes:
      if(height < nn.height): height=nn.height
    for nn in self.nodes:
      
      nn.print_node(height,out)
      #### statistics for each cluster head node ####
      mean_value = 0
      max_value = 0
      min_value =1
      for ii in nn.leaf_eles:
        for jj in nn.leaf_eles:
          if (ii == jj): continue
          this_dist = self.pdist[ii][jj]
          mean_value += this_dist
          if( max_value < this_dist ): max_value = this_dist
          if( min_value > this_dist ): min_value = this_dist
      if( nn.n_ele > 1 ):
        mean_value = mean_value/((nn.n_ele-1)*nn.n_ele )
      print
      with open(outfile,"a") as out:
        print >>out, "mean_value, max_value, min_value, (max_value-min_value)"
        print >>out, mean_value, max_value, min_value, (max_value-min_value)
      print "mean_value, max_value, min_value, (max_value-min_value)"
      print mean_value, max_value, min_value, (max_value-min_value)

  def print_dot(self, out=None):
    if(out is not None):
      output = open(out, 'w')
    else:
      output = sys.stdout
    my_id = self.n + 1
    print >>output, "digraph hclust{\n subgraph %4d {\n  %4d [shape = box, color=red, label=Root] \n  }"%(my_id, my_id)
    parent_id = my_id
    for nn in self.nodes:
      my_id = nn.print_dot(my_id, parent_id, output)
    print >>output, "}"

    if(out is not None):
      output.close()




  def print_neato(self,out='neato.dot',level1=0.90, level2=0.80,level3=0.60):
    n = len(self.pdist)
    distmat=self.pdist
    cmnd = "graph G{"
    for ii in range(n):
      id1 = str(ii+1)
      tot_con=0
      for jj in range(ii+1,n):
        if ii != jj:
          id2=str(jj+1)
          dd = distmat[ii][jj]
          if dd > level1:
            tot_con+=1
            cmnd += "%s -- %s;"%(id1,id2)
          if dd < level1:
            if dd>=level2:
              tot_con+=1
              cmnd += "%s -- %s [len=1.75];"%(id1,id2)
            if dd < level2:
              if dd>level3:
                tot_con+=1
                cmnd += "%s -- %s [len=3.75,style=dotted];"%(id1,id2)
      if tot_con==0:
        cmnd+= "%s;"%id1
    cmnd+="}"
    f = open(out,'w')
    print >> f, cmnd
    f.close()


def test():
   dist_mat = []
   n = 10
   items = flex.random_double(n)
   for ii in range(n):
     dist_mat.append(flex.random_double(n))

   for ii in range(n):
     dist_mat[ii][ii] = 0
     for jj in range(ii):
       dist = abs( items[ii]-items[jj] )
       dist_mat[ii][jj]=dist_mat[jj][ii]=dist

   for ii in range(n):
    print "test"
    print ii,
    for jj in range(n):
      print "%6.3f"%dist_mat[ii][jj],
    print


   cutoff = 0.5

   clusters = hcluster(dist_mat,cutoff)

   for ii in range(n):
     print ii, items[ii]
   clusters.print_hclust()
   dotfile = 'out.dot'
   neatofile = 'neato.dot'
   clusters.print_dot(out=dotfile)
   clusters.print_neato(out=neatofile,level1=0.7,level2=0.5,level3=0.4)


if __name__ == "__main__":
  test()
