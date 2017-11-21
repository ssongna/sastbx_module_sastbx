from scitbx.array_family import flex


class cluster(object):
  def __init__(self,pair_dist,cutoff):
    self.cutoff = cutoff
    self.pdist = pair_dist
    self.clusters=[]
    self.n = pair_dist[0].size()
    self.do_cluster2()
    self.merge(self.cutoff)

  def do_cluster(self):
    cluster=flex.int([0])
    self.clusters.append(cluster)
    for ii in range(1,self.n):
      classified = False
      for ci in range(len(self.clusters)):
        for jj in self.clusters[ci]:
          #print ii,ci,jj
          if(self.pdist[ii][jj] < self.cutoff):
            self.clusters[ci].append(ii)
	    classified = True
	    break
      if (not classified):
        cluster=flex.int([ii])
        self.clusters.append(cluster)

  def do_cluster2(self): # to the center of the cluster
    cluster=flex.int([0])
    self.clusters.append(cluster)
    for ii in range(1,self.n):
      classified = False
      for ci in range(len(self.clusters)):
        dist = 0.0
        for jj in self.clusters[ci]:
          dist += self.pdist[ii][jj]
        if( dist / float(len(self.clusters[ci])) < self.cutoff):
          self.clusters[ci].append(ii)
          classified = True
      if (not classified):
        cluster = flex.int([ii])
        self.clusters.append(cluster)

  def merge(self, cutoff):
    self.centers=[]
    for cluster1 in self.clusters:
      c2_index = 0
      for cluster2 in self.clusters:
        c2_index += 1
        if(cluster1 is not cluster2):
          dist = 0
          count = 0
          for ii in cluster1:
            for jj in cluster2:
              dist += self.pdist[ii][jj]
              count += 1.0
          dist = dist/count
          if(dist < cutoff):
            for jj in cluster2:
              cluster1.append(jj)
            self.clusters.pop(c2_index)

            

  def largest(self):
    max_n = 0
    second_max_n = 0
    largest = self.clusters[0]
    second = largest
    for c in self.clusters:
      if max_n < len(c):
	second = largest
	second_max_n = max_n
        max_n = len(c)
        largest = c
      elif second_max_n < len(c):
	second = c
	second_max_n = len(c)
    return largest, second

  def print_cluster(self):
    for c in self.clusters:
      print list(c)

def test():
   dist_mat = []
   n = 5
   for ii in range(n):
     dist_mat.append(flex.random_double(n))
   for ii in range(n):
     dist_mat[ii][ii] = 0
     for jj in range(ii):
       dist_mat[ii][jj]=dist_mat[jj][ii]

     print list(dist_mat[ii])

   cutoff = 0.2

   clusters = cluster(dist_mat,cutoff)
   clusters.print_cluster()
   largest = clusters.largest()
   print list(largest)

  

if __name__ == "__main__":
  test()
