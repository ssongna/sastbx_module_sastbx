from scitbx.array_family import flex
import os



class bounds(object):
  def __init__(self,n_para,db_file):
    self.n = n_para
    self.bounds = []
    for ii in range(self.n):
      self.bounds.append(flex.double())
    if(os.path.exists(db_file)):
      self.success = self.read_bounds(db_file)
      self.assign()
      if(not self.success):
        self.read_default()
    else:
      self.read_default()
    self.assign()

  def assign(self):
    self.mean=flex.double()
    self.min=flex.double()
    self.max=flex.double()
    self.var=flex.double()
    for ii in range(self.n):
      self.min.append(self.bounds[ii][0])
      self.max.append(self.bounds[ii][1])
      self.mean.append(self.bounds[ii][2])
      self.var.append(self.bounds[ii][3])
    return


  def read_default(self):
    if(self.n > 8):
      for ii in range(self.n):
        self.bounds[ii]=flex.random_double(4)*0
    print "& default constraints are applied"
    if(self.n==4):
       self.bounds[0] = flex.double([-13.6492076687, 4.47249489477, -0.648531500588, 1.00508578883])
       self.bounds[1] = flex.double([-11.7770483671, 3.80867469991, -0.464755522251, 0.805970604608])
       self.bounds[2] = flex.double([-8.26619615161, 2.65922048539, -0.467162236722, 0.684465258322])
       self.bounds[3] = flex.double([-3.88576804298, 1.20291988782, -0.201727072497, 0.321418789576])

    if(self.n==5):  #coef 5
       self.bounds[0] = flex.double([-4.10870916383, 14.5435070049, -0.0237872167003, 1.04148458907])
       self.bounds[1] = flex.double([-4.34557410092, 13.2260007625, -0.0505775399547, 1.04786892835])
       self.bounds[2] = flex.double([-3.02921642669, 9.47967940636, -0.0255496685307, 0.73961261201])
       self.bounds[3] = flex.double([-2.43568258085, 5.45508644475, -0.0428255444417, 0.54548555217])
       self.bounds[4] = flex.double([-1.05904411075, 2.18286060437, -0.015965263936, 0.23150527771])
    if(self.n==6):  #coef 6
       self.bounds[0] = flex.double([-9.98987114005, 5.75104320728, -0.203361142977, 1.2797496454])
       self.bounds[1] = flex.double([-9.0881925762, 5.21730938463, -0.168699949349, 1.12217804255])
       self.bounds[2] = flex.double([-7.65751454288, 4.48953026585, -0.169312093081, 1.01238401591])
       self.bounds[3] = flex.double([-5.15246914011, 3.06792122304, -0.110692433408, 0.671756096776])
       self.bounds[4] = flex.double([-3.13981698125, 2.01054589994, -0.0893939874168, 0.471697111647])
       self.bounds[5] = flex.double([-1.28037493374, 0.876578929774, -0.0363599791812, 0.195285480471])
    if(self.n==7): #coef 7
       self.bounds[0] = flex.double([-20.3241896393, 7.50025274559, -0.0925684548913, 1.90794419639])
       self.bounds[1] = flex.double([-19.271644653, 7.20831756624, -0.0935785997935, 1.87184076462])
       self.bounds[2] = flex.double([-15.8255430056, 5.98498321626, -0.0770293463687, 1.52852426121])
       self.bounds[3] = flex.double([-11.5624283884, 4.56748518877, -0.0646443019148, 1.20681535555])
       self.bounds[4] = flex.double([-7.02554681447, 2.89348936331, -0.0418172018712, 0.751875159574])
       self.bounds[5] = flex.double([-3.64460485882, 1.63779480883, -0.0272670458923, 0.462530654903])
       self.bounds[6] = flex.double([-1.53005568061, 0.649058372991, -0.0112777025715, 0.182214074339])
    if(self.n==8): #coef 8
       self.bounds[0] = flex.double([-15.5696141957, 26.7417847653, -0.429101838698, 2.58512688665])
       self.bounds[1] = flex.double([-14.6649384853, 25.1585349149, -0.379700422122, 2.40814186005])
       self.bounds[2] = flex.double([-12.8078616893, 21.2070970566, -0.380471864759, 2.14461149339])
       self.bounds[3] = flex.double([-9.80561153093, 15.6616904481, -0.291756240601, 1.63458073134])
       self.bounds[4] = flex.double([-6.89022655824, 10.1990490377, -0.254789020103, 1.19989926802])
       self.bounds[5] = flex.double([-4.04173244902, 5.58238627602, -0.156095989293, 0.710323105559])
       self.bounds[6] = flex.double([-2.11816666574, 2.5721593076, -0.112548787807, 0.41572981972])
       self.bounds[7] = flex.double([-0.826062377037, 0.8598617704, -0.0440222198406, 0.158790415553])

  def read_bounds(self,filename):
    file = open(filename,'r')
    start = -1
    for line in file:
      keys = line.split()
      if(start >-1 and start < self.n):
	for ii in range(4):
          self.bounds[start].append(float(keys[ii]))
        start +=1
      if(start == self.n):
        return True

      if(keys[0] == "coef" or keys[0] == "Coef"):
        if(int(keys[1]) == self.n):
	  start = 0
    if(start < 0):
      return False

  def print_coef(self):
    print "Number of coefs = ", self.n
    for ii in range(self.n):
      for jj in range(3):
        print self.bounds[ii][jj],
      print

def test(n):
   filename="test.dat"
   b = bounds(n,filename)
   b.print_coef()

if __name__ == "__main__":
   test(3)
