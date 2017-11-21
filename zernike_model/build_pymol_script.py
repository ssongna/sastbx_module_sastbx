import sys,os

def run(args):
  maps = []
  levels = []
  pdbfile = None
  for ff in args:
    if ('xplor' in ff) or ('ccp4' in ff):
      maps.append( ff )
      levels.append(0.08)
    if 'pdb' in ff:
      pdbfile = ff

  write_pymol_scripts(maps,levels, pdbfile)

def write_pymol_scripts(maps,levels,pdbfile=None,root_name='maps'):
  for ii,map in enumerate(maps):
    ff = open( '%s_%s.pml'%(root_name,ii+1), 'w' )
    print >> ff, 'set bg_rgb,[1,1,1]'
    if pdbfile is not None:
      print >> ff, 'load %s'%pdbfile
      print >> ff, 'hide everything,all'
      print >> ff, 'show cartoon, all'
    print >> ff, 'load %s, map1'%(map)
    print >> ff, 'isomesh m1,map1,%s'%levels[ii]
    print >> ff, 'color blue, m1'
    print >> ff, 'png %s_%s_view1'%(root_name,ii+1)
    print >> ff, 'turn x, 90'
    print >> ff, 'png %s_%s_view2'%(root_name,ii+1)
    print >> ff, 'turn x,-90'
    print >> ff, 'turn y,90'
    print >> ff, 'png %s_%s_view3'%(root_name,ii+1)
    ff.close()


def write_pymol_shapeup(maps,pdbfile=None,root_name='maps'):
  color_list = ["brown","carbon","cyan","deepblue","deepsalmon",\
  "deepteal","gray","lightpink","olive","palecyan"]

  targetfile = root_name+".pml"
  with open(targetfile,"w") as ff:
    print >> ff, 'set bg_rgb,[1,1,1]'
    for ii,xplor in enumerate(maps):
      tmp_name = xplor.split("/")[-1].split("_")[-1].split(".")[0]
      map_name = "map_%s_%s" %(ii+1,tmp_name)
      m_name = "m_%s_%s" %(ii+1,tmp_name)
      if pdbfile is not None:
        print >> ff, 'load %s'%pdbfile
        print >> ff, 'hide everything,all'
        print >> ff, 'show cartoon, all'
      print >> ff, 'load %s, %s'%(xplor,map_name)
      print >> ff, 'isosurface %s,%s,%s' %(m_name,map_name,str(1.0))
      print >> ff, 'color %s, %s' %(color_list[ii],m_name)
       
def write_pymol_superpose(pdblist,targetdir=None):
  color_list = ["brown","carbon","cyan","deepblue","deepsalmon",\
  "deepteal","gray","lightpink","olive","palecyan"]
  targetfile = os.path.join(targetdir,"superpose.pml")
  with open(targetfile,"w") as pml:
    for ii, pdb in enumerate(pdblist):
      pdbname = pdb.split("/")[-1].split(".")[0]
      print >> pml, "load %s,  pdb_%s" %(pdb, pdbname)
      print >> pml, "show_as surface"   
      print >> pml, "set transparency, 0.5"
      print >> pml, "color %s, pdb_%s" %(color_list[ii],pdbname)
    
    print >> pml, "zoom"
    print >> pml, "orient"    


if __name__ == "__main__":
  run(sys.argv[1:])
