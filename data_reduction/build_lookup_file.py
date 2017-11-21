import sys,os

def print_help():
  print "Suggested usuage:"
  print "  sastbx.build_lookup_table <directory> <norma extension> <image extension>  >  aux_info  "
  print "  Very 7.3.3. pilatus specific."
  print

def parse_info_file(file_name):
  ff = open(file_name,'r')
  i0_i1 = []
  for ll in ff:
    try:
      i0_i1.append( float(ll) )
    except: pass
  return i0_i1

def process(directory, info_ext, img_ext, img_dir=None ):
  list_dir = os.listdir( directory )
  if img_dir is None:
    img_dir = directory
  for ff in list_dir:
    n = len(ff)
    m2 = len(info_ext)
    ext2 = ff[ n-m2: ]
    if info_ext == ext2:
      local_img_name = ff[ :n-m2 ]+img_ext
      full_path_name = os.path.join( img_dir, local_img_name )
      i0i1 = parse_info_file( os.path.join(directory,ff) )
      if len(i0i1)==2:
        print full_path_name, i0i1[0], i0i1[1]

if __name__ == "__main__":
  if len(sys.argv)!=4:
    print_help()
  else:
    process( sys.argv[1], sys.argv[2], sys.argv[3] )
