from scitbx.array_family import flex
import sys

def read_pr_from_gnom(filename):
  r= flex.double()
  pr= flex.double()
  error= flex.double()

  file=open(filename, 'r')
  pr_start = False
  for line in file:
    keys = line.split('\n')[0].split()
    if(len(keys) == 3):
      if(not pr_start):
        if(keys[0] == 'R' and (keys[1] == 'P(R)') and (keys[2] == 'ERROR') ):
          pr_start=True
      else:
        r.append( float(keys[0]) )
        pr.append( float(keys[1]) )
        error.append( float(keys[2]) )

  file.close()
  return r, pr, error

def print_pr(r,pr):
  for ri,pri in zip(r,pr):
    print ri, pri

def test(filename):
  r,pr,err = read_pr_from_gnom(filename)
  print_pr(r,pr)


if __name__ == "__main__":
  test(sys.argv[1])
