#!/bin/csh
cat start.pdb > combined.pdb
echo 'ENDMOL' >> combined.pdb
set indx = 1
while(-e $1$indx.pdb)
cat ${1}${indx}.pdb >> combined.pdb
echo 'ENDMOL' >> combined.pdb
@ indx = $indx + 1
end
