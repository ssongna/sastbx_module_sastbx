#!/bin/csh -f

cat << EOF > plot.vmd
#!/usr/local/bin/vmd
# VMD script written by Script
display projection   Orthographic
light 2 on
light 3 on
EOF

set c=0

foreach f ($*)
cat << EOE >> plot.vmd
mol new $f
mol delrep 0 top
mol representation Isosurface 0.15000 0.000000 0.000000 1.000000 1 1
mol color ColorID $c
mol selection {all}
mol material Opaque
mol addrep top
mol representation Cartoon
mol color ColorID $c
mol selection {all}
mol material Opaque
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0
mol scaleminmax top 0 0.000000 0.000000
mol smoothrep top 0 0
mol drawframes top 0 {now}
# done with molecule 0
EOE

@ c++

end

vmd -e plot.vmd
