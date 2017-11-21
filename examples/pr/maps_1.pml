set bg_rgb,[1,1,1]
load m1_3FCN.ccp4, map1
isomesh m1,map1,0.107758460772
color blue, m1
png maps_1_view1
turn x, 90
png maps_1_view2
turn x,-90
turn y,90
png maps_1_view3
