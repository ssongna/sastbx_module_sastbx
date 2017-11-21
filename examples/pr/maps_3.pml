set bg_rgb,[1,1,1]
load m3_2WAD.ccp4, map1
isomesh m1,map1,0.0872160806585
color blue, m1
png maps_3_view1
turn x, 90
png maps_3_view2
turn x,-90
turn y,90
png maps_3_view3
