set bg_rgb,[1,1,1]
load m10_2P1Y.ccp4, map1
isomesh m1,map1,0.0972095968066
color blue, m1
png maps_10_view1
turn x, 90
png maps_10_view2
turn x,-90
turn y,90
png maps_10_view3
