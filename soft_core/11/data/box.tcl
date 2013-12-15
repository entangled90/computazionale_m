set minx 0
set miny 0
set minz 0
set maxx 6.576566e+00
set maxy 6.576566e+00
set maxz 6.576566e+00
draw materials off 
 draw color yellow 
 draw line " $minx $miny $minz " " $maxx $miny $minz"
 draw line "$minx $miny $minz" "$minx $maxy $minz"
 draw line "$minx $miny $minz" "$minx $miny $maxz"
 draw line "$maxx $miny $minz" "$maxx $maxy $minz"
 draw line "$maxx $miny $minz" "$maxx $miny $maxz"
 draw line "$minx $maxy $minz" "$maxx $maxy $minz"
draw line "$minx $maxy $minz" "$minx $maxy $maxz"
 draw line "$minx $miny $maxz" "$maxx $miny $maxz"
 draw line "$minx $miny $maxz" "$minx $maxy $maxz"
 draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
 draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
 draw line "$maxx $maxy $maxz" "$maxx $miny $maxz" 
 mol default style VDW
mol addfile vmd.xyz
set sel [atomselect top all]
$sel set radius 0.1
