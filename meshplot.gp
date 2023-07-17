reset session
set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb 'white' fillstyle solid noborder

set terminal pdfcairo size 300cm,300cm
set autoscale xfix

set border 0
set xzeroaxis lt -1
set yzeroaxis lt -1

plot 'mesh.dat' u 2:3:4 w points pointtype 7 pointsize 2 lt palette
