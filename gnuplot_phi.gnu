#!/usr/bin/gnuplot
set terminal pngcairo size 600, 600 enhanced font 'Verdana, 10' 
set output 'phi2d.png'
#set grid

#GRAFICO DO MEIO

set lmargin at screen 0.1
set rmargin at screen 0.80
set bmargin at screen 0.1
set tmargin at screen 0.9

unset key
set pm3d
unset grid
set size ratio -1
set palette defined ( 0 "black", 0.25 "red", 0.5 "orange", 0.75 "yellow",  1.0 "white")


set view map
set pm3d noborder 
set hidden3d
unset surface
set ytics 0.1 
unset xtics 
set xrange[-1.4 : 1.4]
set yrange[-1.4 : 1.4]
set ylabel "y(nm)" font "Verdana, 18"
#set xlabel "x(nm)" font "Verdana, 18" 
set cblabel "4 {/Symbol p} r^2 |{/Symbol F}(r, l)|^2" font "Verdana, 18"
set y2range[0: 2]
set cbrange[0: 2]


splot "phipolar" u 1:2:3 w pm3d 



