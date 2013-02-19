set terminal png 
set output "`echo $intensity_filename`.png"

set size 3,2
set multiplot

set xrange [ -1000 : 16000 ] 
set yrange [ -1000 : 16000 ] 
set xlabel "C Intensity"
set ylabel "A Intensity"
set style data points

set size 1,1
set origin 0,0
set xlabel "C Intensity"
set ylabel "A Intensity"
plot "`echo $intensity_filename`ca" title "CA"

set size 1,1
set origin 1,0
set xlabel "C Intensity"
set ylabel "G Intensity"
plot "`echo $intensity_filename`cg" title "CG"

set size 1,1
set origin 2,0
set xlabel "C Intensity"
set ylabel "T Intensity"
plot "`echo $intensity_filename`ct" title "CT"

set size 1,1
set origin 0,1
set xlabel "G Intensity"
set ylabel "A Intensity"
plot "`echo $intensity_filename`ga" title "GA"

set size 1,1
set origin 1,1
set xlabel "G Intensity"
set ylabel "T Intensity"
plot "`echo $intensity_filename`gt" title "GT"

set size 1,1
set origin 2,1
set xlabel "A Intensity"
set ylabel "T Intensity"
plot "`echo $intensity_filename`at" title "AT"
