set terminal png 
set output "`echo $averages_filename`.png"

set size 1,1
set origin 0,0
set xlabel "Cycle"
set ylabel "Intensity"
plot "`echo $averages_filename`a" using 1:2 with lines title "A Average Intensity" ,\
     "`echo $averages_filename`t" using 1:2 with lines title "T Average Intensity" ,\
     "`echo $averages_filename`g" using 1:2 with lines title "G Average Intensity" ,\
     "`echo $averages_filename`c" using 1:2 with lines title "C Average Intensity"

