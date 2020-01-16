# gnuplot plot file: quad_eq.plt

set term png
set output "bessel.png"

set title 'Tenth Order Bessel Func'
set xlabel 'x'
set ylabel 'J10(x)'
set xrange [1:100]
set yrange [0:0.15]
set pointsize 0.2   # set the size of the plotted points
set key top left    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "bessel.dat" using 1:3 title 'upward recursion',\
     "bessel.dat" using 1:2 title 'downward recursion'
set out "bessel.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
