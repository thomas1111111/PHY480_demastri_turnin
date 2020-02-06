# gnuplot plot file: quad_eq.plt

set term png
set output "derivative.png"

set title 'Simple Derivative Test'

set xlabel 'log10(h) (mesh size)'
set ylabel 'log10(relative error)'
set xrange [-12:0]
set yrange [-12:0]

set pointsize 0.5   # set the size of the plotted points
set key top left    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "derivative_test_simple.dat" using 1:2 title 'forward difference',\
     "derivative_test_simple.dat" using 1:3 title 'central recursion'
set out "derivative.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
