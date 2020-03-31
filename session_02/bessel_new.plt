# gnuplot plot file: quad_eq.plt
set title 'Bessel_Func Relative Error'
set xlabel 'x'
set ylabel "Error(x)"
set xrange [0:100]
set yrange [0:1.2]
set pointsize 1   # set the size of the plotted points
set key top left    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "bessel_new.dat" using 1:4 title 'relative difference'
set out "bessel.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
