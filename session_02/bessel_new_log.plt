# gnuplot plot file: quad_eq.plt

set term png
set output "bessel_error_log.png"

set logscale
set title 'Bessel_Func Relative Error'
set xlabel 'log(x)'
set ylabel 'log(error(x))'
set xrange [0.1:100]
set yrange [0.0000000000000001:10]
set pointsize 1   # set the size of the plotted points
set key bottom left    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "bessel_new.dat" using 1:4 title 'relative difference'
set out "bessel.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
