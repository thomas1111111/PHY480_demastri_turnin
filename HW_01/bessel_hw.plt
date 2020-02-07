# gnuplot plot file: quad_eq.plt

set term pdf
set output "bessel_hw.pdf"

set title 'Tenth Order Bessel Func Error'
set xlabel 'log10(x)'
set ylabel 'log10(Error(x))'
set xrange [-1:2]
set yrange [-17:0.2]
set pointsize 0.2   # set the size of the plotted points
set key bottom left    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "bessel_hw.dat" using 5:6 title 'error'
set out "bessel_hw.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
