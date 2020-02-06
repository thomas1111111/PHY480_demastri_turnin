# gnuplot plot file: quad_eq.plt

set term png
set output "integ_log.png"

set title 'Integral Test'

set xlabel 'log10(N)'
set ylabel 'log10(relative error)'
set xrange [0.25:2.75]
set yrange [-10:-1]

f(x) = m*x+b
g(x) = n*x+c
h(x) = o*x+d

fit [0.4:2.75] f(x) 'integ.dat' using 1:2 via m,b
fit [0.4:1.25] g(x) 'integ.dat' using 1:3 via n,c
fit [0.4:0.75] h(x) 'integ.dat' using 1:4 via o,d


set pointsize 0.5   # set the size of the plotted points
set key top right    # move the key away from the lines
set timestamp	  # turn on a date/time indicator
plot "integ.dat" using 1:2 title 'trap rule',\
     "integ.dat" using 1:3 title 'simposons',\
     "integ.dat" using 1:4 title 'gauss',\
     m*x+b title 'fit trap',\
     n*x+c title 'fit simp',\
     o*x+d title 'fit gaus'


set out "derivative.ps"	# name of the output postscript file
set terminal postscript	        # switch to postscript mode
replot  		  	# plot to the file
