set term png
set output "fit.png"


set xlabel "t"
set ylabel "y(t)"

set title "Testing GSL Fitting to Exponential decay"

plot "multifit_test.dat" using 1:2:4 with errorbars title "Pseudo-Data",\
     "multifit_test.dat" using 1:3 title "Fit"
