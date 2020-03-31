set xlabel 'Time'
set ylabel 'theta'

set xrange [0:100]
set yrange [-1:1]

set title 'Theta vs. Time'

plot 'diffeq_pendulum0.2.dat' using 1:2


