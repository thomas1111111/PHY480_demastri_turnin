# plot file for diffeq_oscillations
set term png
set output 'diffeq_osc5.png'

set timestamp

set title 'Oscillations'

set yrange [0:0.6]


# plot of kinetic energy and potential energy versus time
set xlabel 't'
set ylabel 'energy'
plot "diffeq_oscillations.dat" using ($1):(($4)+($5)) title 'Total Energy' with lines
 
