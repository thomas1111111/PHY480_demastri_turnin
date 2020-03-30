# plot file for diffeq_oscillations
set term png
set output 'diffeq_osc2.png'
set timestamp

set title 'Oscillations'

# phase space plot
set xlabel 'x(t)'
set ylabel 'v(t)'
plot "diffeq_oscillations.dat" using ($2):($3) title 'phase-space plot' with lines

