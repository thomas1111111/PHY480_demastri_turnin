# plot file for diffeq_test_exp_back
set term png

set output 'diffeq_test_new.png'

set timestamp

set key top right

set xlabel 'log10(h)'
set ylabel 'log10(error)'

set xrange [-14.1:-0.9]
set yrange [-18:0]

set title 'Log(Rel Error) of Differential Equation Algorithms at y=1.'

plot \
  "diffeq_test.dat" using (log10(($1))):(log10(abs(((($2)-($4))/($4))))) title 'Euler Error', \
  "diffeq_test.dat" using (log10(($1))):(log10(abs(((($3)-($4))/($4))))) title 'RK Error',

