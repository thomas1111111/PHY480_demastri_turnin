set term png
set output "integ.png"

set xlabel "log10(Number of Steps)"
set ylabel "log10(Relative Error of Integration Routine)"
set title "Analysis of Error Accumulation in Integration Routines"

f(x) = a*x + b
g(x) = c*x + d

fit [3:6] f(x) 'integ.dat' using 1:2 via a,b
fit [3:5] g(x) 'integ.dat' using 1:3 via c,d

set yrange[-16:2]

plot 'integ.dat' using 1:2 title 'Simps Rule',\
    a*x+b title sprintf("%f*x+%f",a,b),\
    'integ.dat' using 1:3 title 'Milne Rule',\
     c*x+d title sprintf("%f*x+%f",c,d)

