set xrange [-1:201]
set yrange [0:90]

set title "GSL Interpolation on Cross Section Data"
set xlabel "Energy (MeV)"
set ylabel "sigma"

plot "spline.dat" using 1:2 with lines title "theoretical",\
     "spline.dat" using 1:3 with linespoints title "cubic spline",\
     "spline.dat" using 1:4 with linespoints title "linear spline",\
     "spline.dat" using 1:5 with linespoints title "poly spline",\

