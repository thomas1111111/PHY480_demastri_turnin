set term png
set output "derivative_error.png"

set xlabel "log10(h)"
set ylabel "log10(Rel. Error)"

set yrange[-16:0]

set title "Relative Error of Different Differentiation Algorithms"

f1(x) = a1*x + b1
fit [-8:-1] f1(x) "derivative_test.dat" using ($1):($2) via a1,b1
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

f2(x) = a2*x + b2
fit [-5:-1] f2(x) "derivative_test.dat" using ($1):($3) via a2,b2
fit_title2 = sprintf("%-+4.1f*x %-+4.1f",a2,b2)

f3(x) = a3*x + b3
fit [-2.4:-1] f3(x) "derivative_test.dat" using ($1):($4) via a3,b3
fit_title3 = sprintf("%-+4.1f*x %-+4.1f",a3,b3)

f4(x) = a4*x + b4
fit [-1:0.2] f4(x) "derivative_test.dat" using ($1):($5) via a4,b4
fit_title4 = sprintf("%-+4.1f*x %-+4.1f",a4,b4)

plot "derivative_test.dat" using ($1):($2) title 'forward difference', \
     a1*x + b1 title fit_title1, \
     "derivative_test.dat" using ($1):($3) title 'central difference', \
     a2*x + b2 title fit_title2,\
     "derivative_test.dat" using ($1):($4) title 'extrapolation diff',\
     a3*x + b3 title fit_title3,\
     "derivative_test.dat" using ($1):($5) title 'extrapolation diff two',\
     a4*x + b4 title fit_title4
    

