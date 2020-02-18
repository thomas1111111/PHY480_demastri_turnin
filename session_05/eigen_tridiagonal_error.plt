set term png
set output "eigen_error_r4.png"

# record the time and date the graph was generated
set timestamp

# titles and labels
set title 'Error vs. Dimension of Matrix'
set xlabel 'log10(N)'
set ylabel 'log10(error)'

# move the legend to a free space
set key top right

# set the x and y axis scales (already logs)
set xrange [0.29:3.1]
set yrange [-4.2:-.1]

# fit the curve
f1(x) = a1*x + b1
fit [.3:2.0] f1(x) "eigen_tridiagonal.dat" using ($3):($4) via a1,b1
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)



# plot the data as well as the fit, with appropriate titles
plot "eigen_tridiagonal.dat" using ($3):($4) title 'Error', \
     a1*x + b1 title fit_title1

# output the plot to the file derivative_test_plt.ps
set out "eigen_tridiagonal.ps"
set term postscript color enhanced
replot

