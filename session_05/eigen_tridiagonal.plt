#  file: derivative_test.plt 
#
#  Gnuplot plot file for derivative_test output
#  
#  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
# 
#  Revision history
#   2004-01-24  original version for 780.20 session 5
#   2004-01-16  added postscript enhanced and comments for session 4
#
set term png
set output "eigen_tridiagonal.png"

# record the time and date the graph was generated
set timestamp

# titles and labels
set title 'Prob vs. Dist from Origin (E=1.5)'
set xlabel 'r'
set ylabel 'u(r)'

# move the legend to a free space
set key left

# set the x and y axis scales (already logs)
set xrange [0:5]
set yrange [0:0.25]

# plot the data as well as the fit, with appropriate titles 
plot "eigen_tridiagonal.dat" using ($1):($2) title 'E=1.5'

# output the plot to the file derivative_test_plt.ps   
set out "eigen_tridiagonal.ps"
set term postscript color enhanced
replot
