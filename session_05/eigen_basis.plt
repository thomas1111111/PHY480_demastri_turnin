set term png
set output "eigen_basis_error_dim.png"

# record the time and date the graph was generated
set timestamp

# titles and labels
set title 'Error vs. Dimension, b_ho = 0.1585' 
set xlabel 'log10(Dimension)'
set ylabel 'log10(error)'

# move the legend to a free space
set key top right

# set the x and y axis scales (already logs)
set xrange [-0.1:1.8]
set yrange [-7.2:0]


# plot the data as well as the fit, with appropriate titles
plot "eigen_basis.dat" using ($2):($3) title 'Error', 



# output the plot to the file derivative_test_plt.ps
set out "eigen_basis.ps"
set term postscript color enhanced
replot

