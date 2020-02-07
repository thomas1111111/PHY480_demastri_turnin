#  file: sum.plt
#
#  Based off "derivative_test.plt", which was made by Dick Furnstahl, furnstahl.1@osu.edu
# 
#  Created 2/6/20, 8:04 P.M.
#
#  Revision History:
#  -- edited the derivative_test.plt file to plot my sum_up_down data instead
#  -- have it recording a pdf now instead of a png
#
set term pdf
set output "sum_error.pdf"

# record the time and date the graph was generated
set timestamp

# titles and labels
set title 'Comparison of Sum Up and Sum Down Method'
set xlabel 'log10(N)'
set ylabel 'log10(diff/average)'

# move the legend to a free space
set key left

# set the x and y axis scales (already logs)
set xrange [0.9:6]
set yrange [-1.7:0]

# fit the curve
f1(x) = a1*x + b1
fit [2.5:6] f1(x) "sum_up_down.dat" using ($2):($3) via a1,b1 
fit_title1 = sprintf("%-+4.1f*x %-+4.1f",a1,b1)

# plot the data as well as the fit, with appropriate titles 
plot "sum_up_down.dat" using ($2):($3) title 'diff/average', \
     a1*x + b1 title fit_title1

# output the plot to the file derivative_test_plt.ps   
set out "sum_up_down.ps"
set term postscript color enhanced
replot
