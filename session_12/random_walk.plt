set term png
set output "random_walk.png"

set xlabel "x"
set ylabel "y"

set title "Random Walk in Two Dimensions"

plot "random_walk.dat" using 1:2 with lines title "random walk"



