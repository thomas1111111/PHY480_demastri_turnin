set term png

set output "uniform_1_histogram.png"

set xlabel "Bins from 0 to 1"
set ylabel "Number of Hits"
set title "Uniform Histogram 1"

set timestamp

set style data histogram
set style histogram rowstacked
set yrange [0:2500]
set style fill solid
set style fill solid border -1

set xtics autofreq
set xtics rotate by -90
plot "random_histogram.dat" using 3:xtic(2) title "Uniform 1"


set output "uniform_2_histogram.png"

set xlabel "Bins from 0 to 1"
set ylabel "Number of Hits"
set title "Uniform Histogram 2"

set timestamp

set style data histogram
set style histogram rowstacked
set yrange [0:2500]
set style fill solid
set style fill solid border -1

set xtics autofreq
set xtics rotate by -90
plot "random_histogram.dat" using 4:xtic(2) title "Uniform 2"

set output "gaussian_1_histogram.png"

set xlabel "Bins"
set ylabel "Number of Hits"
set title "Gaussian Histogram 1"

set timestamp

set style data histogram
set style histogram rowstacked
set yrange [0:5000]
set style fill solid
set style fill solid border -1

set xtics autofreq
set xtics rotate by -90
plot "random_histogram.dat" using 6:xtic(5) title "Gaussian 1"


set output "gaussian_2_histogram.png"

set xlabel "Bins"
set ylabel "Number of Hits"
set title "Gaussian Histogram 2"

set timestamp

set style data histogram
set style histogram rowstacked
set yrange [0:5000]
set style fill solid
set style fill solid border -1

set xtics autofreq
set xtics rotate by -90
plot "random_histogram.dat" using 7:xtic(5) title "Gaussian 2"

