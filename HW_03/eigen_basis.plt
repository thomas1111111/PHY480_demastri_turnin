set term png
set output "eigen_basis.png"

set xrange [0:10]
set yrange [0:2]

plot 2*x*exp(-x),\
    "eigen_basis_dimension_1_bho_0.158489.dat" using ($1):(-$2) title "dim=1",\
    "eigen_basis_dimension_5_bho_0.158489.dat" using ($1):(-$2) title "dim=5",\
    "eigen_basis_dimension_10_bho_0.158489.dat" using ($1):(-$2) title "dim=10",\
    "eigen_basis_dimension_20_bho_0.158489.dat" using ($1):(-$2) title "dim=20"

