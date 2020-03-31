set xlabel "x(t)"
set ylabel "v(t)"

set title "Phase Diagrams for Different Initial Conditions, Van der Pol"

plot "ode_test_x0_0.1_v0_0.dat" using 2:3 title "x0=0.1,v0=0",\
     "ode_test_x0_1_v0_0.dat" using 2:3 title "x0=1.0,v0=0",\
     "ode_test_x0_-1.5_v0_2.dat" using 2:3 title "x0=-1.5,v0=2.0"
