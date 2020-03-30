# plot file for diffeq_oscillations


# simple plot versus time
set term png   # plot window 1
set output "diffeq_osc1.png"
set timestamp
set title 'Oscillations: Time Dependence'
pi = 3.141592653589793
k = (2.*pi)**2
m = 1
omega = sqrt(k/m)
set xlabel 't'
set ylabel 'y(t)'
plot "diffeq_oscillations.dat" using ($1):($2) title 'p=?' with lines 
#plot "diffeq_oscillations.dat" using ($1):($2) title 'p=?' with lines # , \
#  (1/omega)*sin(omega*x)


# phase space plot
set term png   # plot window 1
set output "diffeq_osc2.png"
set timestamp
set title 'Oscillations: Phase Space'
set xlabel 'x(t)'
set ylabel 'v(t)'
plot "diffeq_oscillations.dat" using ($2):($3) title 'phase-space plot' with lines


# plot the potential
set term png   # plot window 1
set output "diffeq_osc3.png"
set timestamp
set title 'Oscillations: Potential V(x)'
set xlabel 'x'
set ylabel 'V(x)'
plot "diffeq_oscillations.dat" using ($2):($5) title 'V(x)' with lines


# plot of kinetic energy and potential energy versus time
set term png   # plot window 1
set output "diffeq_osc4.png"
set timestamp
set title 'Oscillations: Energy'
set xlabel 't'
set ylabel 'energy'
plot "diffeq_oscillations.dat" using ($1):($4) title 'KE(t)' with lines, \
     "diffeq_oscillations.dat" using ($1):($5) title 'PE(t)' with lines
 
