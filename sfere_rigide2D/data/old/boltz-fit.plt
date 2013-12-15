set autoscale
unset log
set xtic auto
#define fit function
f(x)= N*x*exp(-x*x/2/E)
set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 14 
set output "boltz-fit.png"
N=10^3
E=0.1
fit f(x) 'boltzmann-histo.dat' via N,E
plot f(x), 'boltzmann-histo.dat'

