#ifndef METROPOLIS_H
#define METROPOLIS_H

#ifndef METROPOLIS_C
extern void metropolis_ising( short int *x );
extern double delta_Ham_ising ( short int * configuration,int a, int b);
extern double hamiltonian_ising(short int * configuration);
#endif

#endif
