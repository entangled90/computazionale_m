#ifndef METROPOLIS_H
#define METROPOLIS_H

#ifndef METROPOLIS_C
extern void metropolis_ising( short int *x );
extern double delta_Ham_ising ( short int * configuration,int a, int b);
extern double hamiltonian_ising(short int * configuration);
void spin_init (short int * configuration);
inline void savePPM(short int * x);
inline double magnetization(short int *x);
double sum_row(short int * configuration, int row);

#endif

#endif
