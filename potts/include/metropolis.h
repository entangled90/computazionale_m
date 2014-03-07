#ifndef METROPOLIS_H
#define METROPOLIS_H

#ifndef METROPOLIS_C
extern void metropolis_ising( short int *x ,int N, double BETA);
extern double delta_Ham_ising ( short int * configuration,int a, int b, int N);
extern double hamiltoniana(short int * configuration, int N);
void spin_init (short int * configuration, int N);
inline void savePPM(short int * x, int N);
inline double magnetization(short int *x, int N);
double sum_row(short int * configuration, int row, int N);
double sum_col(short int * configuration, int col, int N);
#endif

#endif
