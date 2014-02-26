#ifndef RACCOLTA_DATI_H

#define RACCOLTA_DATI_H

void binning ( double * ,double *, int ,int );
double meanOfDoubleArray( double *array , int n);
double varianceOfDoubleArray( double *array, int n);
void autocorrelation (double * dati, double * risultato, int N,int corr_max);
void divideByScalar(double * v, double scalar, int N);
#endif