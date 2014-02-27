#ifndef RACCOLTA_DATI_H

#define RACCOLTA_DATI_H
void vec_zeros(double * v, int N);
void binning ( double * ,double *, int ,int );
double meanOfDoubleArray( double *array , int n);
double varianceOfDoubleArray( double *array, int n);
void autocorrelation (double * dati, double * risultato, int N,int corr_max);
void divideByScalar(double * v, double scalar, int N);
void binning_mat ( double * old_vec, double * new_vec,int n_row ,int n_data, int bin_width);
#endif