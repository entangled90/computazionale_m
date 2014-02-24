#ifndef RACCOLTA_DATI_C

#define RACCOLTA_DATI_C

#include "raccolta_dati.h"



double meanOfDoubleArray( double *array , int n){
	int i = 0;
	double sum = 0;
	for (i= 0; i< n; i++){
		sum += array[i];
	}
	sum /= (double ) n ;
	return ( sum );
	}

	double varianceOfDoubleArray( double *array, int n){
	int i = 0;
	double mean = 0;
	double variance =0.0;
	for (i= 0; i< n; i++){
		mean += array[i];
	}
	mean /= (double)(n);
	for ( i =0;i<n;i++){
	  variance += array[i]*array[i];
	}
	variance /=(double)n;
	variance -= (mean*mean);
	return(variance);
	}

// new_vec deve esser già allocato
void binning ( double * old_vec, double * new_vec, int sizeold, int m){
	int i,j ;
	for ( j = 0 ; j< sizeold/m;j++){
			new_vec[j] = 0;
		}
	for (i = 0; i < sizeold/m ; i++ ){
		for ( j = 0 ; j< m;j++){
			new_vec[i] += old_vec[i*m+j]/((double)m);
		}
	}
}



#endif
