#ifndef RACCOLTA_DATI_C

#define RACCOLTA_DATI_C

#include "raccolta_dati.h"
#include <stdio.h>
#include <stdlib.h>

void vec_zeros(double * v, int N){
	int i;
	for (i=0;i<N;i++){
		v[i]=0;
	}
}


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


/* Binning per grandezze derivate dalle osservabili, come il calore specifico e chi.
Esse vengono calcolate in ogni bin 
*/
// new_vec deve esser già allocato
void binning_deriv ( double * old_vec, double * new_vec, int sizeold, int m){
	int i,j ;
	double * tmp = malloc(sizeof(double)*(sizeold/m));
	for ( j = 0 ; j< sizeold/m;j++){
			new_vec[j] = 0;
			tmp[j]=0;
		}
		printf("sizeold/m = %d\n",sizeold/m);
	for (i = 0; i < sizeold/m ; i++ ){
		for ( j = 0 ; j< m;j++){
			new_vec[i] += old_vec[i*m+j]*old_vec[i*m+j];
			tmp[i] += old_vec[i*m+j];
		}
	tmp[i] /= (double)(m);
	new_vec[i] /= (double)(m);
	new_vec[i] -= tmp[i]*tmp[i];
	}

}



/** Fa il binning di dati vettoriali, come la correlazione S_t
n_vec è la lunghezza del vettor
n_data è il numero di iterazioni -> insieme di vettori indipendenti
bin_width è la larghezza dela bin
*/
void binning_mat ( double * old_vec, double * new_vec, int vec_len ,int n_data, int bin_width){
	int i,j ;
	for ( j = 0 ; j< (n_data/bin_width)*vec_len;j++){
		new_vec[j] = 0;
	}
	for(i = 0; i<n_data - n_data%bin_width;i++){
		for(j = 0; j<vec_len;j++){
			new_vec[(i%bin_width)*vec_len+j] += old_vec[i*(vec_len)+j];
		}
	}
	for(i = 0;i< n_data/bin_width;i++){
		for(j=0;j<vec_len;j++){
			new_vec[i*vec_len+j]/=(double)bin_width;
		}
	}
}




void autocorrelation (double * dati, double * risultato, int N, int corr_max){
	int i,j;
	double mean = meanOfDoubleArray(dati,N);
	double var = varianceOfDoubleArray(dati,N);
	printf("Media: %lf \t, Varianza: %lf\n",mean,var);
	for (j = 0; j < corr_max;j++){
		risultato[j] = 0;
	}
	for ( j= 0; j < corr_max;j++){
		for (i =0; i< N-j ;i++){
			risultato[j] += dati[i]*dati[i+j] / (double)(N-j);
		}
	}
	for (j = 0; j<corr_max;j++){
			risultato[j] -= mean*mean;
			risultato[j] /= var;
	}
}
void divideByScalar(double * v, double scalar, int N){
	int i;
	for ( i = 0; i<N;i++){
		v[i] /= scalar;
	}
}

void multByScalar(double * v, double scalar, int N){
	int i;
	for ( i = 0; i<N;i++){
		v[i] *= scalar;
	}
}


#endif
