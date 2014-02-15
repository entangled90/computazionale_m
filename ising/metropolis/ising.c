#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "mtwist.h"
#include "constants-metro.h"
#include "metropolis.h"

double BETA = 1;



int main (int argc, char *argv[]){	
	mt_seed();
	int iteration = 0;
	double mag_abs_mean=0; // Valor medio del modulo della magnetizzazione
	double mag2_mean=0; // Valor medio di magnetizzazione al quadrato
	double mag_mean = 0; // valor medio di magnetizzazione
	double chi = 0;
	double tmp= 0;
	double S_n[N];
	int i = 0;
	// init vettore
	for ( i = 0; i<N;i++){
		S_n[i]=0;
	}
	if (argc>1){
		BETA = atof(argv[1]);
	}
	else{
		printf("Inserire il valore di Beta\n");
		return 0;
	}

//	double mag_prevista = ( 1 - pow(sinh(2*BETA*J),-4));
	short int  * configuration;
	srand(time(NULL));
	configuration = malloc(N*N*sizeof(short int));
	spin_init(configuration);
	printf("hamiltoniana %e\n",hamiltonian_ising(configuration));

	/***** FILENAMES AND FILE OPENING ******/
	char mag_file[64] = "";
	snprintf(mag_file,64,"data/magnetization_mean_%d.dat",N);
	char chi_file[64] = "";
	snprintf(chi_file,64,"data/chi_%d.dat",N); 
	char corr_row_file[64]= "";
	snprintf(corr_row_file,64,"data/corr_row_%dParticle_%lf.3BETA.dat",N,BETA);
	FILE * f_mag_mean = fopen(mag_file,"a");
	FILE * f_chi = fopen(chi_file,"a");	
	FILE * f_corr_row = fopen(corr_row_file,"w");
	while(iteration < ITERATION_THERM){
		if(iteration %200 == 0){
			printf("Iterazione: %d\n",iteration);
		}
		metropolis_ising(configuration);
		iteration++;
	}
	iteration = 0;	
	while(iteration < ITERATION_MAX){
		if(iteration %200 == 0){
			printf("Iterazione: %d\n",iteration);
		}
		metropolis_ising(configuration);
		tmp = magnetization(configuration);
		mag2_mean += tmp*tmp;
		mag_abs_mean += fabs(tmp);
		mag_mean  += tmp ;
		iteration++;
		for ( i = 0; i<N;i++){
		S_n[i] += sum_row(configuration,i);
		}
	}

	for ( i = 0; i<N;i++){
		S_n[i]/=(double)ITERATION_MAX;
		fprintf(f_corr_row, "%d\t%lf\n",i,S_n[i] );
	}
	mag_abs_mean /= (double)(ITERATION_MAX);
	mag2_mean /= (double)(ITERATION_MAX);
	mag_mean /= (double)(ITERATION_MAX);
	chi = (mag2_mean - mag_mean*mag_mean);
	fprintf(f_mag_mean,"%lf\t%lf\n",BETA,mag_abs_mean);
	fprintf(f_chi,"%lf\t%lf\n",BETA,chi);
	fclose(f_mag_mean);
	fclose(f_chi);
	fclose(f_corr_row);
	//printf("Magnetizzazione media: %lf\n",sum/((double) ITERATION_MAX));
	exit(EXIT_SUCCESS);
}
