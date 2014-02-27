#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "mtwist.h"
#include "constants-metro.h"
#include "metropolis.h"
#include "raccolta_dati.h"

#define CORR_MAX 50

int main (int argc, char *argv[]){	
	mt_seed();
	double BETA = 1;
	int N ;
	int iteration = 0;
	double mag2_mean=0; // Valor medio di magnetizzazione al quadrato
	double mag_mean = 0; // valor medio di modulo magnetizzazione
	double chi = 0;
	int i,j ;
	if (argc>1){
		BETA = atof(argv[1]);
		N = atoi(argv[2]);
	}
	else{
		printf("Inserire il valore di Beta\n");
		return 0;
	}
	double * X_n= malloc(sizeof(double)*N);
	double * Y_n=malloc(sizeof(double)*N);

	// init vettore
	for ( i = 0; i<N;i++){
		X_n[i]=0;
		Y_n[i]=0;
	}

//	double mag_prevista = ( 1 - pow(sinh(2*BETA*J),-4));
	short int  * configuration;
	srand(time(NULL));
	configuration = malloc(N*N*sizeof(short int));
	spin_init(configuration,N);
	/***** FILENAMES AND FILE OPENING ******/
/* ------------------MAGN*/
	char mag_filename[64] = "";
	snprintf(mag_filename,64,"data/mag_mean%d.dat",N);
	char chi_filename[64] = "";
	snprintf(chi_filename,64,"data/chi%d.dat",N); 
	char mag_binning_filename[64] = "";
	snprintf(mag_binning_filename,64,"data/binning/mag_N%d__B%.4lf.dat",N,BETA); 
	char mag_autocorr_filename[64] = "";
	snprintf(mag_autocorr_filename,64,"data/mag_corr/mag_autocorrN%d_B%.4lf.dat",N,BETA);
/*----------------------_ENergia*/
	char en_binning_filename[64] = "";
	snprintf(en_binning_filename,64,"data/binning/en_N%d__B%.4lf.dat",N,BETA);
	char en_filename[64] = "";
	snprintf(en_filename,64,"data/en_N%d.dat",N);
	char en_autocorr_filename[64] = "";
	snprintf(en_autocorr_filename,64,"data/en_corr/en_autocorrN%d_B%.4lf.dat",N,BETA);
	char en_temp_filename[64] = "data/en_temp.dat";
	char cv_filename[64]="";
	snprintf(cv_filename,64,"data/cv%d.dat",N);

	char corr_row_filename[64]="";
	snprintf(corr_row_filename,64,"data/corr_row/corr_row_N%dB%.4lf.dat",N,BETA);
/*** I File sono chiamati: f_$(Nomestringa) ****/
	FILE * f_mag = fopen(mag_filename,"a");
	FILE * f_chi = fopen(chi_filename,"a");	
	FILE * f_mag_bin = fopen(mag_binning_filename,"w");
	FILE * f_mag_autocorr = fopen(mag_autocorr_filename,"w");
	FILE * f_mag_temp = fopen("data/mag_temp.dat","w");

	FILE * f_en_bin = fopen(en_binning_filename,"w");
	FILE * f_en_autocorr = fopen(en_autocorr_filename,"w");
	FILE * f_en = fopen(en_filename,"a");
	FILE * f_en_temp = fopen(en_temp_filename,"w");
	FILE * f_cv = fopen(cv_filename,"a");

	FILE * f_corr_row = fopen(corr_row_filename,"w");



/**********************************************************************
 TERMALIZZAZIONE
 **********************************************************************/

	while(iteration < ITERATION_THERM){
		metropolis_ising(configuration,N,BETA);
		iteration++;
	}

	iteration = 0;	
	double mag_tmp;
	double en_tmp;
	double en_mean =0;
	double en2_mean=0;
	double cv=0;
	double * mag_vet_dati ;
	double * mag_vet_binnato;
	double * en_vet_dati;
	double * en_vet_binnato;
	double * en_autocorr = malloc(sizeof(double)*CORR_MAX);
	double * mag_autocorr = malloc(sizeof(double)*CORR_MAX);
	double * S_xt = malloc(sizeof(double)*N);
	double * S_yt = malloc(sizeof(double)*N);
	en_vet_dati = malloc(sizeof(double)*ITERATION_MAX);
	mag_vet_dati = malloc(sizeof(double)*ITERATION_MAX);
	int larghezza_bin;

	while(iteration < ITERATION_MAX){
		metropolis_ising(configuration,N,BETA);
		mag_tmp = magnetization(configuration,N);
		en_tmp = hamiltoniana(configuration,N);
		mag_vet_dati[iteration] = fabs(mag_tmp);
		en_vet_dati[iteration] = en_tmp;
		mag2_mean += mag_tmp*mag_tmp;
		mag_mean += fabs(mag_tmp);
		en_mean+=en_tmp;
		en2_mean += en_tmp*en_tmp;
		iteration++;
		for ( i = 0; i<N;i++){
			X_n[i] = sum_row(configuration,i,N);
			Y_n[i] = sum_col(configuration,i,N);
		}
		for (j = 0; j < N; ++j){
			for (i = 0; i < N;i++){
				S_xt[i]+= X_n[j]*X_n[(i+j)%N];
				S_yt[i]+= Y_n[j]*Y_n[(i+j)%N];
			}
		}
	}
	for (i = 0; i < N;i++){
		S_xt[i]/=(double)N;
		S_yt[i]/=(double)N;
	}

	for ( i = 0; i<N;i++){
		S_xt[i] += S_yt[i];
		S_xt[i] += S_xt[N-1-i];
		S_xt[i] += S_yt[N-1-i];
		S_xt[i]/=3.0*ITERATION_MAX;
	}
	for ( i = 0; i<N/3;i++){
		fprintf(f_corr_row, "%d\t%lf\n",i,S_xt[i] );
	}

	mag_mean /= (double)(ITERATION_MAX);
	mag2_mean /= (double)(ITERATION_MAX);
	en_mean /= (double)(ITERATION_MAX);
	en2_mean /= (double)(ITERATION_MAX);
	//mag_mean /= (double)(ITERATION_MAX);
	chi = (mag2_mean - mag_mean*mag_mean)/(double)(N*N);
	cv = (en2_mean - en_mean*en_mean)/(double)(N*N);
	mag_mean /= (double)(N*N);
	en_mean /= (double)(N*N);
/* Scrivo su file i dati calcolati*/
	fprintf(f_mag,"%lf\t%lf\n",BETA,mag_mean);
	fprintf(f_chi,"%lf\t%lf\n",BETA,chi);
	fprintf(f_en,"%.14e\t%.14e\n",BETA,en_mean);
	fprintf(f_cv,"%lf\t%lf\n",BETA,cv);

/**** BINNING E AUTOCORRELAZIONI */
	divideByScalar(mag_vet_dati,N*N,ITERATION_MAX);
	divideByScalar(en_vet_dati,N*N,ITERATION_MAX);
	mag_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
	en_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
	for ( larghezza_bin = 1; larghezza_bin < 100 ; larghezza_bin+=1){
		binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
		binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
		fprintf(f_mag_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(mag_vet_binnato,ITERATION_MAX/larghezza_bin)));
		fprintf(f_en_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(en_vet_binnato,ITERATION_MAX/larghezza_bin)));
	}

	/** Autocorrelazioni **/
	autocorrelation(en_vet_dati,en_autocorr,ITERATION_MAX,CORR_MAX);
	if(f_en_autocorr){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(f_en_autocorr,"%d\t%.14e\n",i,en_autocorr[i]);
		}
	}
	autocorrelation(mag_vet_dati,mag_autocorr,ITERATION_MAX,CORR_MAX);
	if(f_mag_autocorr){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(f_mag_autocorr,"%d\t%.14e\n",i,mag_autocorr[i]);
		}
	}



	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_en_temp,"%.14e\n",en_vet_dati[i]);
	}
	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_mag_temp,"%.14e\n",mag_vet_dati[i]);
	}
/* Chiusura file */
	fclose(f_mag_bin);
	fclose(f_mag);
	fclose(f_mag_temp);
	fclose(f_mag_autocorr);
	fclose(f_chi);

	fclose(f_en);
	fclose(f_en_temp);
	fclose(f_en_bin);
	fclose(f_en_autocorr);
	fclose(f_cv);

	fclose(f_corr_row);

		/* Free della memoria */
	free(mag_vet_dati) ;
	free(mag_vet_binnato);
	free(en_vet_dati);
	free(en_vet_binnato);
	free(en_autocorr);
	free(mag_autocorr);
	free(configuration);
	//printf("Magnetizzazione media: %lf\n",sum/((double) ITERATION_MAX));
	exit(EXIT_SUCCESS);
}
