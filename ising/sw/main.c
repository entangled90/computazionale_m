#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// my includes
#include "sw.h"
#include "list.h"
#include "constants-sw.h"
#include "mtwist.h"
#include "raccolta_dati.h"

#define CORR_MAX 30
int main ( int argc, char * argv[]) {
	float BETA = 1;
	int N = 40;
	double mag2_mean=0; // Valor medio di magnetizzazione al quadrato
	double mag_mean = 0; // valor medio di magnetizzazione
	double chi = 0;
	int iteration = 0;
	mt_seed();
	int i,j;
	int tau_corr=20;
	int larghezza_bin = 5*tau_corr;
	int n_bin = ITERATION_MAX/larghezza_bin;
	/*Check for command line arguments*/
	if (argc>1){
		BETA = atof(argv[1]);
		N = atoi(argv[2]);
	}
	else{
		printf("Inserire il valore di Beta e N\n");
		exit(1);
	}
	int N_CORR = N/3;
	double * X_n= malloc(sizeof(double)*N);
	double * Y_n=malloc(sizeof(double)*N);
		// init vettore
	for ( i = 0; i<N;i++){
		X_n[i]=0;
		Y_n[i]=0;
	}
	Spin * matrix = (Spin *) malloc(sizeof(Spin)*N*N);
	Node * nodes= (Node *) malloc(sizeof(Node)*N*N);
	/*Check for allocation*/
	if( !matrix || !nodes){
		printf("Cannot call malloc || MAIN \n");
		exit(1);
	}
//	float index_simulation = mt_drand();
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
//	char en_temp_filename[64] = "data/en_temp.dat";
	char cv_filename[64]="";
	snprintf(cv_filename,64,"data/cv%d.dat",N);

	char corr_row_filename[64]="";
	snprintf(corr_row_filename,64,"data/corr_row/corr_row_N%dB%.4lf.dat",N,BETA);

/*** I File sono chiamati: f_$(Nomestringa) ****/
	FILE * f_mag = fopen(mag_filename,"a");
	FILE * f_chi = fopen(chi_filename,"a");	
	FILE * f_mag_bin = fopen(mag_binning_filename,"w");
	FILE * f_mag_autocorr = fopen(mag_autocorr_filename,"w");
//	FILE * f_mag_temp = fopen("data/mag_temp.dat","w");

	FILE * f_en_bin = fopen(en_binning_filename,"w");
	FILE * f_en_autocorr = fopen(en_autocorr_filename,"w");
	FILE * f_en = fopen(en_filename,"a");
//	FILE * f_en_temp = fopen(en_temp_filename,"w");
	FILE * f_cv = fopen(cv_filename,"a");

	FILE * f_corr_row = fopen(corr_row_filename,"w");


	/*Start*/
	spin_init(matrix,nodes,N);

	/* TERMALIZZA*/
	evolve_therm(matrix,nodes,N,BETA);
	/***************/

	double en_mean =0;
	double en2_mean=0;
	double cv=0;
/* Vettori per il binning*/
	double * mag_vet_dati ;
	double * mag_vet_binnato;
	mag_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
	mag_vet_dati = malloc(sizeof(double)*ITERATION_MAX);

	double * en_vet_dati;
	double * en_vet_binnato;
	en_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
	en_vet_dati = malloc(sizeof(double)*ITERATION_MAX);

	double * chi_vet_binnato;
	chi_vet_binnato = malloc(sizeof(double)*(n_bin));

	double * cv_vet_binnato;
	cv_vet_binnato = malloc(sizeof(double)*(n_bin));

	double * en_autocorr = malloc(sizeof(double)*CORR_MAX);
	double * mag_autocorr = malloc(sizeof(double)*CORR_MAX);

	double * S_xt = malloc(sizeof(double)*N);
	double * S_yt = malloc(sizeof(double)*N);
	double * S_med_temp = malloc(sizeof(double)*N_CORR);
	double * S_dati = malloc(sizeof(double)*ITERATION_MAX*N_CORR);
	double * S_dati_binnati= malloc(sizeof(double)*(n_bin)*N_CORR);
	double * S_fin = malloc(sizeof(double)*N_CORR);
	double * S_var_fin=malloc(sizeof(double)*N_CORR);
	/*** Calcolo correlazione su righe e colonne*/
	vec_zeros(S_xt,N);
	vec_zeros(S_yt,N);
	vec_zeros(S_med_temp,N_CORR);
// Non serve, annullo già in binning:	vec_zeros(S_dati_binnati,(n_bin)*N_CORR);
	vec_zeros(S_fin,N_CORR);
	vec_zeros(S_var_fin,N_CORR);


	double mag_tmp;
	double en_tmp;

	/****** CICLO DI EVOLUZIONE: PRENDERE MISURE QUI */
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes,N,BETA);
		mag_tmp = magnetization(matrix,N);
		en_tmp = hamiltoniana(matrix,N);
		mag_vet_dati[iteration] = fabs(mag_tmp);
		en_vet_dati[iteration] = en_tmp;
		mag2_mean += mag_tmp*mag_tmp;
		mag_mean += fabs(mag_tmp);
		en_mean+=en_tmp;
		en2_mean += en_tmp*en_tmp;
		/* Correlazione righe e colonne*/
		for ( i = 0; i<N;i++){
			X_n[i] = sum_row(matrix,i,N);
			Y_n[i] = sum_col(matrix,i,N);
		}
		for (j = 0; j < N_CORR; ++j){
			for (i = 0; i < N;i++){
				S_xt[i]+= X_n[j]*X_n[(i+j)%N];
				S_yt[i]+= Y_n[j]*Y_n[(i+j)%N];
			}
		}
		for (i = 0; i < N;i++){
			S_xt[i]/=(double)N;
			S_yt[i]/=(double)N;
		}
		for (i=0;i<N_CORR;i++){
			S_med_temp[i] = S_xt[i]+S_yt[i]+S_xt[N-1-i]+ S_yt[N-1-i];
			S_med_temp[i] /=(4.0);
		}
		for (i = 0; i < N;i++){
			S_xt[i]=0;
			S_yt[i]=0;
		}
		for(j = 0;j<N_CORR;j++){
			S_dati[iteration*N_CORR+j] = S_med_temp[j];
		}


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
	fprintf(f_mag,"%.14e\t%.14e\n",BETA,mag_mean);
	fprintf(f_chi,"%.14e\t%.14e\n",BETA,chi);
	fprintf(f_en,"%.14e\t%.14e\n",BETA,en_mean);
	fprintf(f_cv,"%lf\t%lf\n",BETA,cv);

/**** BINNING E AUTOCORRELAZIONI */
	divideByScalar(mag_vet_dati,N*N,ITERATION_MAX);
	divideByScalar(en_vet_dati,N*N,ITERATION_MAX);
/* tutto il BINNING della Correlazione fra righe e colonne*/
	binning_mat(S_dati,S_dati_binnati,N_CORR,ITERATION_MAX,larghezza_bin);

	for(j=0;j<N_CORR;j++){
		for(i = 0; i< n_bin;i++){
			S_fin[j] += S_dati_binnati[i*N_CORR+j];
			S_var_fin[j] += S_dati_binnati[i*N_CORR+j]*S_dati_binnati[i*N_CORR+j];
		}
		S_var_fin[j] -= S_fin[j]*S_fin[j]/(double)(n_bin);
		S_fin[j]/= (double)(n_bin);
		S_var_fin[j] /=(double)(n_bin);
		S_var_fin[j] = sqrt((S_var_fin[j]/(n_bin)));
	}
	for ( i = 0; i<N_CORR;i++){
		fprintf(f_corr_row, "%d\t%.14e\t%.14e\n",i,S_fin[i],S_var_fin[i]);
	}



/* Calcolo necessario per stimare cosa scegliere come larghezza del bin!*/
	for ( larghezza_bin = 1; larghezza_bin < CORR_MAX ; larghezza_bin+=1){
		binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
		binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
		fprintf(f_mag_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(mag_vet_binnato,n_bin)/(double)(n_bin)));
		fprintf(f_en_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(en_vet_binnato,n_bin)/(double)(n_bin)));
	}

	/* Binning osservabili scalari*/
	binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
	binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
	fprintf(f_mag,"%.14e\t%.14e\t%.14e\n", BETA, meanOfDoubleArray(mag_vet_binnato,n_bin),
		sqrt(varianceOfDoubleArray(mag_vet_binnato,n_bin)));
	fprintf(f_en,"%.14e\t%.14e\t%.14e\n", BETA,meanOfDoubleArray(en_vet_binnato,n_bin),
		sqrt(varianceOfDoubleArray(en_vet_binnato,n_bin)));

	/* Calcolo autocorrelazione per le due grandezze */
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
/*
	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_en_temp,"%.14e\n",en_vet_dati[i]);
	}
	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_mag_temp,"%.14e\n",mag_vet_dati[i]);
	}
*/


 	/* Chiusura file */
	fclose(f_mag_bin);
	fclose(f_en_bin);
	fclose(f_en_autocorr);
	fclose(f_mag);
	fclose(f_chi);
	fclose(f_en);
	//fclose(f_en_temp);
	fclose(f_mag_autocorr);
	/* Free della memoria */
	free(chi_vet_binnato);
	free(cv_vet_binnato);
	free(mag_vet_dati) ;
	free(mag_vet_binnato);
	free(en_vet_dati);
	free(en_vet_binnato);
	free(en_autocorr);
	free(mag_autocorr);
	free(matrix);
	free(nodes);
	return EXIT_SUCCESS;
}
