/*** Il programma è sostanzialmente identico a Metropolis eccetto per l'algoritmo che genera nuove configurazioni 
(che è appunto Swendsen-Wang)
**/
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

#define CORR_MAX 50
#define CORR_ESTREMO 15
#define BIN_WIDTH_DERIV 3000


int main ( int argc, char * argv[]) {
	double BETA = 1;
	int N = 40;
	int iteration = 0;
	mt_seed();
	int i,j;
	int tau_corr=5;
	int larghezza_bin = 30;
	int n_bin = ITERATION_MAX/larghezza_bin;
	printf("%d\n",n_bin);
	if (n_bin ==0){
		printf("Numero BIN = 0\n");
		exit(1);
	}
	int n_bin_deriv = ITERATION_MAX/BIN_WIDTH_DERIV;

	/*Check for command line arguments*/
	if (argc>1){
		BETA = atof(argv[1]);
		N = atoi(argv[2]);
	}
	else{
		printf("Inserire il valore di Beta e N\n");
		exit(1);
	}
	int N_CORR = N/5;
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
	int rand_sig = (int) (mt_drand()*100000000)%9999;
/* ------------------MAGN*/
	char mag_filename[64] = "";
	snprintf(mag_filename,64,"data/mag_mean%d.dat",N);
	char chi_filename[64] = "";
	snprintf(chi_filename,64,"data/chi%d.dat",N); 
	char mag_binning_filename[64] = "";
	snprintf(mag_binning_filename,64,"data/binning/mag_N%d__B%.8lf.dat",N,BETA); 
	char mag_autocorr_filename[64] = "";
	snprintf(mag_autocorr_filename,64,"data/mag_corr/mag_autocorr%.4dN%d_B%.8lf.dat",rand_sig,N,BETA);
	char mag_tau_filename[64] = "";
	snprintf(mag_tau_filename,64,"data/tau_magN%d.dat",N);

/*----------------------_ENergia*/
	char en_binning_filename[64] = "";
	snprintf(en_binning_filename,64,"data/binning/en_N%d__B%.8lf.dat",N,BETA);
	char en_filename[64] = "";
	snprintf(en_filename,64,"data/en_N%d.dat",N);
	char en_autocorr_filename[64] = "";
	snprintf(en_autocorr_filename,64,"data/en_corr/en_autocorr%.4dN%d_B%.8lf.dat",rand_sig,N,BETA);
//	char en_temp_filename[64] = "data/en_temp.dat";
	char cv_filename[64]="";
	snprintf(cv_filename,64,"data/cv%d.dat",N);
	char en_tau_filename[64] = "";
	snprintf(en_tau_filename,64,"data/tau_enN%d.dat",N);

	char corr_row_filename[64]="";
	snprintf(corr_row_filename,64,"data/corr_row/corr_row_N%dB%.8lf.dat",N,BETA);


/*** I File sono chiamati: f_$(Nomestringa) ****/
	FILE * f_mag = fopen(mag_filename,"a");
	FILE * f_chi = fopen(chi_filename,"a");	
	FILE * f_mag_bin = fopen(mag_binning_filename,"w");

	FILE * f_mag_autocorr = fopen(mag_autocorr_filename,"w");
	FILE * f_mag_tau = fopen(mag_tau_filename,"a");
	FILE * f_mag_temp = fopen("data/mag_temp.dat","w");	FILE * f_en_bin = fopen(en_binning_filename,"w");
	FILE * f_en_autocorr = fopen(en_autocorr_filename,"w");
	FILE * f_en = fopen(en_filename,"a");
//	FILE * f_en_temp = fopen(en_temp_filename,"w");
	FILE * f_cv = fopen(cv_filename,"a");
	FILE * f_en_tau = fopen(en_tau_filename,"a");
	FILE * f_corr_row = fopen(corr_row_filename,"w");


	/*Start*/
	spin_init(matrix,nodes,N);

	/* TERMALIZZA*/
	evolve_therm(matrix,nodes,N,BETA);
	/***************/

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
	chi_vet_binnato = malloc(sizeof(double)*(n_bin_deriv));

	double * cv_vet_binnato;
	cv_vet_binnato = malloc(sizeof(double)*(n_bin_deriv));

	double * en_autocorr = malloc(sizeof(double)*CORR_MAX);
	double * mag_autocorr = malloc(sizeof(double)*CORR_MAX);

	double * S_xt = malloc(sizeof(double)*N);
	double * S_yt = malloc(sizeof(double)*N);
	double * S_med_temp = malloc(sizeof(double)*N_CORR);
	double * S_dati = malloc(sizeof(double)*ITERATION_MAX*N_CORR);
	double * S_dati_binnati= malloc(sizeof(double)*(n_bin)*N_CORR);
	double * S_fin = malloc(sizeof(double)*N_CORR);
	double * S_var_fin=malloc(sizeof(double)*N_CORR);
	double S_test[N_CORR];

	/*** Calcolo correlazione su righe e colonne*/
	vec_zeros(S_xt,N);
	vec_zeros(S_yt,N);
	vec_zeros(S_med_temp,N_CORR);
// Non serve, annullo già in binning:	vec_zeros(S_dati_binnati,(n_bin)*N_CORR);
	vec_zeros(S_fin,N_CORR);
	vec_zeros(S_var_fin,N_CORR);
	vec_zeros(S_test,N_CORR);
	/****** CICLO DI EVOLUZIONE: PRENDERE MISURE QUI */
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes,N,BETA);
		mag_vet_dati[iteration] = fabs(mag_improved(matrix,N));
		en_vet_dati[iteration] = hamiltoniana(matrix,N);
		/* Correlazione righe e colonne*/
		for ( i = 0; i<N;i++){
			X_n[i] = sum_row(matrix,i,N);
			Y_n[i] = sum_col(matrix,i,N);
		}
		//Invarianza per traslazioni
		for (j = 0; j < N; ++j){
			for (i = 0; i < N;i++){
				S_xt[i]+= X_n[j]*X_n[(i+j)%N];
				S_yt[i]+= Y_n[j]*Y_n[(i+j)%N]; 
			}
		}
		for (i = 0; i < N;i++){
		S_xt[i]/=(double)N;
		S_yt[i]/=(double)N;
		}	
		// Invarianza per rotazioni e ciclicità+inversione
		for (i=0;i<N_CORR;i++){
			S_med_temp[i] = S_xt[i]+S_yt[i]+S_xt[N-1-i]+ S_yt[N-1-i];
			S_med_temp[i] /=(4.0);
	//		S_test[i] += S_med_temp[i];
		}
		for(j = 0;j<N_CORR;j++){
			S_dati[iteration*N_CORR+j] = S_med_temp[j];
		}
		for (i = 0; i < N;i++){
			S_xt[i]=0;
			S_yt[i]=0;
		}

	}

	divideByScalar(S_test,ITERATION_MAX,N_CORR);
/**** BINNING E AUTOCORRELAZIONI */
/* tutto il BINNING della Correlazione fra righe e colonne*/
	binning_mat(S_dati,S_dati_binnati,N_CORR,ITERATION_MAX,larghezza_bin);

	for(j=0;j<N_CORR;j++){
		for(i = 0; i< n_bin;i++){
			S_fin[j] += S_dati_binnati[i*N_CORR+j];
			S_var_fin[j] += S_dati_binnati[i*N_CORR+j]*S_dati_binnati[i*N_CORR+j];
		}
		S_fin[j]/= (double)(n_bin);
		S_var_fin[j] /=(double)(n_bin);
		S_var_fin[j] -= S_fin[j]*S_fin[j];
		S_var_fin[j] = sqrt((S_var_fin[j]/(double)(n_bin)));
	}
	for ( i = 0; i<N_CORR;i++){
		fprintf(f_corr_row, "%d\t%.14e\t%.14e\n",i,S_fin[i],S_var_fin[i]);
	}


	/* Binning osservabili scalari*/
	binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
	binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
	binning_deriv(mag_vet_dati,chi_vet_binnato,ITERATION_MAX,BIN_WIDTH_DERIV);
	binning_deriv(en_vet_dati,cv_vet_binnato,ITERATION_MAX,BIN_WIDTH_DERIV);

	divideByScalar(mag_vet_binnato,N*N,n_bin);
	divideByScalar(en_vet_binnato,N*N,n_bin);
	divideByScalar(chi_vet_binnato,N*N,n_bin_deriv);
	divideByScalar(cv_vet_binnato,N*N,n_bin_deriv);

	fprintf(f_mag,"%.8lf\t%.14e\t%.14e\n", BETA, meanOfDoubleArray(mag_vet_binnato,n_bin),
		sqrt(varianceOfDoubleArray(mag_vet_binnato,n_bin)/n_bin));
	fprintf(f_en,"%.8lf\t%.14e\t%.14e\n", BETA,meanOfDoubleArray(en_vet_binnato,n_bin),
		sqrt(varianceOfDoubleArray(en_vet_binnato,n_bin)/n_bin));
	fprintf(f_cv,"%.8lf\t%.14e\t%.14e\n", BETA,meanOfDoubleArray(cv_vet_binnato,n_bin_deriv),
		sqrt(varianceOfDoubleArray(cv_vet_binnato,n_bin_deriv)/n_bin_deriv));
	fprintf(f_chi,"%.8lf\t%.14e\t%.14e\n", BETA,meanOfDoubleArray(chi_vet_binnato,n_bin_deriv),
		sqrt(varianceOfDoubleArray(chi_vet_binnato,n_bin_deriv)/n_bin_deriv));

	divideByScalar(mag_vet_dati,N*N,ITERATION_MAX);
	divideByScalar(en_vet_dati,N*N,ITERATION_MAX);
//	mag_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
//	en_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX));
	for ( larghezza_bin = 1; larghezza_bin < 100 ; larghezza_bin+=1){
		binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
		binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
		fprintf(f_mag_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(mag_vet_binnato,ITERATION_MAX/larghezza_bin)/(double)(ITERATION_MAX/larghezza_bin)));
		fprintf(f_en_bin,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(en_vet_binnato,ITERATION_MAX/larghezza_bin)/(double)(ITERATION_MAX/larghezza_bin)));
	}


	/* Calcolo autocorrelazione per le due grandezze */
	autocorrelation(en_vet_dati,en_autocorr,ITERATION_MAX,CORR_MAX);
	if(f_en_autocorr){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(f_en_autocorr,"%d\t%.14e\n",i,en_autocorr[i]);
		}
	}
	double tau_en = 0.5;
	for (i=0;i<CORR_ESTREMO;i++){
		tau_en+=en_autocorr[i];
	}
	fprintf(f_en_tau, "%.14e\t%.14e\n",BETA,tau_en);

	autocorrelation(mag_vet_dati,mag_autocorr,ITERATION_MAX,CORR_MAX);
	if(f_mag_autocorr){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(f_mag_autocorr,"%d\t%.14e\n",i,mag_autocorr[i]);
		}
	}

	double tau_mag = 0.5;
	for (i=0;i<CORR_ESTREMO;i++){
		tau_mag+=mag_autocorr[i];
	}
	fprintf(f_mag_tau, "%.14e\t%.14e\n",BETA,tau_mag);
/*
	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_en_temp,"%.14e\n",en_vet_dati[i]);
	}
*/
	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(f_mag_temp,"%.14e\n",mag_vet_dati[i]);
	}
 	/* Chiusura file */
	fclose(f_mag_bin);
	fclose(f_en_bin);
	fclose(f_en_autocorr);
	fclose(f_mag);
	fclose(f_chi);
	fclose(f_en);
	fclose(f_mag_temp);
	fclose(f_mag_tau);
	fclose(f_en_tau);
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
