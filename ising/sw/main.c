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
	int i;
	/*Check for command line arguments*/
	if (argc>1){
		BETA = atof(argv[1]);
		N = atoi(argv[2]);
	}
	else{
		printf("Inserire il valore di Beta\n");
		exit(1);
	}
	Spin * matrix = (Spin *) malloc(sizeof(Spin)*N*N);
	Node * nodes= (Node *) malloc(sizeof(Node)*N*N);
	/*Check for allocation*/
	if(!matrix || !nodes){
		printf("Cannot call malloc || MAIN ||\n");
		exit(1);
	}
	float index_simulation = mt_drand();
	/***** FILENAMES AND FILE OPENING ******/
	char mag_file[64] = "";
	snprintf(mag_file,64,"data/mag_mean%d.dat",N);
	char chi_file[64] = "";
	snprintf(chi_file,64,"data/chi%d.dat",N); 
	char mag_file_binning[64] = "";
	snprintf(mag_file_binning,64,"data/binning/mag_N%d__B%.6lf__%.8lf.dat",N,BETA,index_simulation); 
	char en_file_binning[64] = "";
	snprintf(en_file_binning,64,"data/binning/en_N%d__B%.6lf__%.8lf.dat",N,BETA,index_simulation);
	char en_filename[64] = "";
	snprintf(en_filename,64,"data/en_N%d.dat",N);
	char en_corr_filename[64] = "";
	snprintf(en_corr_filename,64,"data/en_autocorrN%d_B%.8lf.dat",N,BETA);
	char en_temp_filename[64] = "data/en_temp.dat";
	char mag_autocorr_file[64] = "data/mag_corr.dat";
	FILE * f_mag_mean = fopen(mag_file,"a");
	FILE * f_chi = fopen(chi_file,"a");	
	FILE * mag_bin_file = fopen(mag_file_binning,"w");
	FILE * en_bin_file = fopen(en_file_binning,"w");
	FILE * en_corr_file = fopen(en_corr_filename,"w");
	FILE * en_file = fopen(en_filename,"a");
	FILE * en_temp_file = fopen(en_temp_filename,"w");
	FILE * f_mag_autocorr = fopen(mag_autocorr_file,"w");
	/*Start*/
	spin_init(matrix,nodes,N);
	evolve_therm(matrix,nodes,N,BETA);
	double mag_tmp;
	double en_tmp;
	double en_mean =0;
	double * mag_vet_dati ;
	double * mag_vet_binnato;
	double * en_vet_dati;
	double * en_vet_binnato;
	double * en_autocorr = malloc(sizeof(double)*CORR_MAX);
	double * mag_autocorr = malloc(sizeof(double)*CORR_MAX);
	en_vet_dati = malloc(sizeof(double)*ITERATION_MAX);
	mag_vet_dati = malloc(sizeof(double)*ITERATION_MAX);
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes,N,BETA);
		mag_tmp = magnetization(matrix,N);
		en_tmp = hamiltoniana(matrix,N)/(double)(N*N);
		mag_vet_dati[iteration] = fabs(mag_tmp);
		en_vet_dati[iteration] = en_tmp;
		mag2_mean += mag_tmp*mag_tmp;
		mag_mean += fabs(mag_tmp);
		en_mean+=en_tmp;
	}
	mag2_mean /= (double)(ITERATION_MAX);
	mag_mean /= (double)(ITERATION_MAX);
	chi = (mag2_mean - mag_mean*mag_mean)/((double)(N*N));
	mag_mean /=(double)(N*N);
	en_mean /= (double)(ITERATION_MAX);
	// BINNING E DATI VARI
	int larghezza_bin;
	mag_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX/2));
	en_vet_binnato = malloc(sizeof(double)*(ITERATION_MAX/2));
	for ( larghezza_bin = 2; larghezza_bin < 100 ; larghezza_bin+=2){
		binning(mag_vet_dati,mag_vet_binnato,ITERATION_MAX,larghezza_bin);
		binning(en_vet_dati,en_vet_binnato,ITERATION_MAX,larghezza_bin);
		fprintf(mag_bin_file,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(mag_vet_binnato,ITERATION_MAX/larghezza_bin)/(double)(ITERATION_MAX/larghezza_bin)));
		fprintf(en_bin_file,"%d\t%.14e\n", larghezza_bin,
			sqrt(varianceOfDoubleArray(en_vet_binnato,ITERATION_MAX/larghezza_bin)/(double)(ITERATION_MAX/larghezza_bin)));
	}
	autocorrelation(en_vet_dati,en_autocorr,ITERATION_MAX,CORR_MAX);
	if(en_corr_file){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(en_corr_file,"%d\t%.14e\n",i,en_autocorr[i]);
		}
	}
	autocorrelation(mag_vet_dati,mag_autocorr,ITERATION_MAX,CORR_MAX);
	if(f_mag_autocorr){
		for (i = 0; i<CORR_MAX;i++){
			fprintf(f_mag_autocorr,"%d\t%.14e\n",i,mag_autocorr[i]);
		}
	}



	for ( i=0;i<ITERATION_MAX;i++){
		fprintf(en_temp_file,"%.14e\n",en_vet_dati[i]);
	}
	printf("BETA: %lf , mag2: %.6e \t mag*mag: %.6e\t mag:%.6e\n",BETA,mag2_mean, mag_mean*mag_mean,mag_mean);
	fprintf(f_mag_mean,"%.14e\t%.14e\n",BETA,mag_mean);
	fprintf(f_chi,"%.14e\t%.14e\n",BETA,chi);
	fprintf(en_file,"%.14e\t%.14e\n",BETA,en_mean);
 	/* Chiusura file */
	fclose(mag_bin_file);
	fclose(en_bin_file);
	fclose(en_corr_file);
	fclose(f_mag_mean);
	fclose(f_chi);
	fclose(en_file);
	fclose(en_temp_file);
	/* Free della memoria */
	free(mag_vet_dati) ;
	free(mag_vet_binnato);
	free(en_vet_dati);
	free(en_vet_binnato);
	free(en_autocorr);
	free(matrix);
	free(nodes);
	return EXIT_SUCCESS;
}
