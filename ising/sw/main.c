#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// my includes
#include "sw.h"
#include "list.h"
#include "constants-sw.h"
#include "mtwist.h"



int main ( int argc, char * argv[]) {
	float BETA = 1;
	int N = 40;
	double mag_abs_mean=0; // Valor medio del modulo della magnetizzazione
	double mag2_mean=0; // Valor medio di magnetizzazione al quadrato
	double mag_mean = 0; // valor medio di magnetizzazione
	double chi = 0;
	int iteration = 0;
	mt_seed();
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

	/***** FILENAMES AND FILE OPENING ******/
	char mag_file[64] = "";
	snprintf(mag_file,64,"data/magnetization_mean_%d.dat",N);
	char chi_file[64] = "";
	snprintf(chi_file,64,"data/chi_%d.dat",N); 
	FILE * f_mag_mean = fopen(mag_file,"a");
	FILE * f_chi = fopen(chi_file,"a");	

	/*Start*/
	spin_init(matrix,nodes,N);
	evolve_therm(matrix,nodes,N,BETA);
	double tmp ;
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes,N,BETA);
		tmp = magnetization(matrix,N);
		mag2_mean += tmp*tmp;
		mag_abs_mean += fabs(tmp);
		mag_mean  += tmp ;
	}
	mag_abs_mean /= (double)(ITERATION_MAX);
	mag2_mean /= (double)(ITERATION_MAX);
	mag_mean /= (double)(ITERATION_MAX);
	chi = (mag2_mean - mag_abs_mean*mag_abs_mean)/((double)(N*N));
	mag_mean /=(double)(N*N);
	mag_abs_mean /=(double)(N*N);
	printf("BETA: %lf , mag2: %.6e \t mag*mag: %.6e\t mag:%.6e\n",BETA,mag2_mean, mag_mean*mag_mean,mag_mean);
	fprintf(f_mag_mean,"%lf\t%lf\n",BETA,mag_abs_mean);
	fprintf(f_chi,"%lf\t%lf\n",BETA,chi);
	fclose(f_mag_mean);
	fclose(f_chi);
	free(matrix);
	free(nodes);
	return EXIT_SUCCESS;
}