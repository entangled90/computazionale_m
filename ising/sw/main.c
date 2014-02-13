#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// my includes
#include "sw.h"
#include "constants.h"
#include "mtwist.h"


int cluster_max = -1;
float BETA = 1;


int main ( int argc, char * argv[]) {
	double mag_mean=0;
	int iteration = 0;
	Spin * matrix = (Spin *) malloc(sizeof(Spin)*N*N);
	Node * nodes = (Node *) malloc(sizeof(Node)*N*N);
	/*Check for allocation*/
	if(!matrix || !nodes){
		printf("Cannot call malloc || MAIN ||\n");
		exit(1);
	}
	mt_seed();
	/*Check for command line arguments*/
	if (argc>1){
		BETA = atof(argv[1]);
	}
	else{
		printf("Inserire il valore di Beta\n");
		return 0;
	}

	/***** FILENAMES AND FILE OPENING ******/
	char mag_file[64] = "";
	snprintf(mag_file,64,"data/magnetization_mean_%d.dat",N);
	char m_file_long[64] = "";
	snprintf(m_file_long,64,"data/magnetization_long_%d.dat",N); 
	FILE * f_mag_mean = fopen(mag_file,"a");
	FILE * f_mag_long = fopen(m_file_long,"a");	

	/*Start*/
	srand(time(NULL));
	spin_init(matrix,nodes);
	evolve_therm(matrix,nodes);
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes);
		fprintf(f_mag_long, "%lf\t%lf\n", BETA, magnetization(matrix) );
		mag_mean += fabs(magnetization(matrix));
	}
	mag_mean /= (double)(ITERATION_MAX);
	fprintf(f_mag_mean,"%lf\t%lf\n",BETA,mag_mean);
	free(matrix);
	free(nodes);
	return 0;
}