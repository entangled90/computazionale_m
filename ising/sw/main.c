#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// my includes
#include "sw.h"
#include "constants.h"


int cluster_max = -1;
float BETAJ = 0;


int main ( int argc, char * argv[]) {

	double mag_mean=0;
	int iteration = 0;
	char mag_file[64] = "";
	snprintf(mag_file,64,"data/magnetization_mean_%d.dat",N); 
	srand(time(NULL));
	Spin * matrix = (Spin *) malloc(sizeof(Spin)*N*N);
	Node * nodes = (Node *) malloc(sizeof(Node)*N*N);
	if(!matrix){
		printf("Cannot call malloc || MAIN ||\n");
		exit(1);
	}
	if (argc>1){
		BETAJ = atof(argv[1]);
	}
	spin_init(matrix,nodes);
	evolve_therm(matrix,nodes);
	for ( iteration=0;iteration<ITERATION_MAX; iteration++){
		evolve(matrix,nodes);
		mag_mean += fabs(magnetization(matrix));
	}
	mag_mean /= (double)(ITERATION_MAX);
	FILE * f_mag_mean = fopen(mag_file,"a");
	fprintf(f_mag_mean,"%lf\t%lf\n",BETAJ,mag_mean);
	free(matrix);
	free(nodes);
	return 0;
}