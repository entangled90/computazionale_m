#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "list.h"
#include <time.h>


#define N 128
#define ITERATION_MAX 50000
#define ITERATION_TEMP 10000

/* Variabili Globali */
int cluster_max=-1;
float BETAJ = 0;

void spin_init ( Spin * matrix, Node * n){
	int i,j;
	double tmp;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			tmp = rand()/(double)RAND_MAX;
			if ( tmp > 0.5){
				matrix[i*N+j].spin = +1;
			}
			else {
				matrix[i*N+j].spin = -1;
			}
			matrix[i*N+j].i = i;
			matrix[i*N+j].j = j;
			matrix[i*N+j].cluster = -1;
			n[i*N+j].data= matrix+i*N+j;
			n[i*N+j].next = NULL;
		}
	}
}
void reset_cluster (Spin * matrix, Node * n){
	int i,j;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			matrix[i*N+j].cluster = -1;
			n[i*N+j].next = NULL;
		}
	}
	cluster_max = -1;
}

/* Ritorna:
1- creare bond
2-non creare bond
*/
int set_bond (Spin * s1, Spin * s2){
	float tmp;
	if( s1->spin == s2->spin){
		tmp = rand()/(double)RAND_MAX;
		if ( tmp < (1- exp(-BETAJ)))
			return 1;
		else
			return 0;
	}
	else
		return 0;
}

void fillCluster( Spin * matrix, Node * nodes, List * l){
	int x,y;
	int i = 0;
	int j = 0;
	int ii,jj;
	while (l->head){
		i = l->head->data->i;
		j = l->head->data->j;
		for ( x = -1; x< 1; x++){
			for ( y = -1; y<1; y++){
				ii = (i + N + x)%N;
				jj = (j + N + y)%N;
				if ( matrix[ii*N+jj].cluster == -1){
					if ( set_bond(matrix+i*N+j, matrix+ii*N+jj) ){
						//La nuova testa viene messa qui
							addToHead( nodes+ii*N+jj,l);
							matrix[ii*N+jj].cluster = cluster_max;
					}
				}
			}
		}
		removeElement(nodes+i*N+j,l);
	}
}


void startClustering (Spin * matrix, Node * nodes){
	int i,j;
	List cluster;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			if ( matrix[i*N+j].cluster == -1){
				cluster_max++;
//				printf("Cluster #%d\n",cluster_max);
				//Crea un nuovo cluster
				cluster = initCluster(nodes+i*N+j,cluster_max);
				//CHIAMA FUNZIONE DEL CLUSTER
				fillCluster(matrix,nodes,&cluster);
			}
		}
	}
}

inline double magnetization(Spin *x){
	int i,j;
	int mag=0;
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			mag += x[i*N+j].spin;
		}
	}
	return mag/((double) (N*N));
}
inline void savePPM(Spin * s)
{
    unsigned char white[3] = {255,255,255};
    unsigned char black[3] = {0,0,0};
   FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    int i,j;
    for(i = N-1; i >= 0; i--)
        for(j = 0; j < N; j++){
            if(s[i*N+j].spin == 1)
                fwrite(white, sizeof(unsigned char), 3, f);
            if(s[i*N+j].spin == -1)
                fwrite(black, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}

void drawCluster(Spin * s){
	int i,j;
	for ( i = 0; i<N;i++){
		for ( j = 0; j<N;j++){
			printf("|\t%d\t", s[i*N+j].cluster);
		}
		printf("|\n");
	}
}

void drawSpin(Spin * s){
	int i,j;
	for ( i = 0; i<N;i++){
		for ( j = 0; j<N;j++){
			if(s[i*N+j].spin == 1){
				printf("|\t0\t");				
			}
			else
				printf("|\t+\t");
		}
		printf("|\n");
	}
}

void flip_spin ( Spin * m){
	int i,j;
	int * flipper;
	flipper = (int *) malloc(sizeof(int)*(cluster_max+1));
	if (!flipper){
		printf("errore in flip spin: cluster_max = %d\n",cluster_max+1);
		exit(1);
	}
	for (i=0 ; i<cluster_max+1 ; i++){
		if ( rand()/(float)RAND_MAX < 0.5)
			flipper[i] = 1;
		else
			flipper[i] = -1;
	}
	for (i=0;i<N;i++){
		for(j = 0; j<N;j++){
			m[i*N+j].spin *= flipper[m[i*N+j].cluster];
		}
	}
	free(flipper);
}

void evolve_therm (Spin * matrix, Node * nodes){
	int iteration = 0;
	for ( iteration = 0 ; iteration < ITERATION_TEMP ; iteration++){
		if (iteration %100 == 0){
			printf("Iteration #%d\n",iteration);
			printf("La magnetizzazione è %lf\n",magnetization(matrix));
			}
		startClustering(matrix,nodes);
		flip_spin(matrix);
	/*	drawCluster(matrix);
		printf("---------------------------------------------------------------------------\n");
		drawSpin(matrix);
	*/
		reset_cluster(matrix,nodes);
	}

}

void evolve( Spin * matrix, Node * nodes){
	startClustering(matrix,nodes);
	flip_spin(matrix);
	//PRENDERE DATI QUI
	reset_cluster(matrix,nodes);
}

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
		mag_mean += magnetization(matrix);
	}
	mag_mean /= (double)(ITERATION_MAX);
	FILE * f_mag_mean = fopen(mag_file,"a");
	fprintf(f_mag_mean,"%lf\t%lf\n",BETAJ,mag_mean);
	free(matrix);
	free(nodes);
	return 0;
}