#ifndef SW_C
#define SW_C

#include "list.h"
#include "sw.h"
#include "constants-sw.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mtwist.h"

int cluster_max = -1;

void spin_init ( Spin * matrix, Node * n, int N){
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
void reset_cluster (Spin * matrix, Node * n, int N){
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
int set_bond (Spin * s1, Spin * s2, float BETA){
	if( s1->spin == s2->spin){
		if ( mt_drand() < (1- exp(-2*BETA)))
			return 1;
		else
			return 0;
	}
	else
		return 0;
}

void fillCluster( Spin * matrix, Node * nodes, List * l, int N, float BETA){
/* Matrice con le possibili direzioni: 4 coppie di vettori bidimensionali
	-1 	0	1 	0
	0 	-1 	0 	1 
*/
	int v[2][4] = { {-1,0,1,0}, {0,-1,0,1} };
	int i,j,k;
	int ii,jj;
	while (l->head){
		i = l->head->data->i;
		j = l->head->data->j;
		for ( k = 0 ; k<4;k++){
				ii = (i + N + v[0][k])%N;
				jj = (j + N + v[1][k])%N;
				if ( matrix[ii*N+jj].cluster < 0){
					if ( set_bond(matrix+i*N+j, matrix+ii*N+jj,BETA) ){
						//La nuova testa viene messa qui
							matrix[ii*N+jj].cluster = l->head->data->cluster;
							addToHead( nodes+ii*N+jj,l);
					}
				}
			}
		removeElement(nodes+i*N+j,l);
	}
}


void startClustering (Spin * matrix, Node * nodes, int N, double BETA){
	int i,j;
	List * cluster = malloc(sizeof(List));
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			if ( matrix[i*N+j].cluster == -1){
				cluster_max++;
//				printf("Cluster #%d\n",cluster_max);
				//Crea un nuovo cluster
				initCluster(nodes+i*N+j,cluster_max,cluster);
				//CHIAMA FUNZIONE DEL CLUSTER
				fillCluster(matrix,nodes,cluster,N,BETA);
			}
		}
	}
}

inline double magnetization(Spin *x, int N){
	int i,j;
	int mag=0;
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			mag += x[i*N+j].spin;
		}
	}
	return mag;
}


inline void savePPM(Spin * s, int N)
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

void drawCluster(Spin * s, int N){
	int i,j;
	for ( i = 0; i<N;i++){
		for ( j = 0; j<N;j++){
			printf("|\t%d\t", s[i*N+j].cluster);
		}
		printf("|\n");
	}
}

void drawSpin(Spin * s, int N){
	int i,j;
	for ( i = 0; i<N;i++){
		for ( j = 0; j<N;j++){
			if(s[i*N+j].spin == 1){
				printf("|\t+\t");				
			}
			else
				printf("|\t0\t");
		}
		printf("|\n");
	}
}

void flip_spin ( Spin * m, int N){
	int i,j;
	int * flipper;
	flipper = (int *) malloc(sizeof(int)*(cluster_max+1));
	if (!flipper){
		printf("errore in flip spin: cluster_max = %d\n",cluster_max+1);
		exit(1);
	}
	for (i=0 ; i<cluster_max+1 ; i++){
		if ( mt_drand() < 0.5)
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
void print_data (Spin * m,int N){
		drawCluster(m,N);
		printf("---------------------------------------------------------------------------\n");
		drawSpin(m,N);		
		printf("---------------------------------------------------------------------------\n\n");

}
void evolve_therm (Spin * matrix, Node * nodes, int N, float BETA){
	int iteration = 0;
	FILE * f_en_therm = fopen("data/en_therm.dat","w");
	FILE * f_mag_therm = fopen("data/mag_therm.dat","w");
	for ( iteration = 0 ; iteration < ITERATION_TEMP ; iteration++){
		startClustering(matrix,nodes,N,BETA);
		fprintf(f_en_therm, "%d\t%.10e\n",iteration, hamiltoniana(matrix,N)/(double)(N*N) );
		fprintf(f_mag_therm, "%d\t%.10e\n",iteration, fabs(magnetization(matrix,N)/(double)(N*N)) );
		flip_spin(matrix,N);
		reset_cluster(matrix,nodes,N);
	}
	fclose(f_en_therm);
	fclose(f_mag_therm);


}

void evolve( Spin * matrix, Node * nodes, int N, float BETA){
	reset_cluster(matrix,nodes,N);
	startClustering(matrix,nodes,N,BETA);
	flip_spin(matrix,N);
	//PRENDERE DATI QUI
}

inline double hamiltoniana( Spin * s, int N){
	double ham=0;
	int a,b;
	for (a = 0; a<N ; a++){
		for ( b= 0; b<N; b++){
			ham += -( s[ ((a+1+N)%N)*N + b].spin+ s[((a-1+N)%N)*N+b].spin
				+ s[a*N+(b+1+N)%N].spin+s[a*N + (b-1+N)%N].spin)*(s[a*N +b].spin); 
		}
	}
	return (ham/2.0);
}

double sum_row(Spin * s, int row, int N){
	double sum = 0;
	int j = 0;
	for (j = 0; j<N ;j++){
		sum += s[row*N+j].spin;
	}
	return sum /= ((double) N) ;
}

double sum_col(Spin * s, int col, int N){
	double sum = 0;
	int j = 0;
	for (j = 0; j<N ;j++){
		sum += s[j*N+col].spin;
	}
	return sum /= ((double) N) ;
}

double mag_improved(Spin *s , int N){
	int i;
	int max=0;
	int spin_in_cluster[cluster_max+1];
	for (i=0; i<cluster_max+1;i++){
		spin_in_cluster[i]=0;
	}
	for(i=0;i<N*N;i++){
		spin_in_cluster[s[i].cluster]++;
	}

	for (i=0; i<cluster_max+1;i++){
		if (spin_in_cluster[i]>max)
			max = spin_in_cluster[i];
	}

	return max;

}

#endif