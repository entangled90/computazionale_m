#ifndef SW_C
#define SW_C

#include "list.h"
#include "sw.h"
#include "constants-potts.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mtwist.h"
#include "cnum.h"

#define PI M_PI
#define BIN_WIDTH_DERIV 1000


int cluster_max = -1;
//Inizializza il reticolo con valori casuali dei 3 possibili stati
void spin_init ( Spin * matrix, Node * n, int N){
	int i,j;
	double tmp;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			tmp = mt_drand();
			if ( tmp < 1/3.0){
				matrix[i*N+j].spin = cNum_create(0);
			}
			else if (tmp < 2/3.0){
				matrix[i*N+j].spin = cNum_create(1);
			}
			else if (tmp <1.0){
				matrix[i*N+j].spin = cNum_create(2);
			}
			else{
				printf("Errore nei numeri casuali\n");
			}
			matrix[i*N+j].i = i;
			matrix[i*N+j].j = j;
			matrix[i*N+j].cluster = -1;
			n[i*N+j].data= matrix+i*N+j;
			n[i*N+j].next = NULL;
		}
	}
}
//Resetta i cluster
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
	if( s1->spin.index == s2->spin.index){
		if ( mt_drand() < (1- exp(-BETA)))
			return 1;
		else
			return 0;
	}
	else
		return 0;
}

//Riempie i cluster
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


void startClustering (Spin * matrix, Node * nodes, int N, float BETA){
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

//Ritorna la magnetizzazione per il numero di siti, senza il modulo. È  ancora un numero complesso.
inline cNum magnetization(Spin *x, int N){
	int i,j;
	cNum mag= { .r= 0.0 ,.i=0.0, .index=0};
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			mag = cNum_sum(mag,x[i*N+j].spin);
		}
	}
	return mag;
}


inline void savePPM(Spin * s, int N, char * filename)
{
    unsigned char r[3] = {255,0,0};
    unsigned char g[3] = {0,255,0};
    unsigned char b[3] = {0,0,255}
;   FILE *f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    int i,j;
    for(i = N-1; i >= 0; i--)
        for(j = 0; j < N; j++){
            if(s[i*N+j].spin.index == 0)
                fwrite(r, sizeof(unsigned char), 3, f);
            if(s[i*N+j].spin.index == 1)
                fwrite(g, sizeof(unsigned char), 3, f);
            if(s[i*N+j].spin.index == 2)
                fwrite(b, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}

//Sceglie in modo casuale la nuova configurazione in modo equiprobabile tra i tre stati per tutto il cluster
void flip_spin ( Spin * m, int N){
	int i,j;
	int * flipper;
	double tmp;
	flipper = (int *) malloc(sizeof(int)*(cluster_max+1));
	if (!flipper){
		printf("errore in flip spin: cluster_max = %d\n",cluster_max+1);
		exit(1);
	}
	for (i=0 ; i<cluster_max+1 ; i++){
		tmp = mt_drand();
		if (  tmp <= 1.0/(double)3.0)
			flipper[i] = 0;
		else if ( tmp <= 2.0/(double)3.0)
			flipper[i] = 1;
		else {
			flipper[i] = 2;
		}
	}
	for (i=0;i<N;i++){
		for(j = 0; j<N;j++){
//			m[i*N+j].spin =  cNum_create(flipper[ m[i*N+j].cluster ]);
		m[i*N+j].spin.r =cos((double)2*PI*flipper[m[i*N+j].cluster]/3.0);
		m[i*N+j].spin.i = sin((double)2*PI*flipper[m[i*N+j].cluster]/3.0);
		m[i*N+j].spin.index = flipper[m[i*N+j].cluster];
		}
	}
	free(flipper);
}
//Evoluzione di termalizzazione
void evolve_therm (Spin * matrix, Node * nodes, int N, float BETA){
	int iteration = 0;
	FILE * f_mag = fopen("data/mag_therm.dat","w");
	FILE * f_en = fopen("data/en_therm.dat","w");
	for ( iteration = 0 ; iteration < ITERATION_TEMP ; iteration++){
		startClustering(matrix,nodes,N,BETA);
		flip_spin(matrix,N);
		reset_cluster(matrix,nodes,N);
		fprintf(f_mag,"%d\t%e\n",iteration,mag_improved(matrix,N)/(double)(N*N));
		fprintf(f_en, "%d\t%e\n",iteration,hamiltoniana(matrix,N)/(double)(N*N));
	}
	fclose(f_mag);
	fclose(f_en);
}


//Evoluzione per la produzione, la raccolta dati avviene nel main.
void evolve( Spin * matrix, Node * nodes, int N, float BETA){
	reset_cluster(matrix,nodes,N);
	startClustering(matrix,nodes,N,BETA);
	flip_spin(matrix,N);
	//PRENDERE DATI QUI

}

//Calcola l'hamiltoniana, partendo da in alto a sinistra e usa i link a destra e in basso.
inline double hamiltoniana( Spin * s, int N){
	double ham=0;
	int a,b;
	for (a = 0; a<N ; a++){
		for ( b= 0; b<N; b++){
		if(s[((a+1+N)%N)*N +b].spin.index == s[a*N +b].spin.index){
				ham -= 1.0;
			}
	
			if(s[a*N+(b+1+N)%N].spin.index == s[a*N +b].spin.index){
				ham -= 1.0;
			}
		}
	}
	return (ham);
}
//Somma sulle righe per correlazione. ritorna un numero complesso
cNum sum_row(Spin * s, int row, int N){
	cNum sum = {0,0,0};
	int j = 0;
	for (j = 0; j<N ;j++){
		sum =  cNum_sum( sum, s[row*N+j].spin);
	}
	return cNum_mul(sum , cNum_create_real( 1.0/(double) N)) ;
}
//Somma sulle colonne per correlazione. ritorna un numero complesso

cNum sum_col(Spin * s, int col, int N){
	cNum sum = {0,0,0};
	int j = 0;
	for (j = 0; j<N ;j++){
		sum =  cNum_sum( sum, s[j*N+col].spin);
	}
	return cNum_mul(sum , cNum_create_real( 1.0/(double) N)) ;
}

//Restituisce la magnetizzazione nella definizione improved. Ossia come frazione di siti nel cluster + grande.
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
void print_state(Spin * s , int N){
	int i=0;
	FILE * f = fopen("data/positions.dat","w");
	for (i=0;i<N*N;i++){
		fprintf(f,"%.14e\t%.14e\n", s[i].spin.r,s[i].spin.i);
	}

	fclose(f);
}
#endif