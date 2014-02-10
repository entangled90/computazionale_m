#include <stdio.h>
#include "list.h"
#include <math.h>

#define N 32
#define BETAJ 1


typedef struct Spin {
	int i,j;
	int spin;
	int checked;
	struct Spin * next;
}

typedef struct list {
	struct Spin * head;
	struct Spin * tail;
}



/* Variabili Globali */
int cluster_max=0;


/* Aggiunge in fondo alla lista*/
void add( Spin * s_new, list * l){
	(l->tail)->next = s_new;
	s_new->next = NULL;
	l->tail = s_new;
}

/*DEVI TENERE CONTO DI CLUSTER_MAX!*/
void join (list * l1, list * l2){
	
}

void spin init ( Spin * matrix){
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
			matrix[i*N+j].checked = 0;
			matrix[i*N+j].next = NULL;
		}
	}
}
void reset_cluster (Spin * matrix){
	int i,j;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			matrix[i*N+j].checked = 0;
			list[i*N+j]->head = NULL;
			list[i*N+j]->tail = NULL;
		}
	}
	cluster_max = 0;

}

/* Ritorna:
1- creare bond
2-non creare bond
*/
int set_bond (Spin * s1, Spin * s2){
	float tmp;
	if( s1->spin == s2->spin){
		tmp = rand()/(double)RAND_MAX;
		if ( tmp < (1- exp(BETAJ)))
			return 1;
		else
			return 0;
	}
	else
		return 0;
}

void clusterize(Spin * matrix, Spin ** cluster){
	int i,j;
	int x,y;
	int ii,jj;
	Spin * s_cycle;
	Spin * s_temp;
	matrix[0].cluster = 1;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			s_cycle = matrix + i*N+j;
			if ( s_cycle->cluster == 0){
				cluster[cluster_max] = s_cycle;
				cluster_max++;
			}
			for ( x = -1; x< 1; x++){
				for ( y = -1; y<1; y++){
					s_temp = matrix +  ((i + N + x)%N)*N + (j + N + y)%N;
					if ( s_temp->checked == 0){
						if (set_bond(s_cycle, s_temp) == 1){
							if(s_temp->cluster == 0){
								s_temp->cluster = s_cycle->cluster;
							}
							else{

							}
						}
					}
				}
			}
		}
	}
}



int main () {
	Spin * matrix = malloc(sizeof(Spin)*N*N);
	Spin * list = malloc(sizeof(list*)*N*N);


}