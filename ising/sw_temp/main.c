#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 32
#define BETAJ 1


typedef struct Spin {
	int i,j;
	int spin;
	int cluster;
} Spin;

typedef struct Node {
	Spin * data;
	struct Node * next;
} Node;


/* Variabili Globali */
int cluster_max=-1;


Node * initCluster(Spin * s){
	Node * t;
	t = malloc(sizeof(Node));
	t->data = s;
	t->next = NULL;
	return t;
}
/* Aggiunge in testa alla lista --> FIFO
Ritorna la nuova testa */
Node * addToHead( Spin * s_new, Node * head){
	Node * t = malloc(sizeof(Node));
	if (t){	
		t->data = s_new;
		t->next = head;
		return t;
	}
	else {
		printf("Non c'è memoria -- add \n");
		exit(1);
	}
}

/* Ritorna la testa (nuova o vecchia che sia */
Node * removeElement (Spin * s , Node * head){
	Node * t = head;
	while (t->next){
		if (t->data == s)
			
		else
			t = t->next;
	}
}


void spin_init ( Spin * matrix){
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
		}
	}
}
void reset_cluster (Spin * matrix){
	int i,j;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			matrix[i*N+j].cluster = -1;
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

void fillCluster( Spin * matrix,Node * c){
	Spin * s_head;
	Spin * s_cycle;
	int x,y;
	int i,j;
	while (c){
		s_head = c->data;
		i = c->data->i;
		j = c->data->j;
		for ( x = -1; x< 1; x++){
			for ( y = -1; y<1; y++){
				s_cycle = matrix +  ((i + N + x)%N)*N + (j + N + y)%N;
				if ( s_cycle->cluster == -1){
					if (set_bond(s_cycle, s_head) == 1){
						//La nuova testa viene messa qui
							c = addToHead(s_cycle,c);
							s_cycle->cluster = cluster_max;
					}
				}
			}
		}

	}
}


void startClustering (Spin * matrix){
	int i,j;
	int x,y;
	int ii,jj;
	Spin * s_cycle;
	Spin * s_temp;
	Node * c;
	for ( i = 0; i< N; i++){
		for ( j = 0; j<N ; j++){
			s_cycle = matrix + i*N+j;
			if ( s_cycle->cluster == -1){
				cluster_max++;
				c = initCluster(s_cycle);
				c->data=NULL;
				//CHIAMA FUNZIONE DEL CLUSTER

			}
		}
	}
}



int main () {
	Spin * matrix = malloc(sizeof(Spin)*N*N);
	return 0;


}