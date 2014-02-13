#ifndef CLUSTER_H
#define CLUSTER_H


#include <stdlib.h>

typedef struct spin_t {
	int i;
	int j;
	int s;
	int cluster;
	struct spin_t * next;
} spin ;

typedef struct cluster_t{
	int id;
	spin * head;
	struct cluster_t * next;
} cluster ;



cluster * get_last_cluster ( cluster * head){
	cluster * c;
	if (  (c = head) ){
		while ( c-> next ){
			c = c->next;
		}		
	return (c);
	}
	else{
		printf("Get_last_cluster null head pointer\n");
		exit(1);
	}
}

void add_cluster (cluster * head, cluster * c_new){
	cluster * c = head;
	c = get_last_cluster(head);
	c->next = c_new;
	c_new->next=NULL;
}


spin * get_last_spin (cluster * c){
	spin * tmp;
	if ( ! (tmp = c->head)){
		while (tmp->next ){
			tmp = tmp->next;
		}		
	}
	return tmp;
}


void add_spin_to_cluster ( cluster * c , spin * s_new){
	spin * s_old = get_last_spin(c);
	s_new->next =NULL;
	s_old->next = s_new;
	s_new->cluster = c->id;
}

/* Trova lo spin con indici i,j nel cluster */
spin * get_spin (cluster * c ,   int i ,   int j){
	spin * tmp = c->head;
	do {
		if ( (tmp->i == i) && (tmp->j == j)){
			return tmp;
		}
		else {
			tmp = tmp->next;
		}
	} while(tmp != NULL);
	return NULL;
}

/* Restituisce il cluster con id n*/
cluster * get_cluster (cluster * head,  int n ){
	cluster * tmp = head;
	int i = 0;
	while ( (i<n)){
		i++;
		tmp = tmp->next;
	}
	return tmp;
}

/*Crea un nuovo cluster che punta a head e lo aggiunge alla lista di cluster*/
cluster *  new_cluster (cluster * c_head , spin * s_head){
	cluster * c;
	if ( ( c = malloc(sizeof(cluster)) )) {
		cluster * tmp = get_last_cluster(c_head);
		c->id = tmp->id+1;
		c->head = s_head;
		s_head->next= NULL; //essendo da sola è anche l'ultima
		tmp->next = c;
		c->next = NULL;
		free(tmp);
	return c;
	}
	else{
		exit(1);
	}
}

void join_cluster (cluster * c1, cluster * c2){
	spin * tmp = get_last_spin(c1);
	tmp->next = c2->head;
	/* l'id del cluster è del primo*/
	c2->id = c1->id;
	}

void read_cluster ( cluster * c){
	spin * cursor;
	if ( (cursor = c->head)){
		while (cursor){
				//something useful
				printf ("Spin (%d,%d) = %d\n Next is %p \n",cursor->i,cursor->j,cursor->s, cursor->next);
				cursor = cursor->next; 	
		}		
	}
	else{
		printf("read_cluster wrong assignment\n");
		exit(1);
	}
}
#endif