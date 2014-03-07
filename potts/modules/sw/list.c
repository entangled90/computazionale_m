#ifndef LIST_C

#define LIST_C

#include <stdio.h>
#include <stdlib.h>
#include <list.h>


/*Il primo elemento dlela lista sarà l'ultimo ed ha NULL come puntatore a next*/
void  initCluster(Node * n, int n_cluster, List * l){
	n->next = NULL;
	l->head = n;
	n->data->cluster = n_cluster;
}
/* Aggiunge in testa alla lista --> FIFO
Ritorna la nuova testa */
void addToHead( Node * n, List * l){	
		n->next = l->head;
		l->head = n;
}

/* Ritorna la testa (nuova o vecchia che sia */
void removeElement (Node * del , List * list){
	Node * tmp = list->head;
	Node * tmp_prev = tmp;
	while (tmp){
		if (tmp == del){
			if ( del == list->head){
				list->head = list->head->next;
			}
			tmp_prev->next = tmp->next;
			return;
		}
		else{
			tmp_prev = tmp;
			tmp = tmp->next;
		}
	}
	printf ("Item not found\n");
	exit(1);
}


#endif