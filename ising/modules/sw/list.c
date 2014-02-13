#ifndef LIST_C

#define LIST_C

#include <stdio.h>
#include <stdlib.h>
#include <list.h>


/*Il primo elemento dlela lista sarà l'ultimo ed ha NULL come puntatore a next*/
List  initCluster(Node * n, int n_cluster){
	List t;
	n->next = NULL;
	t.head = n;
	n->data->cluster = n_cluster;
	return t;
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
			tmp_prev->next = tmp->next;
			if ( tmp == list->head){
				list->head = tmp->next;
			}
			return;
		}
		else{
			tmp = tmp->next;
		}
	}
	printf ("Item not found\n");
	exit(1);
}


#endif