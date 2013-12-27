#include <stdio.h>
#include <stdlib.h>


struct spin {
	unsigned int i,j;
	unsigned int spin;
	int index;
	struct spin * next;
	struct spin * prev;
};



typedef struct spin spin;


// potrebbe servire la lunghezza della lista

void listInit( spin * s , int n){
	int i;
	for ( i = 0 ; i< n; i++){
		s[i].next = s + (i+1)%n;
		s[i].prev = s + (i-1)%n;
	}
}


/* s_new già allocato! s_new viene messo dopo s_old! */
void extend( spin * s_old , spin * s_new ){
	s_new->next = s_old->next;
	s_new->prev = s_old;
	s_old-> next = s_new;
}

/* head sono i puntatori all'inizio delle liste (che concettualmente non ci sono visto che le liste son periodiche)
	ma facendo head.prev ottieni la coda automaticamente.
	Nella memoria le liste son periodiche, nell'utilizzo possono essere utilizzate come orientate!
  attacco la testa della seconda lista alla coda della prima.
*/
void joinList ( spin * head1, spin * head2){
	spin * temp;
	spin * temp2;
	temp = head1->prev; // tail1
	temp2 = head2->prev;
	temp2->next = head1;
	head1->prev = temp2; // nuovo tail1 è tail2
 	head2->prev = temp; //nuovo tail2 è tail1
	temp->next = head2;
}

void flowList ( spin * start, spin * stop ) {
	spin * temp;
	temp = start;
	while ( temp->next != stop){
		printf("prev:%d\t indice: %d\t next:%d\n", temp->prev->index,temp->index,temp->next->index);
		temp = temp->next;
	}
	printf("indece del prossimo è: %d\n", temp->next->index);
}