#include <stdio.h>
#include <stdlib.h>


struct spin {
	unsigned int i,j;
	unsigned int spin;
	struct spin * next;
	struct spin * prev;
};


typedef struct spinList spinList;
typedef struct spin spin;


// potrebbe servire la lunghezza della lista

void spinListInit( spin * s , int n){
	int i;
	for ( i = 0 ; i< n; i++){
		s[i].next = s + (i+1)%n;
		s[i].prev = s + (i-1)%n;
	}
}


// s_new già allocato! s_new viene messo dopo s_old!
void spinExtend( spin * s_old , spin * s_new ){
	s_new->next = s_old->next;
	s_new->prev = s_old;
	s_old-> next = s_new;
}

void flowList ( spin * start, spin * stop ) {
	spin * temp;
	temp = start;
	while ( temp->next != stop){
		temp = temp->next;
		printf("Ho indici (%d,%d)\n", temp->i,temp->j);
	}
}