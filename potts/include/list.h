#ifndef LIST_H
#define LIST_H

//#ifndef LIST_C
#include "cnum.h"
typedef struct Spin {
	int i,j;
	cNum spin;
	int cluster;
} Spin;

typedef struct Node {
	Spin * data;
	struct Node * next;
} Node;

typedef struct List{
	Node * head;
} List;

void  initCluster(Node * n, int n_cluster, List * l);
void addToHead( Node * n, List * l);
void removeElement (Node * del , List * list);

//#endif
#endif