#include <stdio.h>
#include <stdlib.h>


struct Spin {
	unsigned int i,j;
	unsigned int spin;
	int index;
};

struct SpinNode{
	struct Spin * node;
	struct Spin * next;
	struct Cluster * node_cluster;
};

struct Cluster{
	struct Node * head;-
};

struct ClusterNode{
	struct Cluster * node;
	struct Cluster * next;
};

struct ClusterList{
	struct ClusterNode * head;
};


typedef struct Spin Spin;
typedef struct Node Node;
typedef struct Cluster Cluster;
typedef struct ClusterNode ClusterNode;
typedef struct ClusterList ClusterList;

ClusterList * initClusterList (Cluster * c ){
	ClusterList * c_list = (ClusterList *) malloc(sizeof(ClusterList));
	ClusterNode * c_node = (ClusterNode * )malloc(sizeof(ClusterNode));
	c_list->head = c_node;
	c_node->node= &c;
	c_node->next = NULL;
	cluster_t->head = node_t;
	return cluster;
}

void deleteClusterList( ClusterList * head){
	ClusterNode ** temp = &(head->head);
	while (temp){
		free(*temp);
		temp = &((*temp)->next);
	}
	free(head);
}
void addClusterNode (ClusterNode * c_node, ClusterList * c_list){
	ClusterNode * temp = c_list->head;
	while(temp->next){
		temp = temp->next;
	}
	temp->next = c_node;
	c_node->next = NULL;
}

Cluster * initCluster (Node * s ){
	cluster_t = malloc(sizeof(Cluster));
	s->next = NULL;
	s->node_cluster = cluster_t;
	cluster_t->head = s;
	return cluster;
}

/* n_after è l'elemento dopo il quale si aggiunge n_new */
void addNodeAfterElement(Node * n_new, Node * n_after){
	n_new->next = n_after->next;
	n_after->next = n_new;
	n_new->node_t = n_after->node_cluster;
}

/* n_del va tolto dalla lista 
 0- elemento trovato
 1 - elemnto non trovat
 */
void deleteNode(Node * n_del){
	Node * temp = n_del->node_cluster;
	if ( temp == n_del){
		c->head = n_del->next;
		n_del->next = NULL;
	}
	else{
		while (temp->next != n_del !! temp->next == NULL){
			temp = temp->next;
		}
		if(temp->next != NULL){
			temp->next = n_del->next;
			n_del = NULL;
		}
	}
}
