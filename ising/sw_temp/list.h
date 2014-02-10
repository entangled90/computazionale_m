
 struct point {
	int i;
	int j;
	int spin;
	int cluster;
};
typedef struct point point;

 struct Node {
	point p;
	point * next;
};


typedef struct Node Node;
struct List {
	Node firstNode;
};

/*
Inserisce n_new dopo n_old
*/
void insertAfter (Node n_old , Node n_new ){
	if(!n_old){
		n_new.next = n_old.next;
		n_old.next = n_new;
	}
}

/*
Elimina il nodo
*/
void deleteNode (Node n){
	if (!n){
		
	}
}
