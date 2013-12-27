#include <stdio.h>
#include "dynamicList.h"
#include <stdlib.h>


#define  N 32


void init(spin * s){
	int i;
	int j;
	float tmp;
	for (i = 0; i<N ; i++){
		for ( j = 0 ; j<N ; j++){
			s[i*N+j].i = i;
			s[i*N+j].j = j;
			tmp = rand()/(double)RAND_MAX;
			if ( tmp > 0.5){
				s[i*N+j].spin = +1;
			}
			else {
				s[i*N+j].spin = -1;
			}
		}
	}
}

inline void savePPM(spin * x)
{
    unsigned char white[3] = {255,255,255};
    unsigned char black[3] = {0,0,0};
 	FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    int i,j;
    for(i = N-1; i >= 0; i--)
        for(j = 0; j < N; j++){
            if(x[i*N+j].spin == 1)
                fwrite(white, sizeof(unsigned char), 3, f);
            if(x[i*N+j].spin == -1)
                fwrite(black, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}




int main (){
	spin * list;
	spin * list2;
	list = malloc(N*sizeof(spin));
	list2  = malloc (N*sizeof(spin));
	int i ;
	listInit(list,N);
	listInit(list2,N);
	for ( i = 0 ; i< N ; i++){
		list[i].index = i;
		list2[i].index = i+N;
	}
	spin * tmp = list->prev; // culo della prima lista
	spin * tmp2 = list2->prev; // culo della seconda

//	joinList(list, list2);
	if ( list->prev == tmp2 ){
		printf("funge\n");
	}
	for ( i = 0; i<N ; i++){
		extend ( list+N, list2+i);
	}

	spin * temp = list;
	while ( temp->next != list){
		printf("prev:%d\t indice: %d\t next:%d\n", temp->prev->index,temp->index,temp->next->index);
		temp = temp->next;
	}
	printf("indece del prossimo è: %d\n", temp->next->index);

	//init(list);
	//savePPM(list);
	exit (0);
}