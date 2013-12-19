#include <stdio.h>
#include "dynamicList.h"
#include <stdlib.h>


#define  N 256


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
	spinListInit(list,N);
	for ( i = 0 ; i< N ; i++){
		list[i].i = i;
		list[i].j = 0;
		list2[i].i = 0;
		list2[i].j = i;
		spinExtend( list+i, list2+i);
	}
	flowList( list , list);

	//init(list);
	//savePPM(list);
	exit (0);
}