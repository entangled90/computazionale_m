#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "sw.h"
#include "cluster.h"

inline void savePPM(spin *  x)
{
    unsigned char white[3] = {255,255,255};
    unsigned char black[3] = {0,0,0};
   FILE *f = fopen("image.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", N, N);
    int i,j;
    for(i = N-1; i >= 0; i--)
        for(j = 0; j < N; j++){
            if(x[i*N+j].s == 1)
                fwrite(white, sizeof(unsigned char), 3, f);
            if(x[i*N+j].s == -1)
                fwrite(black, sizeof(unsigned char), 3, f);
        }
    fclose(f);
}


/* La head dei cluster va inizializzata con ID*/
int main (int argc, char *argv[]){	
	int i,j;
	spin * matrix;
	cluster *  c_head = malloc (sizeof(cluster));
	cluster * c_tmp;
	spin * s_tmp;
	if ( !( matrix = malloc (N*N*sizeof(spin))) )
		return 1;
	spin_init(matrix);
	c_head->id =-1;
	for (i=0;i<N;i++){
		new_cluster( c_head, &(matrix[i*N]));
//		for ( j = 0; j<N;j++){
		c_tmp = get_cluster(c_head,i);
		printf ("Cluster %d, pointer %p\n",c_tmp->id,c_tmp);
//		s_tmp = search_last_spin(c_tmp);

//			add_spin_to_cluster(get_cluster( c_head,i), matrix+i*N+j);
//		}
	}
	for (i=0;i<N;i++){
	// 	for ( j = 0; j<N;j++){	
		read_cluster(get_cluster(c_head,i));
			
	// 	}
	// }
	}
	//savePPM(matrix);
	exit(EXIT_SUCCESS);
}