#define N 3

inline double r_squared_calc ( particle_s * list_0, particle_s * list_1){
	unsigned int i,k;
	double sum = 0;
	double rdiff[N];
	double distance, min;
	double rdiff2[N];
	int x,y;
	particle_s temp_part;
	for ( i = 0; i< NUMBER_OF_PARTICLES;i++){
		min = DBL_MAX;
		for ( x= -1; x < 2 ; x++){
			for ( y = -1; y<2 ; y++){
				temp_part = list_0[i];
				temp_part.position[0] += x;
				temp_part.position[1] += y;
				diff(list_1[i].position,temp_part.position,rdiff);
				distance = sqrt(scalar_prod(rdiff,rdiff));
				if( distance < min ){
					min = distance;
					for ( k = 0; k<N;k++){
						rdiff2[k] = rdiff[k];
					}
				}
			}
		}
		sum += scalar_prod(rdiff2,rdiff2);
	}
	return sum/NUMBER_OF_PARTICLES;
} 

void r_squared_save ( char * filename){
	FILE *f1 = fopen(filename, "w");
	fclose(f1);
	FILE *f = fopen(filename, "a");
	double sum=0;
	unsigned int delta,init;
	unsigned int count ;
	for ( delta = 1; delta  <  time_counted-1; delta++){
		sum = 0;
		count = 0;
		for ( init = 0; init +delta < time_counted; init++){
			sum += r_squared_calc( time_list+(init+delta)*NUMBER_OF_PARTICLES,time_list + init*NUMBER_OF_PARTICLES);
			count++;
		}
		sum /= (double) count;
//		sum=r_squared_calc(time_list,time_list+delta);
		fprintf(f,"%e\t%e\n",delta*DeltaT, sum);
	}
	fclose(f);
}
