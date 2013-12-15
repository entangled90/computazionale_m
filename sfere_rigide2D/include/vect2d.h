#ifndef VECT2D_H
#define VECT2D_H

#ifndef VECT2D_C
#define D 2


inline void  sum ( double * v1 , double * v2, double * v_sum){
	int i = 0;
	for (i = 0; i < D ; i++){
		v_sum[i] = v1[i] + v2[i];
	}
}


inline void diff (double * v1 , double *v2, double * v_diff){
	int i = 0;
	for (i = 0; i < D ; i++){
		v_diff[i] = v1[i] - v2[i];
	}
}

	
inline double scalar_prod ( double * v1 , double * v2){
	int i = 0;
	double sum  = 0;
	for ( i = 0; i< D ; i++){
		sum += v1[i]*v2[i];
	}
	return (sum);
}
inline void scalar_mult  ( double scalar , double* vec){
	int i = 0;
	for( i = 0 ; i< D ; i++){
		vec[i] *= scalar;
	}
}

#endif
#endif

