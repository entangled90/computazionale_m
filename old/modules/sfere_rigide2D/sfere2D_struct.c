#ifndef sfere2D_struct_C
#define sfere2D_struct_C

typedef struct vector {
	double * v;
	int dim ;
	}

void vector_alloc ( vector * vec){
	vec.v = malloc (vec.dim * sizeof(double));
	}

vector * sum ( vector  * v1 , vector * v2){
	vector result;
	vector_alloc(result);
	int i = 0;
	for (i = 0; i < v1.dim ; i++){
		result.v[i] = v1.v[i] + v2.v[i];
	}
	return ( result);
	}












#endif
