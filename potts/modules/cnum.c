#ifndef CNUM_C 
#define CNUM_C

#include "cnum.h"
#include <math.h>
#include <stdio.h>

#define PI atan(1)

#define cNUM_CREATE(x) 
cNum cNum_create( int i){
	cNum result;
	result.r = cos(2*PI*i/3.0);
	result.i = sin(2*PI*i/3.0);
	result.index = i;
	return (result);
}

cNum cNum_create_real( double a){
	cNum result = {a,0,0};
	return result;
}

cNum cNum_sum (cNum a, cNum b){
	cNum result;
	result.r = a.r + b.r;
	result.i = a.i + b.i;
	return result;
}

cNum cNum_mul(cNum a , cNum b){
	cNum result;
	result.r = a.r*b.r - a.i*b.i;
	result.i = a.i*b.r + a.r*b.i;
	result.index = (a.index+b.index)%3;
	return result;
}

cNum cNum_conj(cNum a){
	cNum result;
	result.r= a.r;
	result.i = -a.i;
	result.index = (3 - a.index)%3;
	return result;
}

double cNum_Re(cNum a){
	return a.r;
}

double cNum_Mod(cNum a){
	return (sqrt(a.r*a.r+ a.i*a.i));
}



#endif