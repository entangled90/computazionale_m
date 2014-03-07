#ifndef CNUM_H
#define CNUM_H
typedef struct cNum{
	double r;
	double i;
	short int index;
} cNum;
cNum cNum_create( int i);
cNum cNum_create_real( double a);
cNum cNum_sum (cNum a, cNum b);
cNum cNum_mul(cNum a , cNum b);
cNum cNum_conj(cNum a);
double cNum_Re(cNum a);
double cNum_Mod(cNum a);
#endif

