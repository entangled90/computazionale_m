#include<iostream>
/*
 * 
 * 
 * mail
 */
using namespace std;

int main (){
	int i;
	double somma;
	double media;
	double v[20] = {12,23,34,1,12,93,8,7,1,2,9,8,3,7,9,1,8,2,3,32};
	somma = 0;
	for ( i = 0; i<20 ; i++){
		somma = somma + v[i];
	}
	media = somma / 20;
	cout << "la media è : "<< media<<endl;
}
