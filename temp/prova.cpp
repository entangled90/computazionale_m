#include<iostream>

using namespace std;

int main (){
	int N,i,somma;
	cout << "Inserisci un numero > 1\n";
	cin >> N;
	if ( N <=1){
		cout<<"Errore: N<=1\n";
	}
	somma = 0;
	for( i=1 ; i< N/2; i++){
		if ( N % i == 0){
			somma = somma + i;
		}
	}
	cout << "Il massimo numero perfetto è :" << somma << endl;
}
