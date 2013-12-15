#include<iostream>

using namespace std;

int main (){
	int a,b,j,i;
	int N = 1e3;
	a = 8;
	b =6;
	j=0;
	for( i=1; i<N; i=i+1 ){
		if (i%a==0 && i%b==0){
			j=j+1;
			cout<<"un multiplo di "<< a<< " e di "<< b<< " e': "<< i<<endl;
		}
	}
	cout<< "i multipli di " << a << " e "<< b << " sono: "<< j<<endl;
}
