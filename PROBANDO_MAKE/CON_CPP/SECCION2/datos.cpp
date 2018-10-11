#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream> 
#include <ctime>
#include <string>
#include <sstream>
using namespace std;
#define PI 3.141516
 

int main(){
	ofstream File; 
	File.open("output.txt");
	int N,i,j;
	
	cout<< "De un valor de Numero de simulaciones N"<< endl;
	cin >> N;
	double suma=0;
	
	for(i=0; i<N+1;i++){
	File << suma << endl;	
	suma+=2*PI/N;
	}
	File.close();
	
	


return 0;
}
