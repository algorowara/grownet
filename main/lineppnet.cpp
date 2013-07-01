#include "../pp1d/pp1d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

	long int n = 1000, m = 2;
	PPGrowingNetwork1D* net = new PPGrowingNetwork1D(n,m,1, .00001, 64);

	ofstream linedata;
	linedata.open("linedata.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		linedata<<net->getNode(i)->position[0]<<endl;

	}
	
	cout<<net->calculatePotential()<<endl;

}
