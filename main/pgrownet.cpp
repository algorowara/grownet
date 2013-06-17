#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

	long int n = 500, m = 3;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.05,.0001,10000);
	ofstream pgrowdata;
	pgrowdata.open("pgrowdata.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		for(long int j = 0; j < DIM; j++){

			pgrowdata<<net->getNode(i)->position[j];

			if(j < DIM-1){

				pgrowdata<<" ";

			}

		}

		pgrowdata<<endl;

	}

	pgrowdata.close();
	cout << net->calculatePotential() << endl;
	
	double* dist = net->degreeDistribution();

	for(long int i = 0; i < n; i++){
	
		cout<<i<<" "<<dist[i]<<endl;
	}

}
