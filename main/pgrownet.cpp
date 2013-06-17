#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

<<<<<<< HEAD
	long int n = 1000, m = 2;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.05,.01,100);
=======
	long int n = 31, m = 2;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n, m, 0.01, 0.00001, 100000);
>>>>>>> 3bb2908de9606b41e5e598c20ad98d10fe2cc128
	ofstream pgrowdata;
	pgrowdata.open("pgrowdata.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		for(long int j = 0; j < DIM; j++){

			pgrowdata<<net->nodes.at(i)->position[j];

			if(j < DIM-1){

				pgrowdata<<" ";

			}

		}

		pgrowdata<<endl;

	}

	pgrowdata.close();
	cout << net->calculatePotential();
}
