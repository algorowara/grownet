#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

	long int n = 31, m = 2;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.01,.00001,100000);
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
}
