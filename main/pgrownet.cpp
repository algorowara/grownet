#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

	long int n = 4, m = 2;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m);
	ofstream pgrowdata;
	pgrowdata.open("main/pgrowdata.txt", ios::out | ios::trunc);

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
