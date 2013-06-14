#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 100, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(10 * n, m, 1.0, 0.001, 1000);
	double init[n][DIM], final[n][DIM];
	
	/*
	for(int k = 100; k < 200; k++){
	
		for(int i = 0; i < n; i++){
			
			for(int j = 0; j < DIM; j++){
				
				init[i][j] = net->nodes.at(i)->position[j];
			
			}
			
		}
		
		net->grow(1);
		
		for(int i = 0; i < n; i++){
			
			for(int j = 0; j < DIM; j++){
				
				final[i][j] = net->nodes.at(i)->position[j];
			
			}
			
		}
		
		for(int i = 0; i < n; i++){
			
			double dx = final[i][0] - init[i][0];
			double dy = final[i][1] - init[i][1];
			double dz = final[i][2] - init[i][2];
			double dr = sqrt(dx * dx + dy * dy + dz * dz);
			
			cout<<dr<<endl;
			
		}
		
	}
	*/
	
	for(int i = 0; i < 10 * n; i++){
		
		for(int j = 0; j < DIM; j++){
			
			cout<<net->nodes.at(i)->position[j]<<" ";
			
		}
		
		cout<<endl;
		
	}
	
}
