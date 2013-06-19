#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 1000, n_max = 1050, m = 3, dim = 3;
	double gamma = 1.0, tolerance = 0.001;
	long int itr = 1000;
	long int sample_size = 20;
	double dr_min = 0.002;
	
	for(long int s = 0; s < sample_size; s++){
		
		GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tolerance, itr);
		
		for(int i = n; i < n_max; i++){
			
			double init[n][dim];
			
			for(long int j = 0; j < n; j++){
				
				for(long int k = 0; k < dim; k++){
					
					init[j][k] = net->getNode(j)->position[k];
					
				}
				
			}
			
			net->grow(1);
			
			for(long int j = 0; j < n; j++){
				
				double dx = net->getNode(j)->position[0] - init[j][0];
				double dy = net->getNode(j)->position[1] - init[j][1];
				double dz = net->getNode(j)->position[2] - init[j][2];
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				
				if(dr > dr_min){
				
					cout<<sqrt(dx * dx + dy * dy + dz * dz)<<endl;
							
				}
				
			}
			
		}
			
		delete net;
			
	}
	
}
