#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 2000, m = 2, dim = 3;
	double gamma = 1.0, tolerance = 0.001;
	long int itr = 1000;
	
	for(long int i = 0; i < 5; i++){
		
		net = new GrowingNetwork3D(n, m, gamma, tolerance, itr);
		double init[n][dim];
		
		for(long int j = 0; j < n; j++){
			
			init[j] = net->getNode(j)->position;
			
		}
		
		net->grow(1);
		
		for(long int j = 0; j < n; j++){
			
			dx = net->getNode(j)->position[0] - init[j][0];
			dy = net->getNode(j)->position[1] - init[j][1];
			dz = net->getNode(j)->position[2] - init[j][2];
			
			cout<<sqrt(dx * dx + dy * dy + dz * dz)<<endl;
			
		}
		
		delete net;
		
	}
	
}
