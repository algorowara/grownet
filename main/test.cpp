#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 1000, m = 3;
	double gamma = 1.0, tolerance = 0.001;
	long int itr = 1000;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tolerance, itr);
	
	double* dist = net->degreeDistribution();
	
	for(long int i = 0; i < n; i++){
	
		cout<<i<<" "<<dist[i]<<endl;
		
	}
	
}
