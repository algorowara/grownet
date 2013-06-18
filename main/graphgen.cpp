#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int m = 3;
	double gamma = 1.0, tolerance = 0.001;
	long int maxItr = 100;
	long int sample_size = 4;
	
	for(long int n = 50; n <= 1000; n += 50){
		
		for(long int i = 0; i < sample_size; i++){
			
			GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tolerance, maxItr);
			cout<<(net->N)<<" "<<(net->averagePathLength())<<endl;
			delete net;
			
		}
		
	}
	
}
