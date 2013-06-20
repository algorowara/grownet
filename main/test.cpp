#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 1000, m = 3;
	GrowingNetwork3D* net;
	double* data;
	long int sample_size = 16;
	double avg[n];
	
	for(long int i = 0; i < n; i++){
		
			avg[n] = 0;
		
	}
	
	for(long int i = 0; i < sample_size; i++){
		
		net = new GrowingNetwork3D(n, m);
		data = net->edgeAgeVsBetweenness();
		
		for(long int j = 0; j < n; j++){
			
			avg[i] += data[i]/sample_size;
			
		}
		
		delete[] data;
		
	}
	
	for(long int i = 0; i < n; i++){
		
		cout<<i<<" "<<avg[i]<<endl;
		
	}
			
}
