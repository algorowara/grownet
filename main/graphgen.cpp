#include "../graph/vertex.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>

int main(){
	
	GrowingNetwork2D* net;
	long int n = 100000, m = 2;
	double average[n];
	int num_reps = 100;
	
	cout<<"Average Degree (k) as a function of Node Start Time"<<endl;
	
	for(int i = 0; i < n; i++){
		
		average[i] = 0;
		
	}
	
	for(int i = 0; i < num_reps; i++){
		
		net = new GrowingNetwork2D(n, m);
		
		for(int j = 0; j < n; j++){
		
			average[net->nodes.at(j)->getStartTime()] += net->K(j);
			
		}
		
		delete net;
		
	}
	
	for(int i = 0; i < n; i++){
		
		average[i] /= num_reps;
		cout<<i<<" "<<average[i]<<endl;
		
	}
	
}
