#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>
#include <ctime>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int n = 1000, m = 2;
	double degree[n];
	long int sample = 16;
	
	for(long int i = 0; i < n; i++){
		
		degree[i] = 0;
		
	}
	
	for(long int i = 0; i < sample; i++){
		
		GrowingNetwork2D* net = new GrowingNetwork2D(n, m);
		
		for(long int i = 0; i < n; i++){
			
			degree[net->getTime() - net->getNode(i)->getStartTime() - 1] += (net->K(i))/((double)sample);
			
		}
		
		delete net;
		
	}
	
	for(long int i = 0; i < n; i++){
		
		cout<<i<<" "<<degree[i]<<endl;
		
	}
	
}
