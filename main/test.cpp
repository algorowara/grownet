#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 500, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
	ofstream points;
	points.open("main/points.txt", ios::out | ios::trunc);
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = 0; j < DIM; j++){
			
			points<<net->nodes.at(i)->position[j];
			
			if(j < DIM-1){
			
				points<<" ";
				
			}
			
		}
		
		points<<endl;
		
	}
	
	points.close();
	
}
