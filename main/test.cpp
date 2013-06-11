#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>

using namespace std;

int main(){
	
	long int n = 50, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = 0; j < DIM; j++){
			
			cout<<net->nodes.at(i)->position[j];
			
			if(j < DIM-1){
			
				cout<<" ";
				
			}
			
		}
		
		cout<<endl;
		
	}
	
}
