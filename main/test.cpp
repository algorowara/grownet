#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 4, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
	
	for(long int i = 0; i < n; i++){
		
		SpatialVertex* no = net->getNode(i);
		cout<<X(no)<<" "<<Y(no)<<" "<<Z(no)<<endl;
		
	}
			
}
