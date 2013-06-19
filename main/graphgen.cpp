#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 2000, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m, 1.0, 0.2, 100);
	
	for(long int i = 0; i < n; i++){
		
		cout<<X(net->getNode(i))<<" "<<Y(net->getNode(i))<<" "<<Z(net->getNode(i))<<endl;
		
	}
	
}
