#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>

using namespace std;

int main(){
	
	GrowingNetwork2D* net = new GrowingNetwork2D(100, 3);
	
	for(long int i = 0; i < 100; i++){
		
		cout<<(net->getNode(i) != NULL);
		
	}
	
	cout<<endl;
	
}
