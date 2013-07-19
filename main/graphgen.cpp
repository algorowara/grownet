#include "../ngraph/nsphere.h"
#include "../ngraph/nball.h"
#include "../growingnetwork3d/growingnetwork3d.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <cstring>
#include <sstream>
#include <sys/time.h>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int n = 200, m = 4, d = 3;
	
	NBall* net = new NBall(m+1, m, d);
	net->baseGam = 4.0;
	net->attractiveForceExp = 2;
	net->forceExp = 3;
	net->grow(n - (m+1));
	float* deg = net->degreeDistribution();
	
	
	for(long int i = 0; i < n; i++){
		
		if(deg[i] > 0){
			
			cout<<i<<" "<<deg[i]<<endl;
			
		}
		
	}
	
}
