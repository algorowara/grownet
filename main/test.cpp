#include "../ngraph/ngraph.h"
#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>
#include <cfloat>
#include <algorithm>
#include <iterator>

using namespace std;

int main(){
	
	for(long int f = 0; f <= 4; f++){
	
		long int n = 200, m = 4, d = 2;
		NSphere* net = new NSphere(m+1, m, d);
		net->baseGam = 0.01;
		net->forceExp = f;
		net->grow(n - (m+1));
		float* deg = net->degreeDistribution();
		
		for(long int i = 0; i < n; i++){
			
			if(deg[i] > 0){
				
				cout<<i<<" "<<deg[i]<<endl;
				
			}
			
		}
		
		cout<<endl;
		
		delete net;
		delete[] deg;
		
	}
	
}
