#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include "../ngraph/ngraph.h"
#include <omp.h>
#include <cfloat>
#include <iostream>

using namespace std;

int main(){
	
	long int n = 500, m = 3;
	long int s = 2;
	float gammin = 1.0, gamstep = 1.0, gammax = 40.0;
	
	for(long int d = 1; d < 10; d++){
		
		float dimmin = FLT_MAX;
		float bestgam;
		
		for(float g = gammin; g <= gammax; g += gamstep){
			
			float avgWeight = 0;
			
			for(long int i = 0; i < s; i++){
				
				NGraph* net = new NBall(m+1, m, d);
				net->baseGam = g;
				net->grow(n - (m+1));
				
				avgWeight += net->iterationWeights/s;
				
			}
			
			if(avgWeight < dimmin){
				
				dimmin = avgWeight;
				bestgam = g;
				
			}
			
		}
		
		cout<<d<<" "<<bestgam<<" "<<dimmin<<endl;
		
	}
	
}
