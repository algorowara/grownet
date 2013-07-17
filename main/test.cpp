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
	
	long int n = 100, m = 4;
	long int minDim = 1, maxDim = 10, stepDim = 1;
	long int minExp = 0, maxExp = 10, stepExp = 1;
	long int sample = 4;
	
	for(long int d = minDim; d <= maxDim; d += stepDim){
		
		for(long int exponent = minExp; exponent <= maxExp; exponent += stepExp){
			
			float averageCoef = 0;
			float averageWeight = 0;
				
			for(long int s = 0; s < sample; s++){
				
				NGraph* net = new NBall(m+1, m, d);
				net->forceExp = exponent;
				net->baseGam = 0.1;
				net->grow(n - (m+1));
				
				averageCoef += net->unweightedClusteringCoefficient()/sample;
				averageWeight += net->iterationWeights/sample;
				delete net;
				
			}		
			
			cout<<d<<" "<<exponent<<" "<<averageCoef<<" "<<averageWeight<<endl;
			
		}
		
	}
	
}
