#include "../ngraph/nball.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <cstring>
#include <sys/time.h>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	NBall* b;
	timeval start, end;
	long int n = 512;
	long int sample = 4;
	long int dmin = 2, dmax = 8;
	double gammin = 0.2, gammax = 5.0, gamstep = 0.2;
	
	for(long int d = 2; d <= dmax; d++){
		
		cout<<"Starting calculations for d = "<<d<<"..."<<endl;
		
		double best_gamma = 0;
		double dimensional_minimum = DBL_MAX;
		
		for(double gamma = gammin; gamma <= gammax; gamma += gamstep){
			
			double local_average = 0;
			
			for(long int i = 0; i < sample; i++){
				
				b = new NBall(n, d+1, d);
				
				local_average += ((double)b->iterationWeights)/sample;
				delete b;
				
			}
			
			if(local_average < dimensional_minimum){
				
				dimensional_minimum = local_average;
				best_gamma = gamma;
				
			}
		
		}
		
		cout<<"Best gamma for d = "<<d<<": "<<best_gamma<<", with an average of "<<dimensional_minimum<<" iteration-weights."<<endl;
				
	}
	
}
