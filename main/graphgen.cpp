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
	
	NBall* bestCurrent;
	NBall* current;
	NBall* base;
	long int nmin = 200, nmax = 2000, nstep = 100;
	long int sample = 16;
	long int d = 3, m = 4;
	double gammin = 0.05, gammax = 1.0, gamstep = 0.05;
	
	base = new NBall(nmin-nstep, m, d);
	
	for(long int n = nmin; n <= nmax; n += nstep){
		
		double best_gamma = 0;
		double worst_gamma = 0;
		double minimum = DBL_MAX;
		double maximum = DBL_MIN;
		
		for(double gamma = gammin; gamma <= gammax; gamma += gamstep){
			
			double local_average = 0;
			
			for(long int i = 0; i < sample; i++){
				
				current = new NBall(base);
				current->baseGamma = gamma;
				current->grow(nstep);
								
				local_average += ((double)(current->iterationWeights - base->iterationWeights))/sample;
				delete current;
				
			}
			
			if(local_average < minimum){
				
				minimum = local_average;
				best_gamma = gamma;
				
			}
			
			else if(local_average > maximum){
				
				maximum = local_average;
				worst_gamma = gamma;
				
			}
		
		}
		
		cout<<n<<" "<<best_gamma<<endl;
		base->baseGamma = best_gamma;
		base->grow(nstep);
		
	}
	
}
