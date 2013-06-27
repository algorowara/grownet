#include "../ngraph/nball.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <cstring>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int nmin = 10, nmax = 101, nstep = 1;
	long int m = 0, d = 3;
	long int samp = 100;
	double minPot[(nmax - nmin)/nstep + 1], maxPot[(nmax - nmin)/nstep + 1];
	memset(&minPot, 0, ((nmax - nmin)/nstep + 1) * sizeof(double));
	memset(&maxPot, 0, ((nmax - nmin)/nstep + 1) * sizeof(double));
	
	for(long int i = 0; i < samp; i++){
		
		NBall* b = new NBall(nmin, m, d);
		
		for(long int n = nmin; n < nmax; n += nstep){
			
			double pot = b->calculatePotential();
			
			if(i == 0){
				
				minPot[(n - nmin)/nstep] = pot;
				maxPot[(n - nmin)/nstep] = pot;
				
			}
			
			else if(pot < minPot[(n - nmin)/nstep]){
				
				minPot[(n - nmin)/nstep]  = pot;
				
			}
			
			else if(pot > maxPot[(n - nmin)/nstep]){
				
				maxPot[(n - nmin)/nstep] = pot;
				
			}
	
			b->grow(nstep);
			
		}
		
		delete b;
		
	}
	
	for(long int i = 0; i < (nmax - nmin)/nstep + 1; i++){
		
		if(minPot[i] != 0){
			
			cout<<(nmin + nstep * i)<<" "<<minPot[i]<<endl;
			
		}
		
	}
	
}
