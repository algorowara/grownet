#include "../ngraph/nball.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>

using namespace std;

int main(){
	
	long int nmin = 1000, nmax = 10000, nstep = 500;
	ofstream file;
	
	for(long int d = 2; d < 12; d += 2){
		
		long int m = d+1;
		NBall* b = new NBall(nmin, m, d);
		
		stringstream sstream;
		sstream<<"c_vs_n_with_error_"<<d<<"ball_m"<<m<<".txt";
			
		file.open(sstream.str().c_str(), ios::out | ios::trunc);
		
		for(long int n = nmin; n < nmax; n += nstep){
			
			double* data = b->averageClusteringCoefficientWithError();
			
			file<<n<<" "<<data[0]<<" "<<data[1]<<endl;
			b->grow(nstep);
			
			if(n + nstep >= nmax){
				
				data = b->averageClusteringCoefficientWithError();
				
				file<<(n + nstep)<<" "<<data[0]<<" "<<data[1]<<endl;
				
			}
			
		}
		
		file.close();
		cout<<"Completed "<<d<<"-Ball data collection."<<endl;
		
	}
	
}
