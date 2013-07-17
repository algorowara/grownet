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
	
	long int nstart = 1000, nstep = 1000, nend = 10000;
	long int m = 6;
	long int sample = 5;
	
	ofstream file[6][3];
	
	for(long int i = 0; i < 6; i++){
		
		for(long int j = 0; j < 3; j++){
			
			stringstream ss;
			
			if(j == 0){
				
				ss<<"l";
				
			}
			
			else if(j == 1){
				
				ss<<"c_unweighted";
				
			}
			
			else{
				
				ss<<"c_weighted";
				
			}
			
			ss<<"_vs_n_";
			
			ss<<(((int)(i/2)) + 1);
			
			if(i%2 == 0){
				
				ss<<"sphere";
				
			}
			
			else{
				
				ss<<"ball";
				
			}
			
			ss<<"_m"<<m;
			ss<<".txt";
		
			file[i][j].open(ss.str().c_str(), ios::out | ios::trunc);
			
		}
		
		for(long int n = nstart; n <= nend; n += nstep){
			
			for(long int s = 0; s < sample; s++){
				
				NGraph* net;
				long int d = ((long int)(i/2)) + 1;
				
				if(i%2 == 0){
					
					net = new NSphere(n, m, d);
					
				}
				
				else{
					
					net = new NBall(n, m, d);
					
				}
				
				file[i][0]<<n<<" "<<net->averagePathLength()<<endl;
				file[i][1]<<n<<" "<<net->unweightedClusteringCoefficient()<<endl;
				file[i][2]<<n<<" "<<net->weightedClusteringCoefficient()<<endl;
				
			}
			
		}
		
		for(long int j = 0; j < 3; j++){
			
			file[i][j].close();
			
		}
		
	}
	
}
