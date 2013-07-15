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
	
	long int nstart = 50, nstep = 50, nend = 1000;
	long int m = 4;
	long int sample = 5;
	
	ofstream file[4][2];
	
	for(long int i = 0; i < 4; i++){
		
		for(long int j = 0; j < 2; j++){
			
			stringstream ss;
			
			if(j == 0){
				
				ss<<"l";
				
			}
			
			else{
				
				ss<<"c";
				
			}
			
			ss<<"_vs_n_";
			
			ss<<(((int)(i/2)) + 1);
			
			if(i%2 == 0){
				
				ss<<"sphere";
				
			}
			
			else{
				
				ss<<"ball";
				
			}
			
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
				
			}
			
		}
		
	}
	
}
