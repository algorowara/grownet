#include "../ngraph/nsphere.h"
#include "../ngraph/nball.h"
#include "../growingnetwork3d/growingnetwork3d.h"
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
	
	long int nmin = 200, ninc = 50, nmax = 1000, m = 3, d = 2;
	NBall* net = new NBall(nmin, m, d);
	
	for(long int n = nmin; n <= nmax; n += ninc){
		
		stringstream ss;
		ofstream file;
		ss<<"../Analysis/Data/Points/Progression of Points on a 2-Ball/position_2ball_n"<<n<<"_p1.txt";
		file.open(ss.str().c_str(), ios::out | ios::trunc);
		
		for(long int i = 0; i < n; i++){
			
			for(long int j = 0; j < d; j++){
				
				file<<net->getNode(i)->position[j];
				
				if(j < d-1){
					
					file<<" ";
					
				}
				
			}
			
			file<<endl;
			
		}
		
		net->grow(ninc);
		
	}
}
