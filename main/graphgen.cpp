#include "../ngraph/nsphere.h"
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
	
	long int nstart = 500, nend = 1000, m = 3, d = 2;
	NBall* nb = new NBall(nstart, m, d);
	nb->equalizationPeriod = 10;
	double init_pos[nstart][d], delta_pos[nstart][d];
	
	for(long int i = 0; i < nstart; i++){
		
		double* nodePos = nb->getNode(i)->position;
		
		for(long int j = 0; j < d; j++){
			
			init_pos[i][j] = nodePos[j];
			
		}
		
	}
	
	nb->grow(nend - nstart);
	
	for(long int i = 0; i < nstart; i++){
		
		double* nodePos = nb->getNode(i)->position;
		
		for(long int j = 0; j < d; j++){
			
			delta_pos[i][j] = nodePos[j] - init_pos[i][j];
			
		}
		
	}
	
	
	for(long int i = 0; i < nstart; i++){
		
		for(long int j = 0; j < d; j++){
			
			cout<<init_pos[i][j]<<" ";
			
		}
		
		for(long int j = 0; j < d; j++){
			
			cout<<delta_pos[i][j]<<" ";
			
		}
		
		cout<<endl;
		
	}
	
}
