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
	
	long int n = 1000, m = 3, d = 2;
	NBall* net = new NBall(n, m, d);
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = 0; j < d+1; j++){
			
			cout<<net->getNode(i)->position[j];
			
			if(j < d){
				
				cout<<" ";
				
			}
			
		}
		
		cout<<endl;
		
	}
	
}
