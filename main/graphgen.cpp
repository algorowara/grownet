#include "../ngraph/nsphere.h"
#include "../ngraph/nball.h"
#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <cstring>
#include <sys/time.h>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int n = 100, m = 3, d = 1;
	NSphere* net = new NSphere(n, m, d);
	
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
