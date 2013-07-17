#include "../ngraph/ngraph.h"
#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>
#include <cfloat>
#include <algorithm>
#include <iterator>

using namespace std;

int main(){
	
	long int n = 100, m = 4, d = 2;
	NBall* net = new NBall(m+1, m, d);
	net->baseGam = 0.1;
	net->forceExp = 3;
	net->grow(n - (m+1));
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = 0; j < d; j++){
			
			cout<<net->getNode(i)->position[j]<<" ";
			
		}
		
		cout<<endl;
		
	}
	
}
