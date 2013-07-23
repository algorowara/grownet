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
	
	for(long int f = 0; f <= 8; f++){
	
		long int n = 1000, m = 12, d = 4;
		NSphere* net = new NSphere(m+1, m, d);
		net->baseGam = 0.01;
		net->forceExp = f;
		net->grow(n - (m+1));
		cout<<f<<" "<<net->clustering()<<endl;
		delete net;
		
	}
	
}
