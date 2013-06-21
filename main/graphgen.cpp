#include "../ngraph/nsphere.h"
#include <iostream>
#include <fstream>
#include <ctime>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int d = 3, m = 1;
	NSphere* s;
	GrowingNetwork3D* net = new GrowingNetwork3D(4, 3);
	
	for(long int n = 2; n < 100; n++){
		
		s = new NSphere(d, n, m);
		cout<<n<<" "<<(s->calculatePotential()/net->calculateMinimumPotential(n, d))<<endl;
		delete s;
		
	}
	
}
