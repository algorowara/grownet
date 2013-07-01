#include "../ngraph/nball.h"
#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>

using namespace std;

int main(){
	
	long int nmin = 16, nmax = 1024, nfactor = 2;
	long int d = 3, m = d+1;
	
	for(long int n = nmin; n <= nmax; n *= nfactor){
		
		NBall* ball = new NBall(n, m, d, NBALL_DEFAULT_RADIUS, NBALL_DEFAULT_ALPHA, NBALL_DEFAULT_GAMMA, 0.01, NBALL_DEFAULT_ITERATIONS);
		cout<<n<<" "<<ball->temp<<endl;
		delete ball;
		
	}
	
	for(long int n = nmin; n <= nmax; n *= nfactor){
		
		NBall* ball = new NBall(n, m, d, NBALL_DEFAULT_RADIUS, NBALL_DEFAULT_ALPHA, NBALL_DEFAULT_GAMMA, 0.0001, NBALL_DEFAULT_ITERATIONS);
		cout<<n<<" "<<ball->temp<<endl;
		delete ball;
		
	}
	
}
