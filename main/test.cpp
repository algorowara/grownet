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
	
	long int n = 1000, m = 3, d = 2;
	NGraph* net = new NBall(m+1, m, d);
	net->baseGam = 10;
	net->grow(n - (m+1));
	cout<<net->clustering();
	
}
