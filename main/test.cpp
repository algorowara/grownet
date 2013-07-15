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
	
	long int n = 1000, m = 11, d = 10;
	NBall* net = new NBall(m+1, m, d);
	net->baseGam = 4.0;
	net->grow(n - (m+1));
	
}
