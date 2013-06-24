#include "../growingnetwork3d/growingnetwork3d.h"
#include "../ngraph/nball.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int d = 3, n = 10, m = 3;
	NBall* b = new NBall(d, n, m);
	delete b;
			
}
