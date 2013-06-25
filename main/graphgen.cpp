#include "../ngraph/nball.h"
#include <iostream>
#include <fstream>
#include <ctime>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int n = 400, m = 0, d = 3;
	
	NBall* nb = new NBall(n, m, d);
	
	delete nb;
	
}
