#include "../growingnetwork3d/growingnetwork3d.h"
#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>

using namespace std;

int main(){
	
	GrowingNetwork3D* net = new GrowingNetwork3D(100, 3);
	long int d = net->DIM;
	
	cout<<(net->DIM)<<endl;
	
}
