#include "../ngraph/nball.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace std;

int main(){
	
	long int d = 3, m = 3;
	long int sample = 10000;
	
	for(long int n = 50; n <= 500; n += 50){
		
		ofstream fstream;
		stringstream sstream;
		sstream<<"force_vs_radius_n"<<n<<".txt";
		string s = sstream.str();
		fstream.open(s.c_str(), ios::out | ios::trunc);
		
		for(long int i = 0; i < sample; i++){
			
			NBall* b = new NBall(n, m, d);
			SpatialVertex* test = new SpatialVertex(d, b->randomLocation(), -1);
			double* force = b->sumForces(test);
			double dot = 0;
			for(long int j = 0; j < d; j++){
				
				dot += force[j] * (test->position[j]) / b->radialDistance(test);
				
			}
			
			if(abs(dot) < 2){
				
				fstream<<b->radialDistance(test)<<" "<<dot<<endl;
			
			}
			
		}
		
		fstream.close();
		
	}
			
}
