#include "../ngraph/nball.h"
#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace std;

int main(){
	
	long int n = 100, m = 0;
	long int samp = 10000000, num_boxes = 1000;
	double avg_force[num_boxes];
	long int count[num_boxes];
	
	for(long int d = 2; d < 6; d++){
		
		stringstream filename;
		filename<<"f_vs_r_"<<d<<"d_n"<<n<<".txt";
		ofstream filestream;
		filestream.open(filename.str().c_str(), ios::out | ios::trunc);
		
		memset(&avg_force, 0, num_boxes * sizeof(double));
		memset(&count, 0, num_boxes * sizeof(long int));	
		double s = samp;
		
		NBall* b = new NBall(n, m, d);
		
		while(s > 0){
			
			SpatialVertex* test = new SpatialVertex(d, b->randomLocation(), -1);
			double* force = b->sumForces(test);
			double radius = test->radialDistance();
			double dot = 0;
			
			for(long int i = 0; i < d; i++){
				
				dot += force[i] * test->position[i] / radius;
				
			}
			
			if(abs(dot) < 2.0){
			
				avg_force[(long int)(radius * num_boxes)] += dot;
				count[(long int)(radius * num_boxes)]++;
				s--;
				
			}
			
			delete test;
			delete force;
			
		}
		
		for(long int i = 0; i < num_boxes; i++){
			
			if(count[i] == 0){
				
				continue;
				
			}
			
			filestream<<((2.0 * i + 1)/(2.0 * num_boxes))<<" "<<(avg_force[i]/count[i])<<endl;
		
		}
		
		filestream.close();
		delete b;
		
	}
	
}
