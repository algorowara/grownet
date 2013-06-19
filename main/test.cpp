#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 1000, n_max = 1050, m = 3, dim = 3;
	double gamma = 1.0, tolerance = 0.001;
	long int itr = 1000;
	long int sample_size = 20;
	long int num_boxes = 1000;
	double dr_box[num_boxes];	// histogram "boxes" from 0 to 1.0, in increments of 0.001
							// range[i] = i/1000 to (i+1)/1000
							// index[value] = (long int)(1000 * value)
	double dr_max = 1.0;
	
	for(long int i = 0; i < num_boxes; i++){
		
		dr_box[i] = 0;
		
	}
	
	for(long int s = 0; s < sample_size; s++){
		
		GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tolerance, itr);
		
		for(int i = n; i < n_max; i++){
			
			double init[n][dim];
			
			for(long int j = 0; j < n; j++){
				
				for(long int k = 0; k < dim; k++){
					
					init[j][k] = net->getNode(j)->position[k];
					
				}
				
			}
			
			net->grow(1);
			
			for(long int j = 0; j < n; j++){
				
				double dx = net->getNode(j)->position[0] - init[j][0];
				double dy = net->getNode(j)->position[1] - init[j][1];
				double dz = net->getNode(j)->position[2] - init[j][2];
				double dr = sqrt(dx * dx + dy * dy + dz * dz);
				
					if(dr < dr_max){
						
						dr_box[(long int)(num_boxes * dr)]++;
						
					}
						
				}
				
			}
			
			delete net;
			
		}
	
	for(long int i = 0; i < num_boxes; i++){
		
		cout<<(2 * i + 1)/((double)2 * num_boxes)<<" "<<dr_box[i]<<endl;
		
	}
			
}
