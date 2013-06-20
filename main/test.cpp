#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 1000, n_max = 1050, m = 3, dim = 3;
	long int itr = 1000;
	long int sample_size = 20;
	long int num_boxes = 1000;
	double dr_box[num_boxes];	// histogram "boxes" from 0 to 0.1, in increments of 0.0001
							// range[i] = dr_max * i/1000 to dr_max * (i+1)/1000
							// index[value] = (long int)(1000 * value/dr_max)
	double dr_max = 0.005;
	long int num_data = 0;
	
	for(long int i = 0; i < num_boxes; i++){
		
		dr_box[i] = 0;
		
	}
	
	for(long int s = 0; s < sample_size; s++){
		
		GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
		
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
						
						dr_box[(long int)(num_boxes * dr/(dr_max - 0.000000001))]++;
						num_data++;
						
					}
						
				}
				
			}
			
			delete net;
			
		}
	
	for(long int i = 0; i < num_boxes; i++){
		
		cout<<(dr_max * (2.0 * i + 1))/(2.0 * num_boxes)<<" "<<dr_box[i]<<endl;
		
	}
			
}
