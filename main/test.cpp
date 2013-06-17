#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	long int n = 500, m = 3;
	double gamma = 1.0, tol = 0.001;
	long int maxItr = 1000;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tol, maxItr);
	double init[n][DIM], final[n][DIM];
	ofstream displacement, end_position, random_position;
	
	/* begin displacement measurements */
	displacement.open("displacement.txt", ios::out | ios::trunc);
	for(int k = n; k < n + 20; k++){
	
		for(int i = 0; i < n; i++){
			
			for(int j = 0; j < DIM; j++){
				
				init[i][j] = net->nodes.at(i)->position[j];
			
			}
			
		}
		
		net->grow(1);
		
		for(int i = 0; i < n; i++){
			
			for(int j = 0; j < DIM; j++){
				
				final[i][j] = net->nodes.at(i)->position[j];
			
			}
			
		}
		
		for(int i = 0; i < n; i++){
			
			double dx = final[i][0] - init[i][0];
			double dy = final[i][1] - init[i][1];
			double dz = final[i][2] - init[i][2];
			double dr = sqrt(dx * dx + dy * dy + dz * dz);
			
			displacement<<dr<<endl;
			
		}
		
	}
	
	delete net;
	displacement.close();
	
	/* end displacement measurements */
	
	/* begin sphere measurements */
	
	n = 1000;
	net = new GrowingNetwork3D(n, m, gamma, tol, maxItr);
	
	end_position.open("end_position.txt", ios::out | ios::trunc);
	
	for(int i = 0; i < n; i++){
		
		for(int j = 0; j < DIM; j++){
			
			end_position<<net->nodes.at(i)->position[j];
			
			if(j < DIM-1){
				
				end_position<<" ";
				
			}
			
		}
		
		end_position<<endl;
		
	}
	
	delete net;
	end_position.close();
	
	/* end sphere measurements */
	
	/* begin random measurements */
	
	random_position.open("random_position.txt", ios::out | ios::trunc);
	
	for(int i = 0; i < n; i++){
		
		double* rand_pos = net->randomLocation();
		
		for(int j = 0; j < DIM; j++){
			
			random_position<<rand_pos[j];
			
			if(j < DIM-1){
			
				random_position<<" ";
				
			}
			
		}
		
		random_position<<endl;
		delete rand_pos;
		
	}
	
	random_position.close();
	
	/* end random measurements */
	
}
