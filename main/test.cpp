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
#include <climits>

using namespace std;

int main(){
	
	long int nstart = 10, nend = 100, m = 3, d = 2;
	NSphere* start = new NSphere(nstart, m, d);
	float* newLoc[nend];
	
	for(long int n = 0; n < nstart; n++){
		
		newLoc[n] = start->getNode(n)->position;
		
	}
	
	for(long int n = nstart; n < nend; n++){
		
		newLoc[n] = start->randomLocation();
		
	}
	
	for(long int f = 0; f <= 3; f++){
		
		for(float gamma = 0.001; gamma > 0.000005; gamma /= 10){
		
			for(long int num = 0; num < 1; num++){
				
				float loc[nend][d+1];
				//ofstream file;
				//stringstream name;
				//name<<"points_2sphere_g"<<gamma<<"_s"<<num<<".txt";
				//file.open(name.str().c_str(), ios::out | ios::trunc);
				
				NSphere* end = new NSphere(nstart, m, d);
				
				for(long int n = 0; n < nstart; n++){
					
					memcpy(loc[n], start->getNode(n)->position, (d+1) * sizeof(float));
					end->nodes.pop_back();
					
				}
				
				for(long int n = 0; n < nstart; n++){
					
					SpatialVertex* newNode = new SpatialVertex(d+1, loc[n], end->getTime());
					end->addNode(newNode);
					
				}
				
				for(long int n = 0; n < nstart; n++){
					
					for(long int o = n+1; o < nstart; o++){
						
						if(start->getNode(n)->hasNeighbor(start->getNode(o))){
							
							end->getNode(n)->addNeighbor(end->getNode(o));
							
						}
						
					}
					
				}
				
				end->baseGam = gamma;
				end->forceExp = f;
				end->baseItr = LONG_MAX;
				
				for(long int n = nstart; n < nend; n++){
					
					for(long int i = 0; i < d+1; i++){
						
						loc[n][i] = newLoc[n][i];
						
					}
					
				}
				
				for(long int n = nstart; n < nend; n++){
					
					SpatialVertex* newNode = new SpatialVertex(d+1, loc[n], end->getTime());
					SpatialVertex** nearNeighbors = end->findMNearestNeighbors(newNode);
		
					for(long int i = 0; i < m; i++){
						
						newNode->addNeighbor(nearNeighbors[i]);
						
					}
		
					end->addNode(newNode);
					end->equalize();
					
				}
				
				for(long int n = 0; n < nend; n++){
					
					for(long int i = 0; i < d+1; i++){
						
						//file<<end->getNode(n)->position[i]<<" ";
						
					}
					
					//file<<endl;
					
				}
				
				cout<<f<<" "<<gamma<<" "<<end->clustering()<<endl;
				//file.close();
				delete end;
				
			}
			
		}
		
	}
	
}
