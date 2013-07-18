#include "../pp1d/pp1d.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include "../ngraph/ngraph.h"
#include "../graph/growingnetwork.h"
#include "../graph/graph.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
int fignum = 12;


int main(){

	if(fignum == 0){	// generate graphs to be used later
		
		long int sampsize = 10;	// the number of graphs of each dimension and topology to be generated
		long int nmin = 500, nmax = 10000, nstep = 500;	// the starting n, ending n, and step size; not independent
		long int m = 4;	// the number of connections formed by each new node
		long int dmin = 1, dmax = 4, dstep = 1;	// the number of dimensions to go through
		
		for(long int d = dmin; d <= dmax; d++){	// for each dimension to be reviewed
			
			for(long int s = 0; s < sampsize; s++){	// generate some number of samples
				
				NBall* net = NULL;
				
				for(long int n = nmin; n <= nmax; n += nstep){	// grow the network at each step
																// and export it to prove n-dependent or asymptotic properties
					
					stringstream filename;
					filename<<"figuredata_"<<s<<"_"<<d<<"ball_m"<<m<<"_n"<<n<<".nball";
					// filename = figuredata_SAMPLENUM_Dball_mM_nN.nball
					
					if(net == NULL){	// if the network does not yet exist
						
						net = new NBall(n, m, d);	// create it
						
					}
					
					else{	// otherwise
						
						net->grow(nstep);	// grow it
						
					}
					
					NBall::exportObject(net, filename.str().c_str());	// export it
					
				}	// end growth of this particular instance
				
				delete net;
				
			}	// end NBall sample generation for this dimension
			
			for(long int s = 0; s < sampsize; s++){	// do the same thing for the NSphere
				
				NSphere* net = NULL;
				
				for(long int n = nmin; n <= nmax; n += nstep){	// grow the network in steps
					
					stringstream filename;
					filename<<"figuredata_"<<s<<"_"<<d<<"sphere_m"<<m<<"_n"<<n<<".nsphere";
					// filename = figuredata_SAMPLENUM_Dball_mM_nN.nsphere
					
					if(net == NULL){	// if the network does not yet exist
						
						net = new NSphere(n, m, d);	// create it
						
					}
					
					else{	// otherwise
						
						net->grow(nstep);	// grow it
						
					}
					
					NSphere::exportObject(net, filename.str().c_str());	// export it
					
				}	// end growth of this particular instance
				
				delete net;
				
			}	// end NSphere sample creation for this dimension
			
		}	// end all sample generation for this dimension
		
	}	// end sample generation
	
	if(fignum == 11){	//degree distribution

		long int n = 10000, m = 4;

		GrowingNetwork2D* T1D = new GrowingNetwork2D(n,m);
		PPGrowingNetwork1D* PP1D = new PPGrowingNetwork1D(n,m,1,.00001,64);
		NSphere* T2D = NSphere::importObject("figuredata_1_2sphere_m4_n10000.nsphere");
		NBall* PP2D = NBall::importObject("figuredata_1_2ball_m4_n10000.nball");
		NSphere* T3D = NSphere::importObject("figuredata_1_3sphere_m4_n10000.nsphere");
		NBall* PP3D = NBall::importObject("figuredata_1_3ball_m4_n10000.nball");
		NSphere* T4D = NSphere::importObject("figuredata_1_4sphere_m4_n10000.nsphere");
		NBall* PP4D = NBall::importObject("figuredata_1_4ball_m4_n10000.nball");
		
		ofstream degreedist;
		degreedist.open("degreedist.txt", ios::out | ios::trunc);
		
		float* T1Ddist = T1D->degreeDistribution();
		float* PP1Ddist = PP1D->degreeDistribution();
		float* T2Ddist = T2D->degreeDistribution();
		float* PP2Ddist = PP2D->degreeDistribution();
		float* T3Ddist = T3D->degreeDistribution();
		float* PP3Ddist = PP3D->degreeDistribution();
		float* T4Ddist = T4D->degreeDistribution();
		float* PP4Ddist = PP4D->degreeDistribution();
		
		for(long int i = 0; i < n; i++){
			
			//the first column is the degree, and then we have:
			//2nd column = T1D
			//3rd column = PP1D
			//4th column = T2D
			//5th column = PP2D
			//6th column = T3D
			//7th column = PP3D
			//8th column = T4D
			//9th column = PP4D
			
			degreedist<<i<<", ";
			degreedist<<T1Ddist[i]<<", ";	
			degreedist<<PP1Ddist[i]<<", ";
			degreedist<<T2Ddist[i]<<", ";
			degreedist<<PP2Ddist[i]<<", ";
			degreedist<<T3Ddist[i]<<", ";
			degreedist<<PP3Ddist[i]<<", ";
			degreedist<<T4Ddist[i]<<", ";
			degreedist<<PP4Ddist[i]<<endl;
		
		}
		
		degreedist.close();
		
	}

	else if (fignum == 12){	//degree vs node age
		
		long int n = 10000, m = 4;
		long int age, pdegree, tdegree;
		
		NSphere* T2D = NSphere::importObject("figuredata_1_2sphere_m4_n10000.nsphere");
		NBall* PP2D = NBall::importObject("figuredata_1_2ball_m4_n10000.nball");
		
		ofstream degvsage;
		degvsage.open("degvsage.txt", ios::out | ios::trunc);
		
		for(long int i = 0; i < n; i++){
			
			age = (T2D->getTime() - T2D->getNode(i)->getStartTime());
			
			tdegree = T2D->K(i);
			pdegree = PP2D->K(i);
			
			//the first column is the age
			//the 2nd column is the degree from the T2D model
			//the 3rd column is the degree from the PP2D model
			degvsage<<age<<", "<<tdegree<<", "<<pdegree<<endl;
			
		}
		
		degvsage.close();
		
	}
	
}



	
