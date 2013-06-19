#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;
double a = 1;

int main(){
  if(a == 1){	
	//make the network
	long int n = 100, m = 3, nodeage = 0;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.05,.0001,10000);
	
	//save the position of the points as validation
	//with this we also want to look at the node ages
	ofstream pgrowdata;
	pgrowdata.open("pgrowdata.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		for(long int j = 0; j < DIM; j++){

			pgrowdata<<net->getNode(i)->position[j];

			if(j < DIM-1){

				pgrowdata<<" ";

			}

		}
		
		nodeage = (net->getTime() - net->getNode(i)->getStartTime()); //get the age of the node (current time - node birthday)
		pgrowdata<<" "<<nodeage<<endl;

	}

	pgrowdata.close();
	cout << net->calculatePotential() << endl;
	
	//get the degree distribution
	double* dist = net->degreeDistribution();

	ofstream pdegdist;
	pdegdist.open("pdegdist.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){
	
		pdegdist<<i<<" "<<dist[i]<<endl;
	
	}
	
	pdegdist.close();

	//get the degree vs. age distribution
	long int age, degree;
	ofstream pagedeg;
	pagedeg.open("pagedeg.txt", ios::out | ios::trunc);	

	for(long int i = 0; i < n; i++){
		
		age = (net->getTime() - net->getNode(i)->getStartTime());
		degree = net->getNode(i)->neighbors.size();

		pagedeg<<age<<" "<<degree<<endl;

	}		
	
	//get the real distance vs. network distance information
	long int netdist;
	double realdist;
	ofstream pdistance;
	pdistance.open("pdistance.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){ //for all the nodes in the network
		
		net->memoize(net->getNode(i));

		for(long int j = 0; j < n; j++){ //for all of the other nodes

			if(i == j){
				continue;
			}
			
			realdist = net->linearDistance(net->getNode(i), net->getNode(j)); //the physical distance
			netdist = net->getNode(j)->distanceFromInitial; //the shortest path
			
			pdistance<<netdist<<" "<<realdist<<endl;			
	
		}
		
		net->clean(net->getNode(i));

	}	
	
	pdistance.close();
		
  }
  //instead grow the network node at a time and get the clustering coefficient and characteristic path length at each timestep
  else { 
	long int n = 1000, m = 3;
	double cc = 0; //clustering coefficient
	double cpl = 0; //characteristic path length
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(4,m,.05,.0001,10000);
	
	ofstream pclust;
	pclust.open("pclust.txt", ios::out | ios::trunc);
	ofstream ppath;
	ppath.open("ppath.txt", ios::out | ios::trunc);	

	for(long int i = 4; i < n; i++){
		
		cc = net->averageClusteringCoefficient();
		pclust<<i<<" "<<cc<<endl;		

		cpl = net->averagePathLength();
		ppath<<i<<" "<<cpl<<endl;	
		
		net->grow(1);

	}
	
	pclust.close();
	ppath.close();

  }
}
