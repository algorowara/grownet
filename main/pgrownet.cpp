#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>

using namespace std;
int a = 4; //1 if we want to grow whole network, 2 if we want step by step, 3 if we want ClustCoeff 
int dcare = 1; //1 if we want distance information, else 0

int main(){
  if(a == 1){	
	//make the network
	long int n = 5000, m = 3, nodeage;
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.05,.001,36);
	
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
	if(dcare == 1){ //if we want this information

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
			
			//realdist = net->linearDistance(net->getNode(i), net->getNode(j)); //the physical distance
			realdist = DISTANCE_2D(net->getNode(i), net->getNode(j));
			netdist = net->getNode(j)->distanceFromInitial; //the shortest path
			
			pdistance<<netdist<<" "<<realdist<<endl;			
	
		}
		
		net->clean(net->getNode(i));

	}	
	
	pdistance.close();

	}
	//get the edge age vs. edge betweenness information
	
	ofstream eagebetw;
	eagebetw.open("eagebetw.txt", ios::out | ios::trunc);

	double* edgeinfo = net->edgeAgeVsBetweenness();

	for(long int i = 0; i < n; i++){

		eagebetw<<i<<" "<<edgeinfo[i]<<endl;

	}

	eagebetw.close();		
		
  }
  //instead grow the network node at a time and get the clustering coefficient and characteristic path length at each timestep
  else if (a == 2) { 
	long int n = 1000, m = 3;
	double cc = 0; //clustering coefficient
	double cpl = 0; //characteristic path length
	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(4,m,.05,.0001,36);
	
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
  else if (a == 3) {
	long int m = 10;
	double cc = 0, cpl = 0;; //clustering coefficient and characteristic path length
	ofstream pclust2;
	pclust2.open("pclust2.txt", ios::out | ios::trunc);		
	ofstream ppath2;
	ppath2.open("ppath2.txt", ios::out | ios::trunc);

	for(long int i = 50; i < 1001; i += 50){ //step from N = 50 to N = 1000

		for(long int j = 0; j < 4; j++){ //four data points at each size	

			PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(i,m,.05,.001,36);
			cc = net->averageClusteringCoefficient();
			pclust2<<i<<" "<<cc<<endl;		

			cpl = net->averagePathLength();
			ppath2<<i<<" "<<cpl<<endl;
	
		}

	}

	pclust2.close();
	ppath2.close();

  }
  else if (a == 4) {	//node betweenness info & tests
	long int m = 3, n = 20;
	double* nodeBetw;	

	PositiveChargeGrowingNetwork2D* net = new PositiveChargeGrowingNetwork2D(n,m,.05,.001,36);

	nodeBetw = net->nodeBetweenness();	//test the node betweenness method

	ofstream pnodebetw;
	pnodebetw.open("pnodebetw.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){
	
		pnodebetw<<net->K(i)<<" "<<nodeBetw[i]<<endl;

	}

	pnodebetw.close();

  }	

}
