#include "../pp1d/pp1d.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>

using namespace std;
long int a = 1;

int main(){
  if(a == 0){

	long int n = 1000, m = 2;
	PPGrowingNetwork1D* net = new PPGrowingNetwork1D(n,m,1, .00001, 64);

	ofstream linedata;
	linedata.open("linedata.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		linedata<<net->getNode(i)->position[0]<<endl;

	}
	
	cout<<net->calculatePotential()<<endl;

  }
  else if(a == 1){	//plots 1, 2, 3, 4 TP1D
	long int m = 9, n = 1000, age, nodeage;
	double* nodeBetw;

	GrowingNetwork2D* net = new GrowingNetwork2D(n, m);

	//get the degree v betweenness and age v betweenness
	nodeBetw = net->nodeBetweenness();

	ofstream ldegbetw;
	ldegbetw.open("ldegbetw.txt", ios::out | ios::trunc);
	ofstream lagebetw;
	lagebetw.open("lagebetw.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		age = (net->getTime() - net->getNode(i)->getStartTime());

		ldegbetw<<net->K(i)<<" "<<nodeBetw[i]<<endl;	//node degree vs betweenness
		lagebetw<<age<<" "<<nodeBetw[i]<<endl;	//node age vs betweenness

	}

	ldegbetw.close();
	lagebetw.close();

	//get the degree distribution
	double* dist = net->degreeDistribution();

	ofstream ldegdist;
	ldegdist.open("ldegdist.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		ldegdist<<i<<" "<<dist[i]<<endl;

	}

	ldegdist.close();

	//get the age vs degree distribution
	long int degree;
	ofstream lagedeg;
	lagedeg.open("lagedeg.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		age = (net->getTime() - net->getNode(i)->getStartTime());
		degree = net->getNode(i)->neighbors.size();

		lagedeg<<age<<" "<<degree<<endl;

	}

	lagedeg.close();

  }
	
}
