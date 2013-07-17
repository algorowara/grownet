#include "../pp1d/pp1d.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include "../ngraph/ngraph.h"
#include "../graph/growingnetwork.h"
#include "../graph/graph.h"
#include <iostream>
#include <fstream>

using namespace std;
int fignum = 12;


int main(){
  if(fignum == 11){	//degree distribution
	long int n = 10000, m = 4;

	GrowingNetwork2D* T1D = new GrowingNetwork2D(n,m);
	PPGrowingNetwork1D* PP1D = new PPGrowingNetwork1D(n,m,1,.00001,64);
	NSphere* T2D = new NSphere(n,m,2);
	NBall* PP2D = new NBall(n,m,2);
	NSphere* T3D = new NSphere(n,m,3);
	NBall* PP3D = new NBall(n,m,3);
	NSphere* T4D = new NSphere(n,m,4);
	NBall* PP4D = new NBall(n,m,4);

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

	NSphere* T2D = new NSphere(n,m,2);
	NBall* PP2D = new NBall(n,m,2);

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



	
