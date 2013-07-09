#include "../pp1d/pp1d.h"
#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <fstream>

using namespace std;
long int a = 4;

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
	long int m = 5, n = 1000, age, nodeage;
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
  else if(a == 2){	//7 & 8, age and path and distance
 	long int n = 1000, m = 3, netdist, agediff;
        double realdist;
        PPGrowingNetwork1D* net = new PPGrowingNetwork1D(n, m, 1, .00001, 64);

        ofstream ldistance;
        ldistance.open("ldistance.txt", ios::out | ios::trunc);
        ofstream lagedistance;
        lagedistance.open("lagedistance.txt", ios::out | ios::trunc);
        ofstream lagepath;
        lagepath.open("lagepath.txt", ios::out | ios::trunc);

        for(long int i = 0; i < n; i++){ //for all the nodes in the network

                net->memoize(net->getNode(i));

                for(long int j = 0; j < n; j++){ //for all of the other nodes

                        if(i == j){
                                continue;
                        }

                        //realdist = net->linearDistance(net->getNode(i), net->getNode(j)); //the physical distance
                        realdist = DISTANCE_1D(net->getNode(i), net->getNode(j));
                        netdist = net->getNode(j)->distanceFromInitial; //the shortest path
                        agediff = abs((net->getNode(i)->getStartTime()) - (net->getNode(j)->getStartTime()));

                        ldistance<<netdist<<" "<<realdist<<endl;
                        lagedistance<<agediff<<" "<<realdist<<endl;
                        lagepath<<agediff<<" "<<netdist<<endl;

                }

                net->clean();

        }

        ldistance.close();
        lagedistance.close();
        lagepath.close();

  }
  else if(a == 3){	//clustering and shortest path length
 	long int m = 3;
        double cc = 0, cpl = 0;; //clustering coefficient and characteristic path length
        ofstream lclust2;
        lclust2.open("lclust2.txt", ios::out | ios::trunc);
        ofstream lpath2;
        lpath2.open("lpath2.txt", ios::out | ios::trunc);

        for(long int i = 10; i < 1001; i += 10){ //step from N = 50 to N = 1000

                for(long int j = 0; j < 5; j++){ //four data points at each size        

                        GrowingNetwork2D* net = new GrowingNetwork2D(i,m);
                        cc = net->averageClusteringCoefficient();
                        lclust2<<i<<" "<<cc<<endl;

                        cpl = net->averagePathLength();
                        lpath2<<i<<" "<<cpl<<endl;

                }

        }

        lclust2.close();
        lpath2.close();

  }
  else if(a == 4){
	long int n = 1000, m = 3;
	bool edge;
	double realdist;
	PPGrowingNetwork1D* net = new PPGrowingNetwork1D(n, m,1,.00001,64);
	
	ofstream ledge;
	ledge.open("ledge.txt", ios::out | ios::trunc);

	for(long int i = 0; i < n; i++){

		for(long int j = 0; j < n; j++){

			edge = 0;

			if(i == j){
				continue;
			}

			realdist = DISTANCE_1D(net->getNode(i), net->getNode(j));
			if(net->getNode(i)->hasNeighbor(net->getNode(j))){
				edge = 1;

			}

			ledge<<realdist<<" "<<edge<<endl;

		}

	}

	ledge.close();

  }

}
