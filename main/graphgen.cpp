#include "../ngraph/nsphere.h"
#include <iostream>
#include <fstream>
#include <ctime>

#define EPSILON pow(2.0, -16.0)

using namespace std;

int main(){
	
	long int n = 1000, m = 3;
	long int sample_size = 36;
	double edgeAgeVsBetweenness[n];
	double edgeAgeVsDistance[n];
	double nodeAgeVsDegree[n];
	ofstream between, distance, degree;
	
	between.open("edge_age_vs_betweenness_n1000_m3.txt", ios::out | ios::trunc);
	distance.open("edge_age_vs_distance_n1000_m3.txt", ios::out | ios::trunc);
	degree.open("node_age_vs_degree_n1000_m3.txt", ios::out | ios::trunc);
	
	for(long int i = 0; i < n; i++){
		
		edgeAgeVsBetweenness[i] = 0;
		edgeAgeVsDistance[i] = 0;
		nodeAgeVsDegree[i] = 0;
		
	}
	
	for(long int i = 0; i < sample_size; i++){
		
		GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
		
		double* localBetweenness = net->edgeAgeVsBetweenness();
		double* localDistance = net->edgeAgeVsLinearDistance();
		double* localDegree = net->nodeAgeVsDegree();
		
		for(long int j = 0; j < n; j++){
			
			edgeAgeVsBetweenness[j] += localBetweenness[j];
			edgeAgeVsDistance[j] += localDistance[j];
			nodeAgeVsDegree[j] += localDegree[j];
			
		}
		
		delete localBetweenness;
		delete localDistance;
		delete localDegree;
		delete net;
		
	}
	
	for(long int i = 0; i < n; i++){
		
		edgeAgeVsBetweenness[i] /= sample_size;
		edgeAgeVsDistance[i] /= sample_size;
		nodeAgeVsDegree[i] /= sample_size;
		
	}
	
	for(long int i = 0; i < n; i++){
		
		between<<i<<" "<<edgeAgeVsBetweenness[i]<<endl;
		distance<<i<<" "<<edgeAgeVsDistance[i]<<endl;
		degree<<i<<" "<<nodeAgeVsDegree[i]<<endl;
		
	}
	
	between.close();
	distance.close();
	degree.close();
	
}
