#include "../ngraph/nball.h"
#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace std;

int main(){
	
	long int nmin = 10, nmax = 6000, nfactor = 2;
	long int d = 3, m = d+1;
	
	ofstream lfile, cfile, kfile;
	stringstream lname, cname, kname;
	
	lname<<"l_vs_n_"<<d<<"ball_n"<<nmin<<"_to_n"<<nmax<<".txt";
	cname<<"c_vs_n_"<<d<<"ball_n"<<nmin<<"_to_n"<<nmax<<".txt";
	
	lfile.open(lname.str().c_str(), ios::out | ios::trunc);
	cfile.open(cname.str().c_str(), ios::out | ios::trunc);
	
	NBall* nb = new NBall(nmin, m, d);
	
	for(long int n = nmin; n < nmax; n *= nfactor){
		
		lfile<<n<<" "<<nb->averagePathLength()<<endl;
		cfile<<n<<" "<<nb->averageClusteringCoefficient()<<endl;
		
		nb->grow(n * (nfactor -1));
		
	}
	
	lfile.close();
	cfile.close();
	
	kname<<"k_dist_"<<d<<"ball_n"<<nmax<<".txt";
	kfile.open(kname.str().c_str(), ios::out | ios::trunc);
	double* deg = nb->degreeDistribution();
	
	for(long int i = 0; i < nmax; i++){
		
		if(deg[i] > 0){
			
			kfile<<i<<" "<<deg[i]<<endl;
			
		}
		
	}
	
	kfile.close();
	delete deg;
	delete nb;
	
}
