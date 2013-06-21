#include "nball.h"

NBall::NBall(long int d, long int n, long int m, double r, double a, double b, double g, double t, long int i) : DIM(d) {
	
	static bool randSeeded = false;	// static variable to check if the pseudorandom number generator has been seeded yet
	
	if(!randSeeded){	// if it has not
		
		srand(time(NULL));	// do so
		randSeeded = true;	// and reflect this by setting the variable to true
		
	}
	
	// initialize the non-constant variables
	this->m = m;
	this->radius = r;
	this->baseAlpha = a;
	this->baseBeta = b;
	this->baseGamma = g;
	this->baseTolerance = t;
	this->baseItr = i;
	
}
