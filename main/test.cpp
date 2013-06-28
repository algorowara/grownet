#include "../ngraph/nball.h"
#include "../pgrownet2d/pgrownet2d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace std;

int main(){
	
	long int n = 10, m = 3, d = 2;
	NBall* b = new NBall(n, m, d);
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = 0; j < d; j++){
			
			cout<<b->getNode(i)->position[j];
			
			if(j < d-1){
				
				cout<<" ";
			
			}
			
		}
		
		cout<<endl;
		
	}
	
}
