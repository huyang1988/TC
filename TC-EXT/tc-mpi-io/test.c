#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> 
#include <stdint.h>
#include <sys/stat.h>

using namespace std;



int main(int argc, char **argv) {

	
	int i;
	int j=0;
	for(int t=0; t<55; t++){
		i = (int)(-0.5 + sqrt(2*t+1));
		j = t - i*(i+1)/2;
		cout<<i<<"-"<<j<<endl;
	}	


	
	
	
	
	return 0;
}

