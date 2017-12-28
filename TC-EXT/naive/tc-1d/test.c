#include "iostream"
#include "omp.h"
using namespace std;


int main(){
//#progma omp parallel num_threads(56)
#pragma omp parallel num_threads(56)
	{
		int ID = omp_get_thread_num();
		cout<<"hello "<<ID<<endl;
	}	
	return 0;
}

