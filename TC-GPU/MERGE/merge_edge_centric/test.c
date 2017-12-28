#include "iostream"
using namespace std;

int BinarySearch(int x, int*A, int bot, int top){
	int r= (bot+top)/2;
//	int result;
	while(top>bot){
		if(x<A[r]){
			top = r;
		}
		else if(x>A[r]){
			bot = r+1;
		}
		else if(x==A[r]){
			break;
		}
		r = (bot+top)/2;
	}
	return r;
}

int main(){
	//phase 1, partition
	int a[8]={17,29,64,73,86,90,95,99};
	for(int i=0;i<8;i++){
		cout<<a[i]<<" ";
	}
	cout<<"\n";

	int *A[4];
	A[2]=a;
	for(int i=0;i<8;i++){
		cout<<A[2][i]<<" ";
	}
	cout<<"\n";

	int x = 50;
	int position = BinarySearch(x,a,0,7);

	cout<<"x= "<<x<<"  position= "<<position<<"  value= "<<a[position]<<"\n";
	return 0;

}

