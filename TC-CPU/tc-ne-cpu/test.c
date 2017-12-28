#include "iostream"
using namespace std;


int main(){
	int mycount=0;
	int m;
	cin>>m;
//	cout<<endl;
	
	int *a = new int [m];
	for(int i=0;i<m;i++){
		a[i]=2*i;
	}
	int *base =a;
	int x;
	cin>>x;
//	cout<<endl;
	int n=m;
	while (n > 1) {
		const int half = n / 2;
		base = (base[half] < x) ? &base[half] : base;
		n -= half;
		cout<<n<<" "<<base[0]<<endl;
	}
	if(a[(*base < x) + base - a]==x && base<(a+m-1)){
		mycount++;
	}
	cout<<"found "<<mycount<<endl;
//	cout<<base<<" "<<a<<" "<<a+m<<endl;
	cout<<"return "<<(*base < x) + base - a<<" a[] = "<<a[(*base < x) + base - a]<<endl;


/*	
	cout<<endl;
	int power = 1;
	int temp = m;
	while(temp){
		temp = temp/2;
		power= power*2;
	}
	for(int i=0; i<m; i++){
		int p = 0;
		int temp = i+1;
		while(temp){
			temp = temp/2;
			p++;
		}
		int q=1;
		for(int j=1; j<p; j++){
			q = q*2;
		}
		
		int index = power * (2*(i-q+2)-1) /(2*q)-1;
		cout<<i<<" "<<index<<endl;
	}
*/
	return 0;

}

