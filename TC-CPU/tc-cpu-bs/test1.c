#include "iostream"
using namespace std;
typedef int vertex_t;
typedef int index_t;
index_t bsp(vertex_t *a, index_t m, vertex_t x){
	 
	vertex_t *base =a;
	index_t n=m;
	while (n > 1) {
		const index_t half = n / 2;
		base = (base[half] < x) ? &base[half] : base;
		n -= half;
//		cout<<n<<" "<<base[0]<<endl;
	}
        //return value is the first position larger than or equal to x.
	return (*base < x) + base - a;
}

int main(){
	int cache[33];
	int part[33];
	int b[128];
//	n=128;

	for(int i=0;i<128;i++){
		b[i]=i;
	}
	for(int i=0;i<33;i++){
		cache[i]=4*i;
	}
	part[0]=0;
	part[32]=128;
	int step = 16;
	while(step>0){
		cout<<"step="<<step;
		int j=0;
		int k=step;
		while(k<32){
//			part[k]=(part[k-step]+part[k+step])/2;
			part[k]=bsp(&b[part[k-step]],part[k+step]-part[k-step],cache[k])+part[k-step];
			cout<<" k="<<k;
			cout<<" part[k]="<<part[k];
			j++;
			k=(2*j+1)*step;
		}
		step /=2;
		cout<<endl;
	}


/*
	int mycount=0;
	int m;
	cin>>m;
//	cout<<endl;
	
	int *a = new int [m];
	for(int i=0;i<m;i++){
		a[i]=2*i+1;
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
*/

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


index_t intersec1(vertex_t *a, vertex_t *b, index_t m, index_t n){
	index_t mycount=0;
	index_t bot;
	index_t top;
	index_t r;

	for(size_t j=0; j<n; j++){
		vertex_t X = b[j];
		vertex_t Y;
		bot = 0;
		top = m-1;
		
		while(top>=bot){
			r = (top+bot)/2;
			Y = a[r];
			if(X==Y){
				mycount++;
				break;
			}
			else if(X<Y){
				top = r-1;
			}
			else if(X>Y){
				bot = r+1;
			}
		}
	}
	return mycount;
}
