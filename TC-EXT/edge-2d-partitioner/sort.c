//sort.cu
#include "comm.h"
#include "iostream"
using namespace std;

void printGraph(vertex_t vertCount, 
		vertex_t* head, 
		vertex_t* adj, 
		index_t* begin){
	for(vertex_t i=0; i<vertCount; i++){
		if(begin[i+1]>begin[i]){
			cout<<begin[i]<<" "<<begin[i+1]-begin[i]<<": ";
		}
//		for(int j=0; j<degree[i]; j++){
		for(vertex_t j=0; j<begin[i+1]-begin[i]; j++){
			cout<<head[begin[i]+j]<<"-"<<adj[begin[i]+j]<<" ";
		}
		if(begin[i+1]>begin[i]){
			cout<<"\n";
		}
	}
}

void quickSort(vertex_t * arr, index_t left, index_t right) {
      index_t i = left, j = right;
      vertex_t tmp;
      vertex_t pivot = arr[(left + right) / 2];
 
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };
 
      /* recursion */
      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}


//sort and trim it to upper triangular
void graph::sort(){
	for(vertex_t i=0; i<vert_count; i++){
		index_t a=beg_pos[i];
		index_t b=beg_pos[i+1]-1;
		quickSort(adj_list,a,b);
	}
	
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
	index_t k=0;
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){	//remove multiple edges
			if(j<(beg_pos[i+1]-1) && adj_list[j]==adj_list[j+1] && head_list[j]==head_list[j+1])
			{
				j++;
				continue;
			}
			if(head_list[j]==adj_list[j]){	//remove self loops
				j++;
				continue;
			}
			
			k++;
			upperBegin[i+1]++;//upperDegree[i]++;
			j++;
		}
	}
	
	upperEdgeCount = k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
	k=0;
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){
			if(j<(beg_pos[i+1]-1) &&adj_list[j]==adj_list[j+1]&&head_list[j]==head_list[j+1])
			{
				j++;
				continue;
			}
			if(head_list[j]==adj_list[j]){
				j++;
				continue;
			}
			upperHead[k] =head_list[j];
			upperAdj[k] = adj_list[j];
			k++;
			j++;
		}
	}
	
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}

/*function to search the begin position to find proper place to cut adjacent list
 return value is the smallest vertex id that begin position is equal or larger then x
 */
vertex_t BinarySearch(index_t x, index_t*A, vertex_t bot, vertex_t top){

//	for(int i=bot;i<=top;i++){
//		cout<<A[i]<<" ";
//	}
//	cout<<"\n";

	vertex_t r= (bot+top)/2;
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

void graph::partition(){
	//step 1, evenly cut the upper CSR by using binary search in upperBegin
	step = (vert_count - 1)/4 + 1;
	cout<<"step ="<<step<<endl;
	rem = vert_count-(vert_count-1)/step*step;
	cout<<"rem ="<<rem<<endl;

	partAdj  = new vertex_t*[PART_NUM];
	partHead = new vertex_t*[PART_NUM];
	partBegin  = new index_t*[PART_NUM];
//	partDegree  = new index_t*[PART_NUM];
	partEdgeCount = new index_t[PART_NUM];
	for (index_t k=0;k<16;k++){
		if((k/4)!=3){
//			partBegin[k] =  new index_t[step+1];
			partBegin[k] =  (index_t*) calloc (step+1,sizeof(index_t));
		}
		else{
//			partBegin[k] =  new index_t[rem+1];
			partBegin[k] =  (index_t*) calloc (rem+1,sizeof(index_t));
		}

//		for(long m=0;m<vert_count+1;m++){
//			partBegin[k][m]=0;
//		}
		partEdgeCount[k]=0;
//		partDegree[i] = new index_t[vert_count];
	}
	index_t i,j;
	for(index_t k=0; k<edge_count; k++){
//		cout<<k<<"::"<<head_list[k]<<"-"<<adj_list[k];
		i = head_list[k]/step;
		j = adj_list[k]/step;
		index_t t = i*4+j;

		partEdgeCount[t]++;
		partBegin[t][head_list[k]%step+1]++;
	}
	
	index_t offset[16];
	for(index_t k=0; k<16; k++){
		cout<<"part count of part "<<k<<" is "<<partEdgeCount[k]<<endl;
		partAdj[k] = new vertex_t[partEdgeCount[k]];
		if((k/4)!=3){
			for(index_t m=0;m<step;m++){
				partBegin[k][m+1] += partBegin[k][m];
			}
		}
		else{
			for(index_t m=0;m<rem;m++){
				partBegin[k][m+1] += partBegin[k][m];
			}
		}
		offset[k]=0;
	}
	
	for(index_t k=0; k<edge_count; k++){
		i = head_list[k]/step;
		j = adj_list[k]/step;
		index_t t = i*4+j;

		partAdj[t][offset[t]]=adj_list[k];
		offset[t]++;
	}

	//print all
/*
	for(index_t k=0; k<16;k++){
		for(vertex_t m=0; m<step; m++){
			cout<< m <<" ::  "<<partBegin[k][m]<<"::";
			for(vertex_t n=partBegin[k][m]; n<partBegin[k][m+1]; n++){
				cout<<partAdj[k][n]<<" ";
			}
			cout<<endl;
		}
		cout<< step <<" ::  "<<partBegin[k][step]<<"::"<<endl;
	}
*/
}






/*
	index_t offset[PART_NUM+1];
	offset[0] = 0;
	offset[PART_NUM] = upperEdgeCount;
	for(int i=1; i<PART_NUM; i++){
		index_t k=i*upperEdgeCount/PART_NUM;
		cout<<"k="<<k<<"\n";
		vertex_t index = BinarySearch(k, upperBegin, 0, vert_count-1);
		offset[i] = upperBegin[index];
		cout<<"part "<<i<<" cut at "<<index<<"\n";
	}
	vertex_t *destAdj = new vertex_t[upperEdgeCount];
	vertex_t *destHead= new vertex_t[upperEdgeCount];
//printGraph(vert_count, upperAdj, upperHead, upperBegin);
//	memcpy(destAdj,upperHead,upperEdgeCount*sizeof(int));
//	memcpy(destHead,upperAdj,upperEdgeCount*sizeof(int));
	//we are going to use space on origin int pointer upperAdj and upperHead as output
	//thus we need two new array to store input data.
	
	//step 2, exchange two end points of each edge, get the lower CSC
	for(int i=0;i<PART_NUM;i++){
		
		vertex_t* tempHead = &upperHead[offset[i]];
		vertex_t* tempAdj  = &upperAdj[offset[i]];
		partHead[i] = &destHead[offset[i]];
		partAdj[i]  = &destAdj[offset[i]];

		partEdgeCount[i] = offset[i+1]-offset[i];
cout<<"part "<<i<<" edge "<<partEdgeCount[i]<< "\n";

		partBegin[i] =  new index_t[vert_count+1];
//		partDegree[i] = new index_t[vert_count];
//		memset(partBegin[i],0,(vert_count+1)*sizeof(index_t));
		for(vertex_t j=0; j<vert_count+1; j++){
			partBegin[i][j]=0;
		}
		index_t *partOffset = new index_t[vert_count];
	//step 3, re-organize 1: go through new head list once to get new lowerDegree
		for(index_t j=0; j<partEdgeCount[i]; j++){
			//partDegree
			vertex_t head = tempHead[j];
			//partDegree[i][head]++; 
			partBegin[i][head+1]++;

		}

	//step 4, re-organize 2: prefix, input is degree, output is begin position
		partBegin[i][0]=0;
		for(vertex_t j=0; j<vert_count; j++){
//			partBegin[i][j] = partBegin[i][j-1] + partDegree[i][j-1];
			partBegin[i][j+1] += partBegin[i][j];
			partOffset[j]=0;
		}
//		memset(partOffset, 0, vert_count*sizeof(index_t));

	//step 5, re-organize 3: go through again and moving data to transfer CSC to CSR
		for(index_t j=0; j<partEdgeCount[i]; j++){
			vertex_t head = tempHead[j];
			partHead[i][partBegin[i][head]+partOffset[head]] = tempHead[j];
			partAdj [i][partBegin[i][head]+partOffset[head]] = tempAdj [j];
			partOffset[head]++;
		}
//		printGraph(vert_count, partHead[i], partAdj[i], partBegin[i]);

	}


}
*/
