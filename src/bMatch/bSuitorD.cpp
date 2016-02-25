#include "bMatching.h"
#include <cstring>
//#include <parallel/algorithm>
using namespace std;

 bool comparator(Edge left, Edge right) 
{ 
    return (left.weight > right.weight || (left.weight==right.weight && left.id > right.id)); 
}

int* endT;
int* startT;
int* bT;
Edge* eT;
int* pivotT;

 bool nodeComp1(int i, int j)
{
    return ((pivotT[i]-startT[i]) > (pivotT[j]-startT[j])); 
}

 bool nodeComp3(int i, int j)
{
    return (eT[endT[i]].weight>eT[endT[j]].weight);
}


 bool nodeComp2(int i, int j)
{
    return (eT[startT[i]].weight>eT[startT[j]].weight);
}
void bSuitorBPQD(CSR* G, int *b, int *nlocks, Node* S, int* start, int* end, char* mark, int* order, int type, int stepM, stack_pcl<Param>* St, bool verbose) 
{

    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;
    int* Q1=(int*)malloc(n*sizeof(int));
    int* Q2=(int*)malloc(n*sizeof(int));
    int* Qb1=(int*)malloc(n*sizeof(int));
    int* Qb2=(int*)malloc(n*sizeof(int));
    int* pivot=(int*)malloc(n*sizeof(int));
    int* Qtemp;
    int Qsize=n,Qindx=0,iter=1;
    int numThreads;
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    

	cout<< "Max degree " << G->maxDeg << endl;
    cout<<"StepM: "<<stepM<<endl; 
    cout<<"type: "<<type<<endl;

    Edge** neighbors_per_thread = (Edge**)malloc(numThreads * sizeof(Edge*)); 	
    for(int i=0;i<numThreads;i++)
	neighbors_per_thread[i] = (Edge*)malloc((G->maxDeg + 1 )* sizeof(Edge));

    double t1,t2,t3,t4,t5;
    long long edgeC=0,heapC=0,kickC=0;
    
    #ifdef LBUF
    int** TQ=(int**)malloc(numThreads*sizeof(int*));
    int* tindx=(int*)malloc(numThreads*sizeof(int));

    for(int i=0;i<numThreads;i++)
    {
        tindx[i]=0;
        TQ[i]=(int*)malloc(BSIZE*sizeof(int));
    }
    #endif
    
    /*int* trav=new int[n];
    ifstream inf;
    inf.open("kron.cbin",ios::in|ios::binary);
    if(inf.is_open())
    {
        inf.read((char*)&trav[0],sizeof(int)*n);
        cout<<"Parts loaded..!!"<<endl;
    }*/
    
    if(verbose)
        t1=omp_get_wtime();

    if(type==3)
    {
        #ifdef DYN
            #pragma omp parallel for schedule(guided,CHUNK)
        #else
            #pragma omp parallel for schedule(static, CHUNK)
        #endif
        for(int i=0;i<n;i++)
            order[i]=stepM*b[i]*(i+1);
    }

#ifdef _NEW_OPT_	

	//cout << "With custom_sort_optimized ongoing" << endl;

#if defined (_ACTIVATE_MIC_CODE_) && defined (__MIC__)
	__m512i vzero_int = _mm512_setzero_epi32();
	//__m512i vb_int    = _mm512_set1_epi32(b); //Change for variable b
	__m512i inc = _mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);

	int vertex_16 = (n >> 4) << 4;
	int vertex_remaining = (n - vertex_16);
	__mmask16 end_mask = ((1 << (vertex_remaining)) - 1);

	int local_chunk_size = pow(2, ceil(log(n / (16 * numThreads)) / log(2.0)));
	//cout << "local_chunk_size " << local_chunk_size << endl;

	#pragma omp parallel for num_threads(numThreads) schedule(dynamic, local_chunk_size) //schedule(static, local_chunk_size)
	for(int i = 0; i < vertex_16; i+=16) 
	{
		__m512i vb_int    = _mm512_set1_epi32(b[i]);
        _MM_STORE_INT(nlocks + i, vzero_int);
		_MM_STORE_INT(Qb1 + i, vb_int); // change for variable b (b+i)
		_MM_STORE_INT(Qb2 + i, vzero_int);
		_MM_STORE_INT(Q1 + i, _mm512_add_epi32(inc, _mm512_set1_epi32(i)));	
	}

	/*{
                
                _MM_STORE_INT_MASK(nlocks + vertex_16, vzero_int, end_mask);
                _MM_STORE_INT_MASK(Qb1 + vertex_16, vb_int, end_mask);
                _MM_STORE_INT_MASK(Qb2 + vertex_16, vzero_int, end_mask);
		_MM_STORE_INT_MASK(Q1 + vertex_16, _mm512_add_epi32(inc, _mm512_set1_epi32(vertex_16)), end_mask);
	}*/
#else
    #pragma omp parallel for num_threads(numThreads) schedule(dynamic, CHUNK) 
    for(int i=0;i<n;i++)  
    {
        nlocks[i]=0;            // Initialize locks
        Q1[i]=i;
        Qb1[i]=b[i];
        Qb2[i]=0;
    }
#endif

    //double t1p1 = omp_get_wtime();	
    //if(verbose)
    //    cout << "init Part 1 taking " << t1p1 - t1 << endl;

    #pragma omp parallel num_threads(numThreads) //schedule(dynamic, 32)
    { 
	int threadid = omp_get_thread_num();
    	#pragma omp for //schedule(dynamic)
    	for(int i=0;i<n;i++)  
    	{
		//prefetch((char *)(verInd + ver[i+3]), _MM_HINT_T0);
		//prefetch((char *)(verInd + ver[i+3] + 16), _MM_HINT_T0);

	        S[i].curSize=0;         // current matching=0
        	S[i].maxSize=b[i];         // maximum matching=b
	        S[i].minEntry.id=-1;
        	S[i].minEntry.weight=0.0;

	        if(type!=1)
        	{
	            start[i]=ver[i];    // Adj List start pointer
        	    //end[i]=custom_sort(verInd,ver[i],ver[i+1],trav[i],order,type,&St[tid]);
	            end[i]=custom_sort_optimized(verInd,ver[i],ver[i+1],stepM*b[i],order,type, neighbors_per_thread[threadid]);
        	}
    	}
    }

    //if(verbose)
    //    cout << "init Part 2 taking " << omp_get_wtime() - t1p1 << endl;
#else
    #ifdef DYN
        #pragma omp parallel for schedule(dynamic,CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(int i=0;i<n;i++)  
    {
        
        int tid=omp_get_thread_num();
        nlocks[i]=0;            // Initialize locks
        S[i].curSize=0;         // current matching=0
        S[i].maxSize=b[i];         // maximum matching=b
        S[i].minEntry.id=-1;
        S[i].minEntry.weight=0.0;
        

        if(type!=1)
        {
            start[i]=ver[i];    // Adj List start pointer
            //end[i]=custom_sort(verInd,ver[i],ver[i+1],trav[i],order,type,&St[tid]);
            if(ptype!=1)
                end[i]=custom_sort(verInd,ver[i],ver[i+1],stepM*b[i],order,type,&St[tid]);
            else
                end[i]=pivoting(verInd,ver[i],ver[i+1],pivotID,pivotW, &pivot[i], stepM*b[i]);
        }
        
        Q1[i]=i;
        Qb1[i]=b[i];
        Qb2[i]=0;
        
    }
#endif

    if(type ==1)
    {
        #ifdef DYN
            #pragma omp parallel for schedule(guided,CHUNK)
        #else
            #pragma omp parallel for schedule(static, CHUNK)
        #endif
        for(int i=0;i<m;i++)
            mark[i]=0;
    }
    //cout<<"Initialization Done...!!"<<endl;
    
    pivotT=pivot;
    endT=end;
    startT=start;
    bT=b;
    eT=verInd;
    
    /*if(ptype==1)
        __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp1);
    else 
    {
        if(ptype==2)
           __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp2);
        else
            if(ptype==3)
                __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp3);
    }*/
    
    /*int ttt=0;
    for(int i=0;i<Qsize;i++)
    {
        int index=Q1[i];
        if(end[index]-pivot[index]>5)
        {    
            ttt++;
            cout<<"I see: "<<index<<" "<<pivot[index]<<" "<<end[index]<<" "<<endl;
        }
        if(ttt==20)
            break;
    }*/

    signed int loopc[16];
    signed int lockc[16];
    for(int i=0;i<16;i++)
    {
        loopc[i]=0;
        lockc[i]=0;
    }
    
    if(verbose)
        t2=omp_get_wtime();
    
    cout << "start opt matching " << Qsize << endl;
    
    while(Qsize>0)
    {
   
        t4=omp_get_wtime();
        //#ifdef DYN
        //    #pragma omp parallel for schedule(guided, CHUNK)
        //#else
        //    #pragma omp parallel for schedule(static, CHUNK)
        //#endif

	#pragma omp parallel num_threads(numThreads) 
	{
		int threadid = omp_get_thread_num();

        	#pragma omp for	schedule(guided,CHUNK)		
		for(int i=0;i<Qsize;i++) 
		{
		    float min_w,heaviest,weight;
		    int min_idx,bVer,current,next_vertex,sold,eold,skip;
		    int partnerIdx,partner,y,j,done,tid=omp_get_thread_num();
		    
		    #ifdef LBUF
			int* LQ=TQ[tid];
			int* lindx=&tindx[tid];
		    #endif
		    
		    current=Q1[i];
		    bVer=Qb1[current];
		    Qb1[current]=0;

		    while(bVer>0){ // For Each vertex we want to find bVer=b partners

			partner=-1;
			heaviest=0.0;
			next_vertex=-1;
			done=1;

			sold=start[current];
			eold=end[current];
			j=sold;
			if(type!=1)
			{
			     while(j<end[current]) // Loop over neighbors of the current vertex 
			    {    
				if(verbose)
				    __sync_fetch_and_add(&edgeC,1);;

				y = verInd[j].id;               // y is the neighbor of the current vertex
				weight=verInd[j].weight;        // weight is w(current,y)
				if(weight<=0)
                {
                    heaviest=-1.0;
                    break;
                }
                min_w=S[y].minEntry.weight;
				
				if(weight <=0 || min_w > weight)
				{
				    j++;
				    continue;
				}
				
				while(__sync_lock_test_and_set(&nlocks[y],1)==1)
                    loopc[tid]++;
                lockc[tid]++;
							   
				min_w=S[y].minEntry.weight;
				min_idx=S[y].minEntry.id;

							    
				if((min_w > weight) || (weight == min_w && current < min_idx))
				{    
				    __sync_lock_release(&nlocks[y]);
				    j++;
				}
				else
				{
				    if(verbose)
					    __sync_fetch_and_add(&heapC,1); // heapC++;

				    heaviest=weight;
				    S[y].AddHeap(heaviest,current);    
								
				    __sync_lock_release(&nlocks[y]); //nlocks[y]--;
				    start[current]=j+1;

				    if(min_idx!=-1 && end[min_idx]!=-1)
				    {    
				    if(verbose)
					    __sync_fetch_and_add(&kickC,1);;
					
                    #ifdef LBUF
					    
					if(__sync_fetch_and_add(&Qb2[min_idx],1)==0)
					{   
					    LQ[(*lindx)]=min_idx;
					    (*lindx)++;
					    if((*lindx)==BSIZE)
					    {
						int start=__sync_fetch_and_add(&Qindx,BSIZE);
						int count=0;
						for(int k=start;k<start+BSIZE;k++)
						{    
						    Q2[k]=LQ[count];
						    count++;
						}
						(*lindx)=0;
					    }
					}
					
					#else
					if(__sync_fetch_and_add(&Qb2[min_idx],1)==0)
					   Q2[__sync_fetch_and_add(&Qindx,1)]=min_idx;
					#endif
				    } 
				    
				    break;
				}
			    }   // while(j<end[current])
			    if(heaviest<=0)
			    {
				    
				start[current]=j;
				if(end[current]<ver[current+1] && heaviest==0)
				{
				    if(type==3)
				    {
					if(start[current]+stepM*b[current]<ver[current+1])
					{
					    end[current]=start[current]+stepM*b[current];
					    sort(verInd+start[current],verInd+start[current]+stepM*b[current],comparator);
					}
					else
#ifdef _NEW_OPT_
					    end[current]=custom_sort_optimized(verInd,start[current],ver[current+1],stepM*b[current],order,type+1, neighbors_per_thread[threadid]);
#else
					    end[current]=custom_sort(verInd,start[current],ver[current+1],stepM*b[current],order,type+1,&St[tid]);
#endif

				    }
				    else
#ifdef _NEW_OPT_
					end[current]=custom_sort_optimized(verInd,start[current],ver[current+1],stepM*b[current],order,type, neighbors_per_thread[threadid]);
#else		
					end[current]=custom_sort(verInd,start[current],ver[current+1],stepM*b[current],order,type,&St[tid]);
#endif
			    
				    done=0;
					
                    #ifdef LBUF
					    
					if(__sync_fetch_and_add(&Qb2[current],bVer)==0)
					{   
					    LQ[(*lindx)]=current;
					    (*lindx)++;
					    if((*lindx)==BSIZE)
					    {
						int start=__sync_fetch_and_add(&Qindx,BSIZE);
						int count=0;
						for(int k=start;k<start+BSIZE;k++)
						{    
						    Q2[k]=LQ[count];
						    count++;
						}
						(*lindx)=0;
					    }
					}
					
					#else
					if(__sync_fetch_and_add(&Qb2[current],bVer)==0)
					   Q2[__sync_fetch_and_add(&Qindx,1)]=current;
					#endif
					
				}
				else
				    end[current]=-1;
			    }
			}
			else /// Unsorted Mode
			{
			    j=ver[current];
			    
			    for(;j<ver[current+1] && j!=-1;j++) { // Loop over neighbors of the current vertex 
				
				  if(verbose)
				    __sync_fetch_and_add(&edgeC,1);;

							
				if(mark[j]==1)
				{
				    //if(verbose)
					//__sync_fetch_and_add(&markC,1);;

				    continue;
				}
				
				y = verInd[j].id;               // y is the neighbor of the current vertex
				weight=verInd[j].weight;        // weight is w(current,y)
				min_w=S[y].minEntry.weight;

				if((weight<heaviest)|| (weight == heaviest && y < partner)||min_w>weight)
				    continue;

				if(min_w==weight) 
				{        
				    skip=0;
				
				    while(__sync_lock_test_and_set(&nlocks[y],1)==1);
				
				    min_w=S[y].minEntry.weight;
				    min_idx=S[y].minEntry.id;
				    if((min_w > weight)||(weight == min_w && current < min_idx))
					skip=1;
				    __sync_lock_release(&nlocks[y]);
				   
				    if(skip==1)
					continue;
				}
			       
				heaviest = weight;      // Store the weight of the heaviest edge found so far
				partner = y;
				partnerIdx=j;
				
				
			    } // loop over neighbors
			    
			    if (heaviest > 0) // True if there is a new partner
			    {    
				if(__sync_lock_test_and_set(&nlocks[partner],1)==0) //Locking partner
				{
				    min_w=S[partner].minEntry.weight;
				    min_idx=S[partner].minEntry.id;

				    if(!S[y].find_id(current) && ((heaviest > min_w)||((heaviest==min_w)&&(current>min_idx)))) // Need to re check again because during time of locking someone else can change the state
				    {
					if(verbose)
					    __sync_fetch_and_add(&heapC,1);;

					next_vertex=min_idx;
					S[partner].AddHeap(heaviest,current);
				    } 
				    else // Got the lock but someone already changed the state so search partner all over again
					next_vertex=current;
				     
				    mark[partnerIdx]=1;
				    __sync_lock_release(&nlocks[partner]); // Unlocking partner
	   
				}
				else // Missed the lock, so someone else may or may not changed the state. Hence, RE DO the partner search
				    next_vertex=current;
		      
				if (next_vertex != -1 && next_vertex!=current)  // True if current vertex just kicked another vertex which is alive
				{  
					#ifdef LBUF

					if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
					{   
					    LQ[(*lindx)]=next_vertex;
					    (*lindx)++;
					    if((*lindx)==BSIZE)
					    {
						int start=__sync_fetch_and_add(&Qindx,BSIZE);
						int count=0;
						for(int k=start;k<start+BSIZE;k++)
						{    
						    Q2[k]=LQ[count];
						    count++;
						}
						(*lindx)=0;
					    }
					}
				
					#else
					if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
					   Q2[__sync_fetch_and_add(&Qindx,1)]=next_vertex;
					#endif
				}
			    }
			    else
				end[current]=-1;    //This means the neighbors list is exhausted, this vertex will never be considered again (dead).!!
			}
			
			if(end[current]==-1)  // if the vertex is alive
			    break;     // Do not forget to decrease bVer ...!!
			else
			    if(next_vertex!=current && done==1)
				bVer--;
                if(type!=1 && done ==0)
                    break;
		    } // while(bVer)
		} // loop over vertices
	}       
 
        #ifdef LBUF


        #pragma omp parallel for
        for(int i=0;i<numThreads;i++)
        { 
            int* LQ=TQ[i];
            int* lindx=&tindx[i];
            if((*lindx)>0)
            {    
                int start=__sync_fetch_and_add(&Qindx,(*lindx));
                int count=0;
                for(int k=start;k<start+(*lindx);k++)
                {        
                    Q2[k]=TQ[i][count];
                    count++;
                }
                (*lindx)=0;
            }
        }

        #endif
        Qtemp=Q1;
        Q1=Q2;
        Q2=Qtemp;
        Qtemp=Qb1;
        Qb1=Qb2;
        Qb2=Qtemp;
        Qsize=Qindx;
        Qindx=0;
        iter++;
        
        /*if(ptype==1)
            __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp1);
        else 
        {
            if(ptype==2)
               __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp2);
            else
                if(ptype==3)
                    __gnu_parallel::sort(&Q1[0],&Q1[Qsize],nodeComp3);
        }*/
        //bT=b;
        endT=end;
        startT=start;
        if(verbose)
            cout<<"Iteration: "<<iter-1<<" "<<Qsize<<" "<<omp_get_wtime()-t4<<endl;
    }

    cout<<"Matching Done....!!"<<endl;

    if(verbose)
    {
        t3=omp_get_wtime();

        /*float avTrav=0,maxTrav=0,vcount=0;
        int indx=-1;
        int* trav=new int[n];
        for(int i=0;i<n;i++)
        {
            if(end[i]==-1)
                trav[i]=ver[i+1]-ver[i];
            else
                trav[i]=start[i]-ver[i];
            if(maxTrav<trav[i])
            {    
                maxTrav=trav[i];
                indx=ver[i+1]-ver[i];
            }
            avTrav+=trav[i];
            if(trav[i]>stepM*b)
                vcount++;
        }*/

        /*ofstream of;
        of.open("kron.cbin",ios::out|ios::binary);
        of.write((char*)&trav[0],sizeof(int)*n);*/

        //cout<<"Total Kicks: "<<sum<<endl;
        //cout<<"Average Kick Length :"<<sum*1.0/count<<endl;
        //cout<<"Longest Kick: "<<t<<endl;
        //cout<<"Graph Info: "<<n<<" "<<m<<" "<<mdeg<<" "<<m*1.0/n<<endl;
        //cout<<"Max Adj List traversed: "<<maxTrav<<" ("<<indx<<")"<<endl;
        //cout<<"Avg Adj List traversed: "<<avTrav*1.0/n<<endl;
        //cout<<"Adj List traversed > stepM*b: "<<vcount<<" "<<vcount/n*100<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        cout<<"# of Iteration: "<<iter-1<<endl;
        cout<<"Matching Time: "<<t3-t2<<endl;
        //cout<<"Matching Time: "<<t3-t2<<" ("<<htime<<" "<<htime/(t3-t2)*100<<"% "
            //<<ntime<<" "<<ntime/(t3-t2)*100<<"% {"<<rtime/ntime*100<<"% "<<ltime/ntime*100<<"%}) "<<endl;
        //cout<<"Segment Timing: "<<max/(t2-t1)*100<<" "<<min/(t2-t1)*100<<" "<<avg/(t2-t1)*100<<endl;
        cout<<"Counts: "<<edgeC<<" "<<" "<<heapC<<" "<<kickC<<endl;
        signed int tloop=0,tlock=0;
        for(int i=0;i<16;i++)
        {
            tloop+=loopc[i];
            tlock+=lockc[i];
        }
        cout<<"HLE: "<<tloop<<" "<<tlock<<" "<<tloop*1.00/tlock<<endl;
        //cout<<"Avg counts: "<<edgeC/n<<" "<<heapC/n<<endl;
        cout<<"Total Time: "<<t3-t1<<endl;
    } 

    for(int i=0;i<numThreads;i++)
        free(neighbors_per_thread[i]);
	
    free(neighbors_per_thread);

}// end of bSuitor
