#include "bMatching.h"
#include "mtxReader.h"
#include "omp.h"
using namespace std;

void localDom(CSR* G, Node* S, int* b, int* start, 
        int* end, int* alive, int stepM, int tsort, bool verbose)
{

    double t1=omp_get_wtime(),t2,t3,t4,t5,t6;    
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;

    int* exhaust=(int*)malloc(n*sizeof(int));
    int matched=0, iter=0;
    long long Ecount=0;
    int numThreads;
    #pragma omp parallel
    numThreads=omp_get_num_threads();

    Edge** neighbors_per_thread = (Edge**)malloc(numThreads * sizeof(Edge*)); 	
    for(int i=0;i<numThreads;i++)
	neighbors_per_thread[i] = (Edge*)malloc((G->maxDeg + 1 )* sizeof(Edge));
    
    #pragma omp parallel for
    for(int i=0;i<n;i++)
    { 
		int threadid = omp_get_thread_num();
        
        alive[i]=1;
        exhaust[i]=0;
        S[i].curSize=0;
        S[i].maxSize=b[i]+2;
        S[i].minEntry.id=-1;
        S[i].minEntry.weight=0.0;

        if(tsort!=1)
        {
            start[i]=ver[i];
            #ifdef _NEW_OPT_
		    end[i]=custom_sort_optimized(verInd,ver[i],ver[i+1],stepM*b[i],NULL,tsort, 
                    neighbors_per_thread[threadid]);
            #else
            end[i]=custom_sort(verInd,ver[i],ver[i+1],stepM*b[i],NULL,tsort,NULL);
            #endif
        }
    }
    
    int flag=1;

    t2=omp_get_wtime();

	double t4_sum = 0, t5_sum = 0; 
 
    if(tsort==1)
    {
        cout<<"Unsorted Mode "<<n<<endl;
        while(true)
        {
            t4=omp_get_wtime();
            flag=0;

            ///// Phase 1: setting pointer to heaviest neighbor
            #pragma omp parallel for
            for(int i=0;i<n;i++)
            {
                float heaviest=0.0,weight;
                int partner=-1,id,pidx;
            
                if(b[i]<=0 || alive[i]==0)
                    continue;
                flag=1;

                heaviest=0.0;
                partner=-1;
                for(int j=ver[i];j<ver[i+1];j++)
                {
                    Ecount++;
                    id=verInd[j].id;
                    weight=verInd[j].weight;

                    if(S[i].find_id(id))
                        continue;
                    if(b[id]>0 && (weight>heaviest || (weight==heaviest && id>partner)))
                    {
                        heaviest=weight;
                        partner=id;
                        pidx=j;
                    }
                }

                if(heaviest>0)
                    S[i].AddHeap(heaviest,partner);
                else
                    alive[i]=0;
            }
            
            t5=omp_get_wtime();
            
            //if(verbose)
                //cout<<"Iteration: "<<t5-t4<<" "<<t6-t5<<" "<<t6-t4<<endl;;
 
            ///// Phase 2: Checking whether neighbor also pointing back
            #pragma omp parallel for
            for(unsigned int i=0;i<n;i++)
            {
                //cout<<"Start,";
                
                if(alive[i]==0)
                    continue;

                int cur=S[i].curSize-1;
                int nbor=S[i].heap[cur].id;

                if(S[nbor].find_id(i))
                    b[i]--;
                else
                    S[i].curSize--;
                //cout<<"End "<<i<<endl;
            }

            iter++;
            t6=omp_get_wtime();
            
            //if(verbose)
                //cout<<"Iteration: "<<t5-t4<<" "<<t6-t5<<" "<<t6-t4<<endl;;
            
            if(flag==0)
                break;
        }
    }
    else
    {
        cout<<"sorting type: "<<tsort<<endl;
        while(true)
        {
            t4=omp_get_wtime();
            flag=0;

            ///// Phase 1: setting pointer to heaviest neighbor
	    #ifdef _LOCAL_DOM_OPT_
            #pragma omp parallel for schedule(dynamic, CHUNK)
            #else
            #pragma omp parallel for
            #endif
            //#pragma omp parallel for
            for(int i=0;i<n;i++)
            {
                float heaviest=0.0;
                int partner=-1,j,id,done=1;
		int threadid = omp_get_thread_num();
                
                if(b[i]<=0 || alive[i]==0)
                    continue;
                flag=1;

                while(true)
                {
                    for(j=start[i];j<end[i];j++)
                    {
                        Ecount++;
                        id=verInd[j].id;
                        if(b[id]>0)
                            break;
                    }
                    
                    exhaust[i]=j+1;
                    if(j<end[i])
                    {
                        heaviest=verInd[j].weight;
                        partner=verInd[j].id;
                        break;
                    }
                    else
                    {
                        if(tsort==4)
                        {
                            if(end[i]<ver[i+1])
                            {
                                start[i]=end[i];
                                #ifdef _NEW_OPT_
                                end[i]=custom_sort_optimized(verInd,start[i],ver[i+1],
                                       stepM*b[i],NULL,tsort,
                                       neighbors_per_thread[threadid]);
                                #else
                                end[i]=custom_sort(verInd,start[i],ver[i+1],stepM*b[i],
                                        NULL,tsort,NULL);
                                #endif
                            }
                            else
                                break;
                        }
                        else
                            break;
                    }
                }

                if(heaviest>0)
		{
                    //S[i].AddHeap(heaviest,partner);

		     S[i].heap[S[i].curSize].weight=heaviest;
		     S[i].heap[S[i].curSize].id=partner;
		     S[i].curSize++;
		}
                else
                    alive[i]=0;
            }
           
	    t5=omp_get_wtime();
            ///// Phase 2: Checking whether neighbor also pointing back
            //int min_size = 1024;	

	    #ifdef _LOCAL_DOM_OPT_			
	    #pragma omp parallel for schedule(dynamic, CHUNK)
	    #else	
            #pragma omp parallel for
	    #endif
            for(int i=0;i<n;i++)
            {
                if(alive[i]==0)
                    continue;
	
		#if 0 //_LOCAL_DOM_OPT_ 

		Node* nd = &S[i]; 	
	
		int cur=nd->curSize-1;
                int nbor=nd->heap[cur].id;

		Node* nd_nbor = &S[nbor];
		int j;

		for(j=0; j < nd_nbor->curSize;j++)
                {	
			if(nd_nbor->heap[j].id==i)
			{
				b[i]--;
				start[i]=exhaust[i];
				break;
			}
		}

		if(j == nd_nbor->curSize)	
		{
			nd->curSize--;
			start[i]=exhaust[i]-1;
		}


		#else
		int cur=S[i].curSize-1;
                int nbor=S[i].heap[cur].id;

                if(S[nbor].find_id(i))
                {    
                    b[i]--;
                    start[i]=exhaust[i];
                }
                else
                {    
                    S[i].curSize--;
                    start[i]=exhaust[i]-1;
                }

		#endif
            }

            iter++;
            t6=omp_get_wtime();
            
            if(verbose)
                cout<<"Iteration "<<iter<<": "<<t5-t4<<" "<<t6-t5<<" "<<t6-t4<<endl;;
            //cout << " " << min_size;    

	t4_sum+= t5-t4;
	t5_sum+= t6-t5;            

            if(flag==0)
                break;
        }
    
    
    }
    t3=omp_get_wtime();

    if(verbose)
    {
        cout<<"Edge Count: "<<Ecount<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        cout<<"Matching Time: "<<t3-t2<<endl;
        cout<<"Total Time: "<<t3-t1<<endl;

	cout << "t4_sum " << t4_sum << " t5_sum " << t5_sum << endl;
    }

}
