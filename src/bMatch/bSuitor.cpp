#include "bMatching.h"
#include <cstring>
using namespace std;

void bSuitorBPQ(CSR* G, int *b, int *nlocks, Node* S, int* start, int* end, char* mark, int* order, int type, int stepM,stack_pcl<Param>*St,bool verbose) 
{

    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;
    int numThreads;
    
    double t1,t2,t3;
    long long edgeC=0,heapC=0,kickC=0;
    
    if(verbose)
    {
        #pragma omp parallel
        numThreads=omp_get_num_threads();
        t1=omp_get_wtime();
    }
    
    if(type ==1 || type == 3 || type==4)
    {
       #ifdef DYN
            #pragma omp parallel for schedule(guided,CHUNK)
        #else
            #pragma omp parallel for schedule(static, CHUNK)
        #endif
        for(int i=0;i<m;i++)  
            mark[i]=0;
    }

    #ifdef DYN
        #pragma omp parallel for schedule(guided,CHUNK)
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
            end[i]=custom_sort(verInd,ver[i],ver[i+1],stepM*b[i],order,type,&St[tid]); 
            
        }
    }
    
    if(verbose)
        t2=omp_get_wtime();
    
    #ifdef DYN
        #pragma omp parallel for schedule(guided, CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(int i=0;i<n;i++) 
    {
        float min_w,heaviest,weight;
        int min_idx,bVer=b[i],done,current,sold,eold;
        int partnerIdx,partner,y,next_vertex,j,skip;
        int tid=omp_get_thread_num();
        
        while(bVer>0){ // For Each vertex we want to find bVer=b partners
            done = 0;
            current = i;

            while (done==0) {
                
                partner=-1;
                heaviest=0.0;
                done=1;
                next_vertex=-1;

                                
                sold=start[current];
                eold=end[current];
                j=sold;
                
                if(sold!=-1 && type!=1)
                {
                    while(j<end[current]) // Loop over neighbors of the current vertex 
                    {    
                        if(verbose)
                            __sync_fetch_and_add(&edgeC,1);

                        y = verInd[j].id;               // y is the neighbor of the current vertex
                        weight=verInd[j].weight;        // weight is w(current,y)
                        min_w=S[y].minEntry.weight;

                        if(weight <=0 || min_w > weight)
                        {    
                            j++;
                            continue;
                        }
 
                        while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                        
                        min_w=S[y].minEntry.weight;
                        min_idx=S[y].minEntry.id;
    
                        if( S[y].find_id(current)||(min_w > weight)||(weight == min_w && current < min_idx))
                        {    
                           __sync_lock_release(&nlocks[y]);
                            j++;
                        }
                        else
                        {
                            if(verbose)
                                __sync_fetch_and_add(&heapC,1);
  
                            heaviest=weight;
                            S[y].AddHeap(heaviest,current);    
                                                        
                            __sync_lock_release(&nlocks[y]);
                           

                           start[current]=j+1;
                            
                            if(min_idx!=-1)
                            {    
                                current=min_idx;
                                done=0;
                            if(verbose)
                                __sync_fetch_and_add(&kickC,1);
                            } 
                            
                            break;
                            
                        }
                        
                    }   // while(j<end[current])
                    
                    if(heaviest==0)
                    {
                        if(type==2)
                            start[current]=j;
                        else
                        {
                            done=0;
                            start[current]=-1;
                        }
                    }
                }
                else
                {
                    if(type==1) 
                        j=ver[current];
                    else
                        j=end[current];
                    
                    for(;j<ver[current+1] && j!=-1;j++) { // Loop over neighbors of the current vertex 
                        
                        if(verbose)
                            __sync_fetch_and_add(&edgeC,1);

                                                
                        if(mark[j]==1)
                        {    
                            //if(verbose)
                                //__sync_fetch_and_add(&markC,1);
                            continue;
                        }
                        
                        y = verInd[j].id;               // y is the neighbor of the current vertex
                        weight=verInd[j].weight;        // weight is w(current,y)
                        min_w=S[y].minEntry.weight;
                        
                        if((weight<heaviest)||(weight == heaviest && y < partner) || min_w > weight)
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
                        while(__sync_lock_test_and_set(&nlocks[partner],1)==1); //Locking partner
                        
                        min_w=S[partner].minEntry.weight;
                        min_idx=S[partner].minEntry.id;

                        if(mark[partnerIdx]==0 &&((heaviest > min_w)||((heaviest==min_w)&&(current >  min_idx)))) // Need to re check again because during time of locking someone else can change the state
                        {
                            if(verbose)    
                                __sync_fetch_and_add(&heapC,1);

                            next_vertex=min_idx;
                            S[partner].AddHeap(heaviest,current);
                        } 
                        else // Got the lock but someone already changed the state so search partner all over again
                            next_vertex=current;
                        
                        mark[partnerIdx]=1;
                        __sync_lock_release(&nlocks[partner]); // Unlocking partner
                            
                        if (next_vertex != -1)  // True if current vertex just kicked another vertex which is alive
                        {       
                            current = next_vertex;      // Pick up the old suitor and thread will continue with it (EAGER UPDATE)
                            done = 0;
                            if(verbose)
                                __sync_fetch_and_add(&kickC,1);
                        }
                    }
                    else
                        end[current]=-1;    //This means the neighbors list is exhausted, this vertex will never be considered again (dead).!!
                }
            } // while(!done) loop
            
            if(end[i]==-1)  // if the vertex is dead
                break;     
            else
                bVer--;     // if Alive, decrease bVer.
        } // while(bVer)

    } // loop over vertices

    cout<<"Matching Done....!!"<<endl;
    
    if(verbose)
    {
        t3=omp_get_wtime();

        //cout<<"Total Kicks: "<<sum<<endl;
        //cout<<"Average Kick Length :"<<sum*1.0/count<<endl;
        //cout<<"Longest Kick: "<<t<<endl;
        //cout<<"Graph Info: "<<n<<" "<<m<<" "<<mdeg<<" "<<m*1.0/n<<endl;
        //cout<<"Max Adj List traversed: "<<vkick-1<<" ("<<deg<<")"<<endl;
        //cout<<"Avg Adj List traversed: "<<avtrav*1.0/n<<endl;
        //cout<<"Adj List traversed > stepM*b: "<<avkick*1.0/n*100<<"%"<<" "<<avtemp*1.0/n*100<<"%"<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        cout<<"Matching Time: "<<t3-t2<<endl;
        //cout<<"Matching Time: "<<t3-t2<<" ("<<htime<<" "<<htime/(t3-t2)*100<<"% "
            //<<ntime<<" "<<ntime/(t3-t2)*100<<"% {"<<rtime/ntime*100<<"% "<<ltime/ntime*100<<"%}) "<<endl;
        //cout<<"Segment Timing: "<<max/(t2-t1)*100<<" "<<min/(t2-t1)*100<<" "<<avg/(t2-t1)*100<<endl;
        cout<<"Edge Count and Heap Count: "<<edgeC<<" "<<heapC<<" "<<kickC<<endl;
        //cout<<"Avg counts: "<<edgeC/n<<" "<<heapC/n<<endl;
        cout<<"Total Time: "<<t3-t1<<endl;
    } 
}// end of bSuitor
