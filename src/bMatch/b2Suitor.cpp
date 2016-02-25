#include "bMatching.h"
#include <cstring>
using namespace std;


//#define LOAD
//#define STATISTICS
#define CHUNK 64
#define DYN

void b2Suitor(CSR* G, int b, int *nlocks, Node* S, int* start, int* end, char* mark, int type) 
{

    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;
    int numThreads;
    #ifdef STATISTICS 
        double t1,t2,t3;
        #pragma omp parallel
        numThreads=omp_get_num_threads();

        double* heapTime = new double[16*numThreads];
        double* nborTime = new double[16*numThreads];
        double* redoTime = new double[16*numThreads];
        double* nborLock = new double[16*numThreads];
        int* sumKick=new int[16*numThreads];
        int* countKick=new int[16*numThreads];
        int* longestKick=new int[16*numThreads];
        int* verKick= new int[n];
        
        #pragma omp parallel for
        for(int i=0;i<numThreads;i++)
        {
            longestKick[16*i]=0;
            sumKick[16*i]=0;
            countKick[16*i]=0;
            heapTime[16*i]=0.0;
            nborTime[16*i]=0.0;
            redoTime[16*i]=0.0;
            nborLock[16*i]=0.0;
        }
        #pragma omp parallel for
        for(int i=0;i<n;i++)
            verKick[i]=0;
        
        t1=omp_get_wtime();
    #endif

    
    #ifdef LOAD
        #pragma omp parallel
        numThreads=omp_get_num_threads();

        double* mlbTime=new double[16*numThreads];
        double* ilbTime=new double[16*numThreads];
        #pragma omp parallel for
        for(int i=0;i<numThreads;i++)
        {
            mlbTime[i*16]=0.0;
            ilbTime[i*16]=0.0;
        }
    #endif

    #ifdef DYN
        #pragma omp parallel for schedule(guided,CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(int i=0;i<n;i++)  
    {
        #ifdef LOAD
            double lbtemp=omp_get_wtime();
        #endif
        nlocks[i]=0;            // Initialize locks
        S[i].curSize=0;         // current matching=0
        S[i].maxSize=b;         // maximum matching=b
        if(type!=1)
        {
            start[i]=ver[i];    // Adj List start pointer
            end[i]=custom_sort(verInd,ver[i],ver[i+1],5*b,type); 
            
        }
        for(int j=0;j<b;j++)
        {
            S[i].heap[j].weight=-1.0;
            S[i].heap[j].id=-1;
        }

        #ifdef LOAD
            ilbTime[16*omp_get_thread_num()]+=(omp_get_wtime()-lbtemp);
        #endif

    }
    if(type ==1 || type == 3 || type==4)
    {
        //memset(mark,0,m*sizeof(char));
       #ifdef DYN
            #pragma omp parallel for schedule(guided,CHUNK)
        #else
            #pragma omp parallel for schedule(static, CHUNK)
        #endif
        for(int i=0;i<m;i++)  
        {    
            #ifdef LOAD
                double lbtemp=omp_get_wtime();
            #endif
            mark[i]=0;
            #ifdef LOAD
                ilbTime[16*omp_get_thread_num()]+=(omp_get_wtime()-lbtemp);
            #endif
        }
    }

#ifdef STATISTICS
        t2=omp_get_wtime();
    #endif
    
    #ifdef DYN
        #pragma omp parallel for schedule(guided, CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(int i=0;i<n;i++) 
    {
        float min_w,heaviest,weight;
        int min_idx,bVer=b,done,current,sold,eold;
        int partnerIdx,partner,y,next_vertex,j,skip;
        
        #ifdef LOAD
            double lbtemp=omp_get_wtime();
        #endif

        /*#ifdef STATISTICS
            int localKick;
            bool redo=false;
        #endif*/
        
        while(bVer>0){ // For Each vertex we want to find bVer=b partners
            done = 0;
            current = i;

            /*#ifdef STATISTICS
                localKick=0;
            #endif*/

            while (done==0) {
                
                partner=-1;
                heaviest=0.0;
                //cout<<"Done: "<<done<<endl;
                done=1;
                next_vertex=-1;

                #ifdef STATISTICS
                    double temp1=omp_get_wtime();
                #endif
                
                sold=start[current];
                eold=end[current];
                j=sold;
                if(sold!=-1 && type!=1)
                {
                    while(j<end[current]) // Loop over neighbors of the current vertex 
                    {    
                       #ifdef STATISTICS
                            double ltemp=omp_get_wtime();
                        #endif

                        y = verInd[j].id;               // y is the neighbor of the current vertex
                        weight=verInd[j].weight;        // weight is w(current,y)
                        min_w=S[y].heap[0].weight;

                        if(weight <=0)
                        {    
                            j++;
                            continue;
                        }
                        // Check if w(current,y) is the best so far, and if it is a better than the worst partner
                        // Also check that current is not already a suitor of y
                        
                        //while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                        
                        while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                        
                        //min_w=S[y].min_weight();
                        //min_idx=S[y].min_id();
                        min_w=S[y].heap[0].weight;
                        min_idx=S[y].heap[0].id;
                        
                        #ifdef STATISTICS
                            ltemp=omp_get_wtime()-ltemp;
                            nborLock[16*omp_get_thread_num()]+=ltemp;
                        #endif
                            
                        if( S[y].find_id(current) || (min_w > weight)||(weight == min_w && current < min_idx))
                        {    
                            __sync_lock_release(&nlocks[y]);
                            j++;
                        }
                        else
                        {
                            #ifdef STATISTICS
                                temp1=omp_get_wtime()-temp1;
                                nborTime[16*omp_get_thread_num()]+=temp1;
                                /*if(redo)
                                    redoTime[16*omp_get_thread_num()]+=temp1;
                                redo=false;*/
                            #endif

                            
                            #ifdef STATISTICS
                            /*if(next_vertex!=-1)
                            {    
                                localKick++;
                                __sync_fetch_and_add(&verKick[next_vertex],1);
                            }*/
                            double temp=omp_get_wtime();
                            #endif
                            
                            heaviest=weight;
                            S[y].Add(heaviest,current);    
                                                        
                            __sync_lock_release(&nlocks[y]);
                            //__sync_bool_compare_and_swap(&start[current],sold,j+1);
                           start[current]=j+1;
                            
                            if(min_idx!=-1)
                            {    
                                current=min_idx;
                                done=0;
                            } 
                            
                            #ifdef STATISTICS
                                heapTime[16*omp_get_thread_num()]+=(omp_get_wtime()-temp);
                            #endif
  
                            break;
                            
                        }
                        
                    }   // while(j<end[current])
                    if(heaviest==0)
                    {
                        if(type==2)
                            start[current]=j;
                            //__sync_bool_compare_and_swap(&start[current],sold,j);
                        else
                        {
                            done=0;
                            //__sync_bool_compare_and_swap(&start[current],sold,-1);
                            start[current]=-1;
                        }
                        #ifdef STATISTICS
                            temp1=omp_get_wtime()-temp1;
                            nborTime[16*omp_get_thread_num()]+=temp1;
                        #endif

                    }
                }
                else
                {
                    if(type==1) 
                        j=ver[current];
                    else
                        j=end[current];
                    
                    //cout<<"Start: "<<j<<" "<<ver[current]<<endl;
                    for(;j<ver[current+1] && j!=-1;j++) { // Loop over neighbors of the current vertex 
                        
                        
                        if(mark[j]==1)
                            continue;
                        
                        y = verInd[j].id;               // y is the neighbor of the current vertex
                        weight=verInd[j].weight;        // weight is w(current,y)
                        //min_w=S[y].min_weight();
                        min_w=S[y].heap[0].weight;

                        if((weight<heaviest)||(weight == heaviest) && (y < partner))
                            continue;

                        // Check if w(current,y) is the best so far, and if it is a better than the worst partner
                        // Also check that current is not already a suitor of y
              
                        //while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                        if(min_w==weight) 
                        {        
                            skip=0;
                            while( __sync_lock_test_and_set(&nlocks[y],1)==1);
                        
                            min_w=S[y].heap[0].weight;
                            min_idx=S[y].heap[0].id;
                            //min_w=S[y].min_weight();
                            //min_idx=S[y].min_id();
                            
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
                    //cout<<"End"<<endl;
                    #ifdef STATISTICS
                        temp1=omp_get_wtime()-temp1;
                        nborTime[16*omp_get_thread_num()]+=temp1;
                        /*if(redo)
                            redoTime[16*omp_get_thread_num()]+=temp1;
                        redo=false;*/
                    #endif


                    if (heaviest > 0) // True if there is a new partner
                    {            
                        while(__sync_lock_test_and_set(&nlocks[partner],1)==1); //Locking partner
                        
                        min_w=S[partner].heap[0].weight;
                        min_idx=S[partner].heap[0].id;

                        if(mark[partnerIdx]==0 &&((heaviest > min_w)||((heaviest==min_w)&&(current >  min_idx)))) // Need to re check again because during time of locking someone else can change the state
                        {
                            next_vertex=min_idx;
                            
                            #ifdef STATISTICS
                                /*if(next_vertex!=-1)
                                {    
                                    localKick++;
                                    __sync_fetch_and_add(&verKick[next_vertex],1);
                                }*/
                                double temp=omp_get_wtime();
                            #endif

                            if(heaviest>S[partner].heap[1].weight || (heaviest==S[partner].heap[1].weight &&current>S[partner].heap[1].id))
                            {
                                S[partner].heap[0].weight=S[partner].heap[1].weight;
                                S[partner].heap[0].id=S[partner].heap[1].id;
                                S[partner].heap[1].weight=heaviest;
                                S[partner].heap[1].id=current;
                            }
                            else
                            {
                                S[partner].heap[0].weight=heaviest;
                                S[partner].heap[0].id=current;
                            }

                            //S[partner].Add(heaviest,current);

                            #ifdef STATISTICS
                                heapTime[16*omp_get_thread_num()]+=(omp_get_wtime()-temp);
                            #endif
                        } 
                        else // Got the lock but someone already changed the state so search partner all over again
                            next_vertex=current;
                        
                        mark[partnerIdx]=1;
                        __sync_lock_release(&nlocks[partner]); // Unlocking partner
                            
                        
              
                        if (next_vertex != -1)  // True if current vertex just kicked another vertex which is alive
                        {       
                            current = next_vertex;      // Pick up the old suitor and thread will continue with it (EAGER UPDATE)
                            done = 0;
                        }
                    }
                    else
                        end[current]=-1;    //This means the neighbors list is exhausted, this vertex will never be considered again (dead).!!
                }
            } // while(!done) loop
            
            /*#ifdef STATISTICS
                if(localKick>0)
                {
                    int tid=omp_get_thread_num();
                    sumKick[16*tid]+=localKick;
                    countKick[16*tid]++;
                    if(localKick>longestKick[16*tid])
                            longestKick[16*tid]=localKick;
                }
            #endif*/

            if(end[i]==-1)  // if the vertex is alive
                break;     // Do not forget to decrease bVer ...!!
            else
                bVer--;
        } // while(bVer)
        #ifdef LOAD
            mlbTime[16*omp_get_thread_num()]+=(omp_get_wtime()-lbtemp);
        #endif

    } // loop over vertices

    cout<<"Matching Done....!!"<<endl;
    #ifdef LOAD
        double ilbavg=0.0,mlbavg=0.0,ilbstd=0.0,mlbstd=0.0;;
        for(int i=0;i<numThreads;i++)
        {
            ilbavg+=ilbTime[16*i];
            mlbavg+=mlbTime[16*i];
        }
        ilbavg=ilbavg/numThreads;
        mlbavg=mlbavg/numThreads;

        for(int i=0;i<numThreads;i++)
        {
            ilbstd+=pow((ilbavg-ilbTime[16*i]),2);
            mlbstd+=pow((mlbavg-mlbTime[16*i]),2);
        }
        ilbstd=sqrt(ilbstd/numThreads);
        mlbstd=sqrt(mlbstd/numThreads);

        cout<<"Init Threads Times: "<<ilbavg<<" "<<ilbstd<<" "<<ilbstd/ilbavg*100<<"%"<<endl;
        cout<<"matching Threads Times: "<<mlbavg<<" "<<mlbstd<<" "<<mlbstd/mlbavg*100<<"%"<<endl;
    #endif
    #ifdef STATISTICS
        t3=omp_get_wtime();
        int t=0,sum=0,count=0;
        double htime=0.0,ntime=0.0,rtime=0.0,ltime=0.0;
        
        for(int i=0; i<numThreads;i++)
        {    
            if(longestKick[16*i]>t)
                t=longestKick[16*i];
            
            /*if(heapTime[16*i]>htime)
                htime=heapTime[16*i];
            if(nborTime[16*i]>ntime)
                ntime=nborTime[16*i];
            if(redoTime[16*i]>rtime)
                rtime=redoTime[16*i];
            if(nborLock[16*i]>ltime)
                ltime=nborLock[16*i];*/
            htime+=heapTime[16*i]/numThreads;
            ntime+=nborTime[16*i]/numThreads;
            ltime+=nborLock[16*i]/numThreads;

            sum=sum+sumKick[16*i];
            count=count+countKick[16*i];
        }
        
        int vkick=0,avkick=0,vtemp,deg=0,mdeg=0,dtemp,avtemp=0, avtrav=0;;
        for(int i=0;i<n;i++)
        {    
            vtemp=start[i]-ver[i];
            dtemp=ver[i+1]-ver[i];
            avtrav+=vtemp;
            if(dtemp>mdeg)
                mdeg=dtemp;

            if(vtemp>vkick)
            {
                vkick=start[i]-ver[i];
                deg=ver[i+1]-ver[i];               
            }
            if(vtemp>5*b)
                avkick++;
            if(dtemp>5*b)
                avtemp++;
            /*if(verKick[i]>vkick)
                vkick=verKick[i];
            if(verKick[i]>3*b)
                avkick++;*/
        }
        cout<<"Total Kicks: "<<sum<<endl;
        cout<<"Average Kick Length :"<<sum*1.0/count<<endl;
        cout<<"Longest Kick: "<<t<<endl;
        cout<<"Graph Info: "<<n<<" "<<m<<" "<<mdeg<<" "<<m*1.0/n<<endl;
        cout<<"Max Adj List traversed: "<<vkick-1<<" ("<<deg<<")"<<endl;
        cout<<"Avg Adj List traversed: "<<avtrav*1.0/n<<endl;
        cout<<"Adj List traversed > 5*b: "<<avkick*1.0/n*100<<"%"<<" "<<avtemp*1.0/n*100<<"%"<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        //cout<<"Matching Time: "<<t3-t2<<" ("<<htime<<" "<<htime/(t3-t2)*100<<"%)"<<endl;
        cout<<"Matching Time: "<<t3-t2<<" ("<<htime<<" "<<htime/(t3-t2)*100<<"% "
            <<ntime<<" "<<ntime/(t3-t2)*100<<"% {"<<rtime/ntime*100<<"% "<<ltime/ntime*100<<"%}) "<<endl;
        
        cout<<"Total Time: "<<t3-t1<<endl;
    #endif
    
}// end of bSuitor
