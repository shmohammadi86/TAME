#include "bMatching.h"
#include "mtxReader.h"
using namespace std;

static int Ecount;
void DP_walk(CSR* G, Node* S, Walker* W, int* mark, int u, int l, int pc)
{

    //cout<<"DP: "<<u<<endl;
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;

    float heaviest=0.0,weight,temp;
    int partner=-1,id,pidx;
    float * tp;
    for(int j=ver[u];j<ver[u+1];j++)
    {
        Ecount++;
        
        weight=verInd[j].weight;
        id=verInd[j].id;

        if(mark[id]==pc || mark[id]==-1)
            continue;
        if(S[id].find_id(u) || S[u].find_id(id))
            continue;
        
        if(S[id].curSize == S[id].maxSize)
            continue;

        if(weight> heaviest || (weight==heaviest && id>partner))
        {
            heaviest=weight;
            partner=id;
            pidx=j;
        }
    }
    
    if(heaviest>0)
    {
        mark[partner]=pc;
        if(l==1)
        {
            W->W1=heaviest;
            W->M1[W->midx1]=u;
            W->M1[W->midx1+1]=partner;
            W->M1[W->midx1+2]=heaviest;
            W->midx1+=3;
        }
        else
        {
            if(W->W2 + heaviest > W->W1)
            {
                temp=W->W2;
                W->W2=W->W1;
                W->W1=temp+heaviest;

                W->M2[W->midx2]=u;
                W->M2[W->midx2+1]=partner;
                W->M2[W->midx2+2]=heaviest;
                W->midx2+=3;

                tp=W->M2;
                W->M2=W->M1;
                W->M1=tp;

                pidx=W->midx2;
                W->midx2=W->midx1;
                W->midx1=pidx;

            }
            else
            {
                W->W2=W->W1;
                for(int i=0;i<W->midx1;i++)
                    W->M2[i]=W->M1[i];
                W->midx2=W->midx1;
            }
        }

        DP_walk(G,S,W,mark,partner,l+1,pc);
    }
}


double PGDP(CSR* G, Node* S, int* b, Walker* W, int* mark, bool verbose)
{
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;
    float MatchedWeight=0.0;
    double t1,t2,t3,dpt=0.0,ht=0.0,tt,ft=0.0,tt1;
    int pathCount=0;
    int * Q=(int*)malloc(n*sizeof(int));
    int numThreads;
    #pragma omp parallel
    numThreads=omp_get_num_threads();

    if(verbose)
        t1=omp_get_wtime();

    for(int i=0;i<n;i++)
    {
        S[i].curSize=0;         // current matching=0
        S[i].maxSize=*b;         // maximum matching=b
        S[i].minEntry.id=-1;
        S[i].minEntry.weight=0.0;
        mark[i]=0;

        Q[i]=i;
    }

    //srand(time(NULL));
    //random_shuffle(&Q[0],&Q[n]);

    int flag=1;
    float u,v,w;
    //cout<<"Looping"<<endl;
    Ecount=0;
    if(verbose)
        t2=omp_get_wtime();

    int iter=0;
    while(true)
    {
        flag=0;

        iter++;

        if(verbose)
            tt1=omp_get_wtime();
        
        for(int k=0;k<n;k++)
        {
            int i=Q[k];
            if(mark[i]==-1)
                continue;

            if(S[i].curSize<S[i].maxSize)
            {
                flag=1;
                
                pathCount++;

                W->W1=0.0;
                W->W2=0.0;
                W->midx1=0;
                W->midx2=0;
                
                //cout<<"Calling with: "<<i<<" "<<S[i].curSize<<endl;
                if(numThreads==1)
                    mark[i]=pathCount;
                else
                    mark[i]=iter;

                
                if(verbose)
                    tt=omp_get_wtime();

                if(numThreads==1)
                    DP_walk(G,S,W,mark,i,1,pathCount);
                else
                    DP_walk(G,S,W,mark,i,1,iter);
                
                if(verbose)
                    dpt+=(omp_get_wtime()-tt);
                
                if(verbose)
                    tt=omp_get_wtime();

                if(W->midx1>0)
                {
                    MatchedWeight+=W->W1;
                    for(int j=0;j<W->midx1;j+=3)
                    {
                      u=W->M1[j];
                      v=W->M1[j+1];
                      w=W->M1[j+2];  
                      S[(int)u].AddHeap(w,(int)v);
                      S[(int)v].AddHeap(w,(int)u);
                    }
                }
                else
                    mark[i]=-1;

                if(verbose)
                    ht+=(omp_get_wtime()-tt);
            }
            else
                mark[i]=-1;
        }
        if(verbose)
        {
            tt=omp_get_wtime();
            cout<<"Iteration "<<iter<<":"<<tt-tt1<<endl;
        }

        if(flag==0)
            break;
    }
    if(verbose)
    {
        t3=omp_get_wtime();
        //cout<<"DP Time: "<<dpt<<endl;
        //cout<<"Heap Time: "<<ht<<endl;
        cout<<"Number of Paths: "<<pathCount<<endl;
        cout<<"Counts: "<<Ecount<<endl;
        cout<<"Initialization Time: "<<t2-t1<<endl;
        cout<<"Matching Time: "<<t3-t2<<endl;
        cout<<"Total Time: "<<t3-t1<<endl;
    }
    return MatchedWeight;
}
