#include "bMatching.h"
#include "mtxReader.h"
using namespace std;

static long long Ecount;
void find_walk(CSR* G, Node* S, Node* S1, int* b, float* M, char* mark, int* curDeg, int u, int l)
{
    if(b[u]>0)
        b[u]--;
    else
        return;
    
    if(curDeg[u]<=0)
        return;

    /*if(l%2==0)
        cout<<"Find walk: "<<b[u]<<" "<<M[0]<<" "<<curDeg[u]<<" "<<u<<" "<<l<<endl;
    else
        cout<<"Find walk: "<<b[u]<<" "<<M[1]<<" "<<curDeg[u]<<" "<<u<<" "<<l<<endl;
    */
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;

    float heaviest=0.0,weight;
    int partner=-1,id,pidx;
    for(int j=ver[u];j<ver[u+1];j++)
    {
        Ecount++;
        if(mark[j]==1)
            continue;
        
        weight=verInd[j].weight;
        id=verInd[j].id;

        if(S[id].find_id(u))
            continue;
        
        if(b[id]>0 && (weight> heaviest || (weight==heaviest && id>partner)))
        {
            heaviest=weight;
            partner=id;
            pidx=j;
        }
    }
    
    if(heaviest>0)
    {
        mark[pidx]=1;
        /*for(int j=ver[partner];j<ver[partner+1];j++)
        {    
            Ecount++;
            if(verInd[j].id==u)
            {
                mark[j]=1;
                break;
            }
        }*/
        curDeg[partner]--;
        curDeg[u]--;
        if(l%2==1)
        {
            M[1]+=heaviest;
            S1[u].AddHeap(heaviest,partner);
            S1[partner].AddHeap(heaviest,u);
        }
        else
        {
            M[0]+=heaviest;
            S[u].AddHeap(heaviest,partner);
            S[partner].AddHeap(heaviest,u);

        }
        if(b[u]==0)
        {
            curDeg[u]=0;
            /*for(int j=ver[u];j<ver[u+1];j++)
            {
                Ecount++;
                mark[j]=1;
                curDeg[u]=0;
            }*/
        }
        find_walk(G,S,S1,b,M,mark,curDeg,partner,l+1);
    }
}


int PG(CSR* G, Node* S, Node* S1, int* b, float* M, char* mark, int* curDeg)
{
    int n=G->nVer;
    int m=G->nEdge;
    int* ver=G->verPtr;
    Edge* verInd=G->verInd;
    double t1,t2,t3;

    t1=omp_get_wtime();

    M[0]=0.0;
    M[1]=0.0;

    for(int i=0;i<n;i++)
    {
        S[i].curSize=0;         // current matching=0
        S[i].maxSize=b[i]+1;         // maximum matching=b
        S[i].minEntry.id=-1;
        S[i].minEntry.weight=0.0;
        
        S1[i].curSize=0;         // current matching=0
        S1[i].maxSize=b[i]+1;         // maximum matching=b
        S1[i].minEntry.id=-1;
        S1[i].minEntry.weight=0.0;

        curDeg[i]=ver[i+1]-ver[i];
        for(int j=ver[i];j<ver[i+1];j++)
            mark[j]=0;
    }

    int flag=1;

    Ecount=0;

    t2=omp_get_wtime();

    while(true)
    {
        flag=0;

        for(int i=0;i<n;i++)
        {
            if(b[i]>0 && curDeg[i]>0)
            {
                flag=1;
                find_walk(G,S,S1,b,M,mark,curDeg,i,0);
            }
        }

        if(flag==0)
            break;
    }

    t3=omp_get_wtime();

    cout<<"Edge Count: "<<Ecount<<endl;
    cout<<"Initialization Time: "<<t2-t1<<endl;
    cout<<"Matching Time: "<<t3-t2<<endl;
    cout<<"Total Time: "<<t3-t1<<endl;

    if(M[0]>M[1])
        return 0;
    else
        return 1;
}
