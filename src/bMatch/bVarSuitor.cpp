// parallel weighted algorithm
//
// Computes greedy b-matching and performs dynamic programming
// Dynamic programming is done in a separate routine 
// One routine for the paths, and one routine for the cycles


#define DYN 1
#define CHUNK 1


typedef struct _Node 
{
    int maxSize;
    int curSize;
    int* id;
    double* weight;
    
} Node;

void print(Node* n)
{
    int i;
    printf("Printing the heap array\n");
    for(i=0;i<n->curSize;i++)
    {
        printf("%d %f\n",n->id[i],n->weight[i]);
    }
}

int min_id(Node* n)
{
    if(n->curSize>0 && n->curSize==n->maxSize)
        return n->id[0];
    else return 0;
}

double min_weight(Node* n)
{
    if(n->curSize>0 && n->curSize==n->maxSize)
        return n->weight[0];
    else return 0.0;
}

int find_id(Node* n, int id)
{
    int i;
    for(i=0;i<n->curSize;i++)
        if(n->id[i]==id)
            return 1;
    
    return 0;
}

void Add(Node* n, double wt, int id )
{
    if(n->curSize==n->maxSize)
    {
        n->weight[0]=wt;
        n->id[0]=id;

        int small,ri,li,pi=0;
        int done=0;
        while(!done)
        {
            li=2*pi+1;
            ri=2*pi+2;
            small=pi;

            if(li <n->maxSize && (n->weight[li]< n->weight[small] || (n->weight[li] == n->weight[small] && n->id[li] < n->id[small] )))
                small=li;
            if(ri <n->maxSize && (n->weight[ri]< n->weight[small] || (n->weight[ri] == n->weight[small] && n->id[ri] < n->id[small] )))
                small=ri;

            if(small != pi)
            {
                wt=n->weight[pi];
                id=n->id[pi];

                n->weight[pi]=n->weight[small];
                n->id[pi]=n->id[small];

                n->weight[small]=wt;
                n->id[small]=id;
            }
            else done=1;

            pi=small;
        }
    }
    else
    {
        if(n->curSize==0)
        {
            n->weight[0]=wt;
            n->id[0]=id;
        }
        if(n->curSize==1)
        {
            if(n->weight[0]< wt || (n->weight[0] == wt && n->id[0] < id))
            {    
                n->weight[1]=wt;
                n->id[1]=id;
            }
            else
            {
                n->weight[1]=n->weight[0];
                n->id[1]=n->id[0];
                
                n->weight[0]=wt;
                n->id[0]=id;
            }
        }
        
        if(n->curSize==2)
        {
            if(n->weight[0]< wt || (n->weight[0] == wt && n->id[0] < id))
            {    
                n->weight[2]=wt;
                n->id[2]=id;
            }
            else
            {
                n->weight[2]=n->weight[0];
                n->id[2]=n->id[0];
                
                n->weight[0]=wt;
                n->id[0]=id;
            }
        }
        if(n->curSize>2)
        {
            n->weight[n->curSize]=wt;
            n->id[n->curSize]=id;
            
            int pi,ci=n->curSize;
            int done=0;
            
            while(!done)
            {
                pi=(int)floor((ci-1)/2.0);

                if(n->weight[pi]< n->weight[ci] || (n->weight[pi] == n->weight[ci] && n->id[pi] < n->id[ci] ))
                    done=1;
                else
                {
                    wt=n->weight[pi];
                    id=n->id[pi];

                    n->weight[pi]=n->weight[ci];
                    n->id[pi]=n->id[ci];

                    n->weight[ci]=wt;
                    n->id[ci]=id;
                }

                if(pi==0)
                    break;

                ci=pi;
            }
        }
        
        n->curSize++;
    }

}

void bweightBPQ(int n,int m, int *ver,int *edges, double* weight, int b, Node* S, int *nlocks, double * weight1) 
{

// Start of matching algorithm
    int i,js;

/*#ifdef DYN
    #pragma omp for schedule(dynamic, CHUNK) private(i)
#else
    #pragma omp for schedule(static, CHUNK) private(i)
#endif*/

    double t1,t2,t3,t4,pr=0.0;
    //b=10;
    //Node *S;
    #pragma omp master
    t1=omp_get_wtime();
    
    #pragma omp for schedule(static) private(i)
    for(i=0;i<=n;i++)  
    {
        weight1[i]=0.0;
        nlocks[i]=0; // Initialize locks
        S[i].curSize=0;
        S[i].maxSize=b;
    }
  
//  count = 0;

  
/*#ifdef DYN
    #pragma omp for schedule(dynamic, CHUNK) private(i)
#else
    #pragma omp for schedule(static, CHUNK) private(i)
#endif*/

    #pragma omp master
    t2=omp_get_wtime();
    
    #pragma omp for schedule(static) private(i)
    for(i=1;i<=n;i++) {
      
      int j,k,x;
      double min_w,t1,diff;
      int min_idx;

       #pragma omp master
       t3=omp_get_wtime();
      
      for(k=1;k<b+1;k++) {
        int done = 0;
        int current = i;
        while (!done) {
            int partner = 0;                  // No point in trying a partner worse than the best suitor
            double heaviest = 0.0;
        

            for(j=ver[current];j<ver[current+1];j++) { // Loop over neighbors of the current vertex 
                int y = edges[j];            // y is the neighbor of the current vertex

                // Check if w(current,y) is the best so far, and if it is a better option for y
                // Also check that current is not already a suitor of y
          
                min_w=min_weight(&S[y]);
                min_idx=min_id(&S[y]);

                //weight1[i]=weight1[i]+(omp_get_wtime()-t1);

                if ((weight[j] < heaviest) || min_w > weight[j])
                    continue;
                if(find_id(&S[y],current))
                    continue;
                if ((weight[j] == heaviest) && (y < partner))
                    continue;

                if (weight[j] == min_w)  
                {
                    int skip = 0;
                    if(__sync_lock_test_and_set(&nlocks[y],1)==0) //Locking Y
                        if ((weight[j] == min_w) && (current < min_idx))            // Must have higher index to beat previous suitor when equal weight
                            skip = 1;
            
                    __sync_lock_release(&nlocks[y]); // Unlocking y
            
                    if (skip)
                        continue;
                }
          
                heaviest = weight[j];      // Store the weight of the heaviest edge found so far
                partner = y;
        } // loop over neighbors

        //cout<<"For ends"<<endl;
        //cout<<"Partner: "<<js<<" "<<partner<<" "<<T[js]<<endl;
        done = 1;
        int next_vertex=0;
        if (heaviest > 0) {            // True if there is a new partner

            if(__sync_lock_test_and_set(&nlocks[partner],1)==0) //Locking partner
            {
                t1=omp_get_wtime();
                next_vertex=min_id(&S[partner]);
                
                // Check if current can become the first suitor of partner
                // Add the priority queue here to handle generalized b candidates
            
                
                //cout<<"JJ: "<<js<<" "<<T[js]<<endl;
                // Update S
                
                Add(&S[partner],heaviest,current);
                weight1[i]=weight1[i]+(omp_get_wtime()-t1);
                
                //t1=omp_get_wtime();
                //cout<<"Assign: "<<next_vertex<<endl;
            }
            //weight1[i]=weight1[i]+(omp_get_wtime()-t1);

          __sync_lock_release(&nlocks[partner]); // Unlocking partner
         
            //cout<<"Make Suitor done"<<endl;
          
          if (next_vertex != 0) {       // True if partner already had a suitor
            //cout<<"Next Vertex: "<<next_vertex<<endl;
            current = next_vertex;      // Pick up the old suitor and continue
            done = false;
          }
          //printf("current = %d, s1(%d) = %d, s2(%d) = %d \n",current,partner,s[partner],partner,s2[partner]); 
        }
      } // while loop
    } // for k=1:2
    #pragma omp master
    t4=omp_get_wtime();
    //cout<<"Iteration "<<i<<": "<<t4-t3<<endl;
  } // loop over vertices

  /*double twt=0.0;
  node t;
  for(i=1;i<=n;i++)
  {
    if(S[i].size()>0 && S[i].size()<3)
    {
        while(!S[i].empty())
        {
            t=S[i].top();
            if(t->weight>0)
                twt+=t->weight;
            else cout<<"WTF .....!!"<<endl;
            S[i].pop();
        }
    }
    else cout<<"Trouble node: "<<i<<endl;
  }
  cout<<"Weight: "<<twt<<endl;*/


}
