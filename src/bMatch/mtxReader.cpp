#include "mtxReader.h"
#include <cstring>
using namespace std;

bool CSR::printDense() {
	register int i, j, nnz = 0;
	//edge_idx = 0, row_idx = 1, total_count = 0;
	for(i = 0; i < sVer; i++) {
		//printf("i = %d\n", i); fflush(stdout);
		for(j = 0; j < nVer-sVer; j++) {
			//printf("\tj = %d (%d)\n", j, verInd[nnz].id); fflush(stdout);
			while(j < verInd[nnz].id && j <  nVer-sVer) {
				printf("0 ");
				j++;
			}
			//printf("j = %d, nnz = %d, vertPtr[i] = %d\n", j, nnz, verPtr[i+1]);
			if(j == nVer-sVer) {
				break;
			}
			if(nnz == verPtr[i+1]) {
				while(j <  nVer-sVer) {
					printf("0 ");
					j++;
				}
				break;
			}
			printf("%f ", verInd[nnz].weight);
			nnz++;
		}
		printf("\n");
	}

	return true;
}

bool CSR::DenseMat2CSR(double** X, int m, int n) {
	register int i, j;
	long nnz = 0;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			if (i == j) continue;
			if(std::numeric_limits<double>::epsilon() < fabs(X[i][j])) {
				nnz++;
			}
		}
	}

	nVer = m + n;
	sVer = m;
	nEdge = 2*nnz;
    bipartite=true;

    vector<vector<int> > graphCRSIdx(nVer);
    vector<vector<double> > graphCRSVal(nVer);
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			if (i == j) continue;
			if(std::numeric_limits<double>::epsilon() < fabs(X[i][j])) {
	            graphCRSIdx[i].push_back(j+sVer);
	            graphCRSVal[i].push_back(X[i][j]);

	            // Sym
	            graphCRSIdx[j+sVer].push_back(i);
	            graphCRSVal[j+sVer].push_back(X[i][j]);
			}
		}
	}


	int count;
    verPtr=new int[nVer+1];
    verInd=new Edge[nEdge];

    verPtr[0]=0;
    int max=0,offset;
    for(int i=1;i<=nVer;i++)
    {
        offset=graphCRSIdx[i-1].size();
        verPtr[i]=verPtr[i-1]+offset;
        count=verPtr[i-1];
        //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
        for(int j=0;j<offset;j++)
        {
            verInd[count].id=graphCRSIdx[i-1][j];
            verInd[count].weight=graphCRSVal[i-1][j];
            count++;

            //cout<<verInd[count-1]<<" ";
        }
        //cout<<endl;
        if(offset>max)
            max=offset;
    }

    assert(count==nEdge);
    maxDeg=max;

    return true;
}

bool CSR::DenseVec2CSR(double* x, int m, int n) {
	register int i, j, idx;
	long nnz = 0;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			if (i == j) continue;
			idx = i*m + j;
			if(std::numeric_limits<double>::epsilon() < fabs(x[idx])) {
				nnz++;
			}
		}
	}

	nVer = m + n;
	sVer = m;
	nEdge = 2*nnz;
    bipartite=true;

    vector<vector<int> > graphCRSIdx(nVer);
    vector<vector<double> > graphCRSVal(nVer);
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			if (i == j) continue;
			idx = i*m + j;

			if(std::numeric_limits<double>::epsilon() < fabs(x[idx])) {
	            graphCRSIdx[i].push_back(j+sVer);
	            graphCRSVal[i].push_back(x[idx]);

	            // Sym
	            graphCRSIdx[j+sVer].push_back(i);
	            graphCRSVal[j+sVer].push_back(x[idx]);
			}
		}
	}


	int count;
    verPtr=new int[nVer+1];
    verInd=new Edge[nEdge];

    verPtr[0]=0;
    int max=0,offset;
    for(int i=1;i<=nVer;i++)
    {
        offset=graphCRSIdx[i-1].size();
        verPtr[i]=verPtr[i-1]+offset;
        count=verPtr[i-1];
        //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
        for(int j=0;j<offset;j++)
        {
            verInd[count].id=graphCRSIdx[i-1][j];
            verInd[count].weight=graphCRSVal[i-1][j];
            count++;

            //cout<<verInd[count-1]<<" ";
        }
        //cout<<endl;
        if(offset>max)
            max=offset;
    }

    assert(count==nEdge);
    maxDeg=max;

    return true;
}

bool CSR::readMtxB(char* filename)
{
    int count=0,i,j;
    int inp, m1, sym, edgecnt_;
    int numRow, numCol, nonZeros, numEdges;
    double f;
    string s;
    ifstream inf;

    inf.open(filename, ios::in);

    if(inf.is_open())
    {
        size_t found1, found2, found3;
        getline(inf,s);
        found1 = s.find("pattern");
        if (found1 != string::npos)
            m1 = 2;
        else
            m1 = 3;
        found1 = s.find("symmetric");
        found2 = s.find("hermitian");
        found3 = s.find("skew-symmetric");
        if (found1 != string::npos || found2 != string::npos || found3 != string::npos)
            sym = 1;
        else
            sym = 0;
        while(inf.peek()=='%')
            getline(inf,s);
 
        inf>>inp;
        numRow=inp;
        inf>>inp;
        numCol=inp;
        inf>>inp;
        nonZeros=inp;

        bipartite=true;
        nVer=numRow+numCol;
        sVer=numRow;

        count=inp;
        
        vector<vector<int> > graphCRSIdx(nVer);
        vector<vector<double> > graphCRSVal(nVer);
        int diag=0;
        
        while(count>0) 
        {     
            inf>>i; 
            inf>>j;

            j+=sVer; //adjusting for the right hand vertices
            
            if(m1==3) 
                inf>>f; 
            else
                f=1.0;

            if(i==j)
                diag++;
            else
            {
                graphCRSIdx[i-1].push_back(j-1); 
                graphCRSVal[i-1].push_back(f);
                if(sym) 
                {     
                    graphCRSIdx[j-1].push_back(i-1); 
                    graphCRSVal[j-1].push_back(f); 
                } 
            }   
            count--; 
        }     
        inf.close(); 
     
        numEdges = nonZeros;  
        if(sym == 1) //symmetric matrix 
            numEdges = nonZeros*2 - 2*diag; 
     
        nEdge=numEdges;
        
        verPtr=new int[nVer+1];
        verInd=new Edge[nEdge];

        verPtr[0]=0;
        int max=0,offset; 
        for(int i=1;i<=nVer;i++)
        {
            offset=graphCRSIdx[i-1].size();
            verPtr[i]=verPtr[i-1]+offset;
            count=verPtr[i-1];
            //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
            for(int j=0;j<offset;j++)
            {
                verInd[count].id=graphCRSIdx[i-1][j];
                verInd[count].weight=graphCRSVal[i-1][j];
                count++;

                //cout<<verInd[count-1]<<" ";
            }
            //cout<<endl;
            if(offset>max)
                max=offset;
        }
        
        assert(count==nEdge);
        maxDeg=max;

    }
    else return false;
   
   return true;
}

bool CSR::readMtxG(char* filename)
{
    int count=0,i,j;
    int inp, m1, sym, edgecnt_;
    int numRow, numCol, nonZeros, numEdges;
    double f;
    string s;
    ifstream inf;

    inf.open(filename, ios::in);

    if(inf.is_open())
    {
        size_t found1, found2, found3;
        getline(inf,s);
        found1 = s.find("pattern");
        if (found1 != string::npos)
            m1 = 2;
        else
            m1 = 3;
        found1 = s.find("symmetric");
        found2 = s.find("hermitian");
        found3 = s.find("skew-symmetric");
        if (found1 != string::npos || found2 != string::npos || found3 != string::npos)
            sym = 1;
        else
            sym = 0;
        while(inf.peek()=='%')
            getline(inf,s);
 
        inf>>inp;
        numRow=inp;
        inf>>inp;
        numCol=inp;
        inf>>inp;
        nonZeros=inp;

        if(numRow==numCol)
            bipartite=false;
        else
            return false;

        count=inp;
        
        vector<vector<int> > graphCRSIdx(numRow);
        vector<vector<double> > graphCRSVal(numRow);
        int diag=0;
        
        while(count>0) 
        {     
            inf>>i; 
            inf>>j; 
            if(m1==3) 
                inf>>f; 
            else
                f=1.0;

            if(i==j)
                diag++;
            else
            {
                graphCRSIdx[i-1].push_back(j-1); 
                graphCRSVal[i-1].push_back(f);
                if(sym) 
                {     
                    graphCRSIdx[j-1].push_back(i-1); 
                    graphCRSVal[j-1].push_back(f); 
                } 
            }   
            count--; 
        }     
        inf.close(); 
     
        numEdges = nonZeros;  
        if(sym == 1) //symmetric matrix 
            numEdges = nonZeros*2 - 2*diag; 
     
        nVer=numRow; 
        sVer=0;
        nEdge=numEdges;
        
        verPtr=new int[nVer+1];
        verInd=new Edge[nEdge];

        verPtr[0]=0;
        int max=0,offset; 
        for(int i=1;i<=nVer;i++)
        {
            offset=graphCRSIdx[i-1].size();
            verPtr[i]=verPtr[i-1]+offset;
            count=verPtr[i-1];
            //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
            for(int j=0;j<offset;j++)
            {
                verInd[count].id=graphCRSIdx[i-1][j];
                verInd[count].weight=graphCRSVal[i-1][j];
                count++;

                //cout<<verInd[count-1]<<" ";
            }
            //cout<<endl;
            if(offset>max)
                max=offset;
        }
        
        assert(count==nEdge);
        maxDeg=max;

    }
    else return false;
   
   return true;
}


bool CSR::mtxB2csrbin(char* filename, char* outfile)
{
    int count=0,i,j;
    int inp, m1, sym, edgecnt_;
    int numRow, numCol, nonZeros, numEdges;
    double f;
    string s;
    ifstream inf;

    inf.open(filename, ios::in);

    if(inf.is_open())
    {
        size_t found1, found2, found3;
        getline(inf,s);
        found1 = s.find("pattern");
        if (found1 != string::npos)
            m1 = 2;
        else
            m1 = 3;
        found1 = s.find("symmetric");
        found2 = s.find("hermitian");
        found3 = s.find("skew-symmetric");
        if (found1 != string::npos || found2 != string::npos || found3 != string::npos)
            sym = 1;
        else
            sym = 0;
        while(inf.peek()=='%')
            getline(inf,s);
 
        inf>>inp;
        numRow=inp;
        inf>>inp;
        numCol=inp;
        inf>>inp;
        nonZeros=inp;

        bipartite=false;
        nVer=numRow+numCol;
        sVer=numRow;
        count=inp;
 
        vector<vector<int> > graphCRSIdx(nVer);
        vector<vector<double> > graphCRSVal(nVer);
        int diag=0;
	    
        while(count>0) 
        {     
            inf>>i; 
            inf>>j; 
            j+=sVer; // Adjust for right hand side of biaprtite graph

            if(m1==3) 
                inf>>f; 
            else
                f=1.0;

            if(i==j)
                diag++;
            else
            {
                graphCRSIdx[i-1].push_back(j-1); 
                graphCRSVal[i-1].push_back(f);
                if(sym) 
                {     
                    graphCRSIdx[j-1].push_back(i-1); 
                    graphCRSVal[j-1].push_back(f); 
                } 
            }   
            count--; 
        }  
        inf.close(); 
 
     
        numEdges = nonZeros;  
        if(sym == 1) //symmetric matrix 
            numEdges = nonZeros*2 - 2*diag; 
     
        nEdge=numEdges;
    
        int* _verPtr=new int[nVer+1];
        int* _verInd=new int[nEdge];
        double* _verWt=new double[nEdge];

        _verPtr[0]=0;
        int max=0,offset;
        for(int i=1;i<=nVer;i++)
        {
            offset=graphCRSIdx[i-1].size();
            _verPtr[i]=_verPtr[i-1]+offset;
            count=_verPtr[i-1];
            //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
            for(int j=0;j<offset;j++)
            {
                _verInd[count]=graphCRSIdx[i-1][j];
                _verWt[count]=graphCRSVal[i-1][j];
                count++;

                //cout<<verInd[count-1]<<" ";
            }
            //cout<<endl;
            if(offset>max)
                max=offset;
        }
        
        assert(count==nEdge);

        ///////////////// reading and csr format done //////////////
        // So write as binary file
        
        ofstream of;
        of.open(outfile,ios::out|ios::binary);
        of.write((char*)&nVer, sizeof(int));
        of.write((char*)&sVer, sizeof(int));
        of.write((char*)&nEdge, sizeof(int));
        of.write((char*)&max, sizeof(int));
        of.write((char*)&_verPtr[0], sizeof(int) * (nVer+1));
        of.write((char*)&_verInd[0], sizeof(int) * nEdge);
        of.write((char*)&_verWt[0], sizeof(double) * nEdge);
        of.close();
    }
    return true;
}

bool CSR::mtxG2csrbin(char* filename, char* outfile)
{
    int count=0,i,j;
    int inp, m1, sym, edgecnt_;
    int numRow, numCol, nonZeros, numEdges;
    double f;
    string s;
    ifstream inf;

    inf.open(filename, ios::in);

    if(inf.is_open())
    {
        size_t found1, found2, found3;
        getline(inf,s);
        found1 = s.find("pattern");
        if (found1 != string::npos)
            m1 = 2;
        else
            m1 = 3;
        found1 = s.find("symmetric");
        found2 = s.find("hermitian");
        found3 = s.find("skew-symmetric");
        if (found1 != string::npos || found2 != string::npos || found3 != string::npos)
            sym = 1;
        else
            sym = 0;
        while(inf.peek()=='%')
            getline(inf,s);
 
        inf>>inp;
        numRow=inp;
        inf>>inp;
        numCol=inp;
        inf>>inp;
        nonZeros=inp;

        if(numRow==numCol)
            bipartite=false;
        else return false;
 
        count=inp;
 
 
        vector<vector<int> > graphCRSIdx(numRow);
        vector<vector<double> > graphCRSVal(numRow);
        int diag=0;
	    
        while(count>0) 
        {     
            inf>>i; 
            inf>>j; 
            if(m1==3) 
                inf>>f; 
            else
                f=1.0;

            if(i==j)
                diag++;
            else
            {
                graphCRSIdx[i-1].push_back(j-1); 
                graphCRSVal[i-1].push_back(f);
                if(sym) 
                {     
                    graphCRSIdx[j-1].push_back(i-1); 
                    graphCRSVal[j-1].push_back(f); 
                } 
            }   
            count--; 
        }  
        inf.close(); 
 
     
        numEdges = nonZeros;  
        if(sym == 1) //symmetric matrix 
            numEdges = nonZeros*2 - 2*diag; 
     
     
        nVer=numRow; 
        nEdge=numEdges;
    
        int* _verPtr=new int[nVer+1];
        int* _verInd=new int[nEdge];
        double* _verWt=new double[nEdge];

        _verPtr[0]=0;
        int max=0,offset;
        for(int i=1;i<=nVer;i++)
        {
            offset=graphCRSIdx[i-1].size();
            _verPtr[i]=_verPtr[i-1]+offset;
            count=_verPtr[i-1];
            //cout<<i-1<<" "<<verPtr[i-1]<<" "<<verPtr[i]<<": ";
            for(int j=0;j<offset;j++)
            {
                _verInd[count]=graphCRSIdx[i-1][j];
                _verWt[count]=graphCRSVal[i-1][j];
                count++;

                //cout<<verInd[count-1]<<" ";
            }
            //cout<<endl;
            if(offset>max)
                max=offset;
        }
        
        assert(count==nEdge);

        ///////////////// reading and csr format done //////////////
        // So write as binary file
        
        ofstream of;
        of.open(outfile,ios::out|ios::binary);
        of.write((char*)&nVer, sizeof(int));
        of.write((char*)&sVer, sizeof(int));
        of.write((char*)&nEdge, sizeof(int));
        of.write((char*)&max, sizeof(int));
        of.write((char*)&_verPtr[0], sizeof(int) * (nVer+1));
        of.write((char*)&_verInd[0], sizeof(int) * nEdge);
        of.write((char*)&_verWt[0], sizeof(double) * nEdge);
        of.close();
    }
    return true;
}

bool CSR::readCSRbin(char* filename, int opt)
{
    ifstream inf;
    inf.open(filename,ios::in|ios::binary);
    if(inf.is_open())
    {
        inf.read((char*)&nVer,sizeof(unsigned int));
        inf.read((char*)&nEdge,sizeof(unsigned int));
        inf.read((char*)&maxDeg,sizeof(unsigned int)); 
        cout<<"n: "<<nVer<<" m: "<<nEdge<<" del: "<<maxDeg<<endl;
        
        #ifdef NEWDEL
        verPtr=new unsigned int[nVer+1];
        verInd=new Edge[nEdge];
        unsigned int* _verInd=new unsigned int[nEdge];
        float* _verWt=new float[nEdge];
        #else
        verPtr=(int*)malloc((nVer+1)*sizeof(int));
        verInd=(Edge*)malloc(nEdge*sizeof(Edge));
        unsigned int* _verInd=(unsigned int*)malloc(nEdge*sizeof(unsigned int));
        float* _verWt=(float*)malloc(nEdge*sizeof(float));
        #endif

        cout<<"Memory Alloc"<<endl;
        
        
        inf.read((char*)&verPtr[0],sizeof(unsigned int)*(nVer+1));
        inf.read((char*)&_verInd[0],sizeof(unsigned int)*nEdge);
        inf.read((char*)&_verWt[0],sizeof(float)*nEdge);
        inf.close();


        
        
        if(opt==1)/// Save edge weight and degree dist
        {
            unsigned int pos=strrchr(filename, '.')-filename+1;
            filename[pos]='\0';
            char edges[100],degree[100];
            strcpy(edges,filename);
            strcpy(degree,filename);
            strcat(edges,"edges");
            strcat(degree,"degree");
            
            ofstream of1,of2;
            of1.open(degree, ios::out);
            of2.open(edges, ios::out);
            for(unsigned int i=0;i<nVer;i++)
                of1<<(verPtr[i+1]-verPtr[i])<<endl;
            for(unsigned int i=0;i<nEdge;i++)
                of2<<_verWt[i]<<endl;
            of1.close();
            of2.close();
        }
        if(opt==2) /// Do METIS STUFF
        {
            ;
        }

        cout<<"Sizes: "<<sizeof(Edge)<<" "<<sizeof(unsigned int)<<" "<<sizeof(int)<<endl;
        #pragma omp parallel for schedule(dynamic,64)
        for(unsigned int i=0;i<nEdge;i++)
        {
            verInd[i].id=_verInd[i];
            verInd[i].weight=(float)_verWt[i];
            //if(verInd[i].weight>0)
                //verInd[i].ewt=verInd[i].weight/2.0;
            //else verInd[i].ewt=0.0;
        }

        #ifdef NEWDEL
        delete _verInd;
        delete _verWt;
        #else
        free(_verInd);
        free(_verWt);
        #endif
        return true;
    }
    else
    return false;

}
