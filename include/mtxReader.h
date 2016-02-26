#ifndef MTXREADER_H
#define MTXREADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <limits>
#include <math.h>

using namespace std;

struct EdgeE
{
    int head;
    int id;         // Edge tail
    float weight;  // Edge weight
};
struct Edge
{
    int id;         // Edge tail
    float weight;  // Edge weight
};

class CSR
{
    public:
    int nVer;       // number of vertices
    int nEdge;      // number of edges
    int maxDeg;
    int* verPtr;    // vertex pointer array of size nVer+1
    Edge* verInd;   // Edge array
    bool bipartite; // general graph or bipartite graph
    int sVer;       // The number of vertices on left for bipartite graph;
    
    bool readMtxG(char * filename); // reading as a general graph
    bool readMtxB(char * filename); // reading as a bipartite graph
    bool mtxG2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool mtxB2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool readCSRbin(char* filename, int opt); // reading binary file of csr format
    bool printDense();
    bool DenseMat2CSR(double** X, int m, int n);
    bool DenseVec2CSR(double* x, int m, int n);

    CSR():nVer(0),nEdge(0),verPtr(NULL),verInd(NULL),bipartite(false){}
    ~CSR()
    {
        if(verPtr!=NULL)
            free(verPtr);

        if(verInd!=NULL)
            free(verInd);
    }

};

#endif //MTXREADER_H
