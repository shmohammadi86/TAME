#pragma once
#include <generics.h>


class Graph {
public:
    Graph();
    Graph(unsigned int n, unsigned int ne, unsigned short int *adj, vector<string> &id2gene);
	void setEdge(int i, int j, int val);
	int getEdge(int i, int j);
	const char * getVertexName(int v);
	int getVertexID(const char *name);
	int getVertexDegree(int v);
	void pruneVertices(vector<int> fixed_set);
	
	
    ~Graph();
	unsigned int n, ne;	
    
private:
	unsigned short int *adj;	
	unsigned short int *deg;
	
	vector<string> id2gene;	
	map <string, int> gene2id;
	
	void clear();
};
