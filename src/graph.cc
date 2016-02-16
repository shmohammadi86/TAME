#include <graph.h>

Graph::Graph() {
    adj = NULL;
    n = ne = 0;
}

Graph::Graph(unsigned int n, unsigned int ne, unsigned short int *adj, vector<string> &id2gene) {
	this->n = n;
	this->ne = ne;
	long mem_size = n*n*sizeof(unsigned short int);	
	this->adj = (unsigned short int*)calloc(1, mem_size);
	if(this->adj == NULL) {
		fprintf(stderr, "Graph constructor:: Error allocating memory for adjacency matrix.\n");
		this->n = 0;
		this->ne = 0;
		return;
	}
	memcpy(this->adj, adj, mem_size);
	this->id2gene = id2gene;
	for(register unsigned int i = 0; i < id2gene.size(); i++) {
		this->gene2id[id2gene[i]] = i;
	}
	
	this->deg = (unsigned short int*)calloc(n, sizeof(unsigned short int));
	if(this->deg == NULL) {
		fprintf(stderr, "Graph constructor:: Error allocating memory for degree vector.\n");
		this->n = 0;
		this->ne = 0;
		return;		
	}

	
	int row_idx;
	for(register unsigned int i = 0; i < this->n; i++) {
		row_idx = n*i;
		for(register unsigned int j = 0; j < this->n; j++) {
			if(i == j)
				continue;
				
			if(adj[row_idx+j]) {
				this->deg[i] ++;
			}
		}
	}	
}

void Graph::pruneVertices(vector<int> fixed_set) {
	unsigned register int i, j, current_vertex, removed_edges = 0;
	for(i = 0; i < fixed_set.size(); i++) {
		current_vertex = fixed_set[i];
		this->deg[current_vertex] = 0;
		for(j = 0; j < this->n; j++) {
			if(adj[(n*j) + current_vertex]) {
				removed_edges++;
				adj[(n*j) + current_vertex] = adj[(n*current_vertex) + j] = 0;
			}
		}
	}
	this->ne -= removed_edges;	
	printf("Number of edges after pruning: %u\n", this->ne);
}

void Graph::setEdge(int i, int j, int val) {
	adj[(n*i) + j] = val;
}

int Graph::getEdge(int i, int j) {
	return adj[n*i+j];
}

int Graph::getVertexID(const char *name) {
	map <string, int>::iterator it;
	if( (it = gene2id.find(name)) != gene2id.end()) {
		return it->second;
	}
	else
		return -1;
}

const char * Graph::getVertexName(int v) {
	return (id2gene[v]).c_str();
}


inline void Graph::clear() {
	n = ne = 0;
	free(this->adj);
	this->adj = NULL;
}

int Graph::getVertexDegree(int v) {
	if(v < 0 || v >= (int)n) {
		printf("getVertexDegree:: Invalid vertex id %d\n", v);
		return 0;
	}
	
	return this->deg[v];
}


Graph::~Graph() {
	clear();
}
