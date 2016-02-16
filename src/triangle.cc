#include <triangle.h>

/*
 * triangle_counts: For each vertex, counts the number of triangles incident to it, i.e., diag(G^3)) ./ 2;
 * edge_triangle_counts: For each edge, counts how many triangles are incident to it. For each triangle, we need to update all three edges, and to consider symmetry,
 * we should also consider both <i, j> and <j, i>.
 * Triangles: A vector of triangles in form of <i, j, k> tuple.
 * counter: total number of triangles in the graph
 * triangle_adjacency_list: A special datastructure to hold all triangles that are incident to a vertex. The first dimension is index by vertex indices,
 * and the each entry i, which has size of 2*number of adjacent triangles + 1, holds all triangles that are incident to vertex i.
 * The first element is the total number of triangles incident to i * 2, and the rest of the the array is filled with <j, k> such that <i, j, k> is a triangle in G.
 */

TriangleCounts::TriangleCounts(Graph *g) {
	vector<counter_type> tuple(3);	
	register unsigned int i, j, k;
	this->num_vertices = g->n;
	
	
	triangle_counts = (counter_type *)calloc(this->num_vertices, sizeof(counter_type));
	if(triangle_counts == NULL) {
		fprintf(stderr, "TriangleCounts:: Can't allocate memory for triangle_counts\n");
		exit(-1);
	}	
	

	edge_triangle_counts = (counter_type **)calloc(this->num_vertices, sizeof(counter_type*));
	if(edge_triangle_counts == NULL) {
		fprintf(stderr, "TriangleCounts:: Can't allocate memory for edge_triangle_counts\n");
		exit(-1);
	}	
	for(i = 0; i < this->num_vertices; i++) {
		edge_triangle_counts[i] = (counter_type *)calloc(this->num_vertices, sizeof(counter_type));
		if(edge_triangle_counts == NULL) {
			fprintf(stderr, "TriangleCounts:: Can't allocate memory for edge_triangle_counts\n");
			exit(-1);
		}		
	}

	
	for(i = 0; i < this->num_vertices; i++) {
		tuple[0] = i;
		for(j = i+1; j < this->num_vertices; j++) {			
			if(!g->getEdge(i, j))
				continue;
			tuple[1] = j;
			for(k = j+1; k < this->num_vertices; k++) {
				if(g->getEdge(i, k) && g->getEdge(j, k)) {
					tuple[2] = k;
					Triangles.push_back(tuple);
					triangle_counts[i]++; triangle_counts[j]++; triangle_counts[k]++;
					
					edge_triangle_counts[i][j]++; edge_triangle_counts[i][k]++; edge_triangle_counts[j][k] ++;
					edge_triangle_counts[j][i]++; edge_triangle_counts[k][i]++; edge_triangle_counts[k][j] ++;
				}
			}			
		}		
	}
	
	double counter = 0;
	for(i = 0; i < this->num_vertices; i++) {
		counter += triangle_counts[i];
	}	
	
	printf("Mean tries = %f\n", (counter / (double)this->num_vertices));
	
	triangle_adjacency_list = (counter_type **) calloc(g->n, sizeof(counter_type *));
	if(triangle_adjacency_list == NULL) {
		fprintf(stderr, "TriangleCounts:: Can't allocate memory for triangle_adjacency_list\n");
		exit(-1);
	}	
	
	// First entry holds the length, the rest are linearize tuples ( <j1, k1>, <j2, k2> , ... ) and is init to 0 by def.
	int array_size;
	for(i = 0; i < num_vertices; i++) {
		array_size = 2*triangle_counts[i] + 1;
		triangle_adjacency_list[i] = (counter_type *) calloc(array_size, sizeof(counter_type)); 
		if(triangle_adjacency_list[i] == NULL) {
			fprintf(stderr, "TriangleCounts:: Can't allocate memory for triangle_adjacency_list[%d], size = %d\n", i, 2*triangle_counts[i] + 1);
			exit(-1);			
		}
	}
		


	for(register unsigned int index = 0; index < Triangles.size(); index++) {
		i = Triangles[index][0];
		j = Triangles[index][1];
		k = Triangles[index][2];
		
		triangle_adjacency_list[i][++triangle_adjacency_list[i][0]] = j;
		triangle_adjacency_list[i][++triangle_adjacency_list[i][0]] = k;
		
		triangle_adjacency_list[j][++triangle_adjacency_list[j][0]] = i;
		triangle_adjacency_list[j][++triangle_adjacency_list[j][0]] = k;


		triangle_adjacency_list[k][++triangle_adjacency_list[k][0]] = i;
		triangle_adjacency_list[k][++triangle_adjacency_list[k][0]] = j;			
	}	
}

unsigned int TriangleCounts::size() {
	return Triangles.size();
}

unsigned int TriangleCounts::pruned_size() {
	unsigned register int i, sum = 0;
	for(i = 0; i < num_vertices; i++) {
		sum += triangle_adjacency_list[i][0];
	}
	return sum/6;
}

vector< vector<counter_type> > TriangleCounts::getTriangles() {
	return this->Triangles;
}

counter_type **TriangleCounts::get_triEdgeDegs() {
	return this->edge_triangle_counts;	
}


counter_type *TriangleCounts::getAdjTriangles(int i) {
	return this->triangle_adjacency_list[i];
}

counter_type *TriangleCounts::get_triDegs() {
	return triangle_counts;
}

// prunes triangles adjacent to matched vertices
void TriangleCounts::prune(vector<int> excluded_vertices) {		
	printf("Pruning triangles ...\n");
	register unsigned int i, j, k;
	set<int> excluded_vertices_set(excluded_vertices.begin(), excluded_vertices.end());
	
	for(i = 0; i < excluded_vertices.size(); i++) {
		triangle_adjacency_list[excluded_vertices[i]][0] = 0;
	}


	counter_type *N;
	for(i = 0; i < this->num_vertices; i++) {
		N = getAdjTriangles(i);
		k = 1; // Index of first free slot in the current list
		printf("Pruning vertex %d (size = %d)\n", i, N[0]);
		for(j = 1; j <= N[0]; j+=2) { // index of first <j, k> pair to investigate on the tuple list
			printf("\t<%d, %d> ->", N[j], N[j+1]);
			if( (excluded_vertices_set.find(N[j]) != excluded_vertices_set.end()) || (excluded_vertices_set.find(N[j+1]) != excluded_vertices_set.end()) ) {
				printf("pruned\n");
				continue;
			}
			else {
				printf("kept (k = %d)\n", k);
				N[k] = N[j];
				N[k+1] = N[j+1];
				k += 2;
			}			
			
		}
		N[0] = k-1;
	}
}

void TriangleCounts::prune2(bool *isMatched) {		
	register unsigned int i, j, k;
	
	for(i = 0; i < num_vertices; i++) {
		if(isMatched[i]) {
			triangle_adjacency_list[i][0] = 0;			
		}
	}


	counter_type *N;
	for(i = 0; i < this->num_vertices; i++) {
		N = getAdjTriangles(i);
		k = 1; // Index of first free slot in the current list
		for(j = 1; j <= N[0]; j+=2) { // index of first <j, k> pair to investigate on the tuple list
			if(isMatched[N[j]] || isMatched[N[j+1]]) {
				continue;
			}
			else {
				N[k] = N[j];
				N[k+1] = N[j+1];
				k += 2;
			}						
		}
		N[0] = k-1;
	}
}

TriangleCounts::~TriangleCounts() {	
	if(triangle_adjacency_list != NULL) {
		for(int i = 0; i < this->num_vertices; i++) {
			free(triangle_adjacency_list[i]);
		}
		free(triangle_adjacency_list);	
		triangle_adjacency_list = NULL;
	}	

	free(triangle_counts);	
	
	for(int i = 0; i < this->num_vertices; i++) {
		free(edge_triangle_counts[i]);
	}
	free(edge_triangle_counts);	
		
	Triangles.clear();	
}
