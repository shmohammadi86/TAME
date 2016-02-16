#pragma once

#include <generics.h>
#include <graph.h>

typedef unsigned int counter_type;

using namespace std;
class TriangleCounts {
public:
	TriangleCounts(Graph *g = NULL);
	
	vector< vector<counter_type> > getTriangles();
	counter_type *getAdjTriangles(int i);
	unsigned int size(); // Number of triangles
	unsigned int pruned_size(); // Number of triangles	, after 0 or more rounds of pruning
	void prune(vector<int> excluded_vertices); // prunes triangles adjacent to matched vertices
	void prune2(bool *isMatched); // prunes triangles adjacent to matched vertices
	counter_type *get_triDegs(); // 1/2 * sum(T(i, :, :)) -> Vector of size num_vertices
	counter_type **get_triEdgeDegs(); // sum(T(:, j, k)) -> Matrix of size num_vertices*num_vertices
	~TriangleCounts();
	
private:
	counter_type *triangle_counts;
	counter_type **edge_triangle_counts;
	
	unsigned short int num_vertices;
	counter_type **triangle_adjacency_list; // List of triangles adjacent to each vertex in form #<j, k>^*, where # (first element) is 2X the number of triangles and the rest of the list are pairs of vertices
	vector< vector<counter_type> > Triangles; // linear list (list of 3-tuples) of triangles
};

