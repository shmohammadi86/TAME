#pragma once

#include <generics.h>
#include <triangle.h>
#include <io.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include "bMatching.h"
#include "mtxReader.h"
#include "pcl_stack.h"

enum graph_idx {first_graph, second_graph};

#define MAXRAND 1000

struct eigen {
	int flag; // Status of computation
	int its; // Actual number of iterations used
	double lambda; // Trace of lambda (the leading eigen value) in each iteration
	double *x; // Final X, in linearized form
	double *trace;
};

struct stats {
	long T; // Conserved triangles
	long M; // Conserved edges
	unsigned int num_matches;
	double matching_score;
	double locality;
	double correctness;
};

struct alignment {
	vector<int> left_match; // left_match[e] is the vertex in G for the e^th match in the alignment
	vector<int> right_match; // right_match[e] is the vertex in G for the e^th match in the alignment

	int *left_project; //  For every vertex i' in H, it holds the projection of i' into G. In otherwords, it is either a vertex id i in G, if i' is matched, or -1 if it is not matched
	int *right_project; //  For every vertex i in G, it holds the projection of i into H. In otherwords, it is either a vertex id i' in H, if i is matched, or -1 if it is not matched

	vector<int> *PrefG; // For every vertex i' in H, PrefG[i'] is the set of potential matches for i' in G
	vector<int> *PrefH; // For every vertex i in G, PrefH[i] is the set of potential matches for i in H


	unsigned long int  match_no;
	unsigned long int conserved_edges;
	unsigned long int conserved_triangles;
	double seqsim;
	int ortho_count;
	double NSim;

	double expected_tri, expected_NSim;
	double expected_tri_avg, expected_NSim_avg;
	double triWeight, NSimWeight;
};

struct Delta {
	long triangle;
	long edge;
	double seqsim;
	int ortho_count;
	double NSim;
	double connectedNess;

	double score;
};

struct Move { // Encodes even-sized augmenting paths
	unsigned move_no;
	vector<int> m_id; // index of match edges to be swapped
	vector< vector<int> > e; // new edge to replace each

	Delta delta;
};

struct Match_deg {
	int match_id;
	long edge_deg;
	long tri_deg;
	double score_deg;
};

double normalize(double *v_norm, double *v, unsigned long n);
double *initUniformX(long int n);
double *initRandomX(long int n);
	
using namespace std;

class ProdTensor {
public:
	ProdTensor(Graph *G, TriangleCounts* T_G, Graph *H, TriangleCounts* T_H, double *w, Sparsity_Type sparsity_flag, char *output_path, char *prefix);
	void impTTV(double *new_x, double *x);
	eigen* issHOPM(int max_it, double alpha, double beta, double epsilon, double *w, double *x0, int init_type);
	alignment* postprocess(double *x_final, int max_iter, int topoDeg, int seqDeg, double seqSim_threshold);
	void addPref(alignment* align, double *weights, int B);
	
	stats evalX(double *x);
	int InitX(double *x, int init_type, double *w = NULL);

	unsigned int m, n1, n2, n;
	~ProdTensor();
	
private:
	int sparsity_type;
	int bin_type;
	unsigned int subspace_dim;
	double relative_weight;
	vector< vector<int> > prior_alignments;
	double rho;

	void Shift(int shift_type);

	int initMatchingDatasets();
	unsigned int setupBipartiteGraph(double *x_in);
	double norm(double *v, unsigned long n);
	double normalize(double *v_norm, double *v, unsigned long n);


	long countTrianglesUnderAlignment(vector<int> mi, vector<int> mj);
	long countEdgesUnderAlignment(vector<int> mi, vector<int> mj);
	Delta Delta_removeMatch(vector<int> mi, vector<int> mj, unsigned int i);
	Delta Delta_addMatch(vector<int> mi, vector<int> mj, unsigned int i, vector<int> e);
	double 	evaluateMove(Move &new_move, alignment *align, int cycle);
	void copyMove(Move &dst, Move& src);
	void applyMove(Move &best_move, alignment* align);
/*	void computeMatchDeg(vector<int> mi, vector<int> mj);*/

	TriangleCounts *T_G, *T_H;
	Graph *G, *H; 

	eigen *res;
	char prefix[1024], output_path[1024];
	double unmatched_weight;
	unsigned int Round;
	
	vector<int> matched_mi, matched_mj;
/*	vector<int> matching_deg;
	vector<int> matching_TriDeg;*/

	double nonzeros_estimate(vector<double> v);
	double opt_tri_stat(double *x, double *mi, double *mj, vector<int> perm);

	double alpha, beta, sigma;
	vec x_vec, x_hat_vec, w_vec, best_x, mixed_x;
	mat X;
	double *pruned_w;
	vector<int> unmatched_rows, unmatched_cols, unmatched_pairs_row, unmatched_pairs_col;

};

