#pragma once

#include <generics.h>
#include <triangle.h>
#include <io.h>

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

struct block {
	unsigned int size;
	double match_score;
	vector<int> rows, cols, unmatched_rows, unmatched_cols;	
	double *B, *V, *H;
	long tries, total_tries;
	double locality;
};

struct block_stats {
	vector<double> tri_contrib;
	vector<double> w;
};

double normalize(double *v_norm, double *v, unsigned long n);
double *initUniformX(long int n);
double *initRandomX(long int n);
	
class ProdTensor {
public:
	ProdTensor(Graph *G, TriangleCounts* T_G, Graph *H, TriangleCounts* T_H, double *w, Sparsity_Type sparsity_flag, char *output_path, char *prefix);
	void symProdTensorVecMul1(double *new_x, double *x);
	void impTTV(double *new_x, double *x);
	void scaled_impTTV(double *new_x, double *x);
	eigen* issHOPM(int subspace_dim, int max_it, double alpha, double beta, double epsilon, double *w, double *x0, int init_type, int shift_type);
	eigen* MLPR(double alpha, double *w, double *x0, int max_it, double epsilon, bool verbose=false);
	
	block Binarize(double *x_in, int method = 0);
	int TAME(double alpha, double beta, double *w, int max_it, double epsilon, int shift_type, int init_type, int num_rounds, double locality_threshold);
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
	
	void ScaleX(double *x, int n);
	void Orthogonalize(mat M);
	vec combineB(mat B);
	void Project(int k);
	//void Shift(double *x_hat, double *x, double *w, double alpha, double beta, int shift_type);
	void Shift(int shift_type);

	void CombineB(mat B, double *w, int k);

	int initMatchingDatasets();
	unsigned int setupBipartiteGraph(double *x_in);
	long countTrianglesUnderAlignment(vector<int> mi, vector<int> mj);
	double norm(double *v, unsigned long n);
	double normalize(double *v_norm, double *v, unsigned long n);
	block_stats computeWeights(double rel_weight);

	// These two only work on the unmatched pairs in v. Ensures the norm of unmatched portion contributes with given weight to the overal norm
	double unmatched_norm(double *v, unsigned long n);
	double unmatched_normalize(double *v_norm, double *v, unsigned long n, double weight);

	TriangleCounts *T_G, *T_H;
	Graph *G, *H; 

	eigen *res;
	vector<int> unmatched_rows, unmatched_cols, unmatched_pairs_row, unmatched_pairs_col;
	vector<block> Blocks;
	char prefix[1024], output_path[1024];
	double unmatched_weight;
	unsigned int Round;
	
	vector<int> matched_mi, matched_mj;
	
	double nonzeros_estimate(vector<double> v);
	double opt_tri_stat(double *x, double *mi, double *mj, vector<int> perm);
	vec Combine_basis(mat B);
	
	double alpha, beta, sigma;
	mat M, B, X, Q, R, S;
	vec x_vec, x_hat_vec, w_vec;

};

