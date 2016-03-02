#include <tensor.h>
//#define VERBOSE_DEBUG

extern char *InitType_names[];

bool *G_isMatched, *H_isMatched;
double *rows, *cols, *edge_weights, *mi, *mj; // matching aux. vectors

bool triDeg_cmp (Match_deg x, Match_deg y) { return (y.tri_deg < x.tri_deg ); }
bool edgeDeg_cmp (Match_deg x, Match_deg y) { return (y.edge_deg < x.edge_deg ); }
bool scoreDeg_cmp (Match_deg x, Match_deg y) { return (y.score_deg < x.score_deg ); }




/**
 * n the number of nodes
 * m the number of nodes
 * nedges the number of edges
 * vv1 is the source for each of the nedges 
 * vv2 is the target for each of the nedges
 * weight is the weight of each of the nedges
 * out1 is a vector of length at most min(n,m),
 * out2 is a vector of length at most min(n,m),
 * noutedges is the number of out edges
 */

double match(int n, int m, int nedges, double *vv1, double *vv2, double *weight, double *out1, double *out2, int *noutedges) {	
	
/*	double **adj = (double **)calloc(n, sizeof(double *));
	for(register int ii = 0; ii < n; ii++) {
		adj[ii] = (double *) calloc(m, sizeof(double));
	}
	for(register int ii = 0; ii < nedges; ii++) {
		adj[(int)(vv1[ii])][(int)(vv2[ii])] = weight[ii];
	}
	
	char path[1024];
	sprintf(path, "%s/match_%d.log", output_path, ++counter);
	FILE *log = fopen(path, "w");
	for(register int ii = 0; ii < n; ii++) {
		for(register int jj = 0; jj < m; jj++) {
			fprintf(log, "%e\t", adj[ii][jj]);
		}
		fprintf(log, "\n");
	}	
	fclose(log);
*/
	
	double ret, al;
	double *l1, *l2, *w;
	int *match1, *match2, *v1, *v2;
	int i, j, k, p, q, r, t1, t2;
	int *s, *t, *deg, *offset, *list;

	l1 = new double[n];
	l2 = new double[n+m];
	v1 = new int[nedges];
	v2 = new int[nedges];
	s = new int[n];
	t = new int[n+m];
	match1 = new int[n];
	match2 = new int[n+m];
	offset = new int[n];
	deg = new int[n];
	list = new int[nedges + n];
	w = new double[nedges + n];
		

	for (i = 0; i < nedges; i++) {
		v1[i] = (int)(vv1[i] + .5);
		v2[i] = (int)(vv2[i] + .5);
	}
	for (i = 0; i < n; i++) {
		offset[i] = 0;
		deg[i] = 1;
	}
	for (i = 0; i < nedges; i++) deg[v1[i]]++;
	for (i = 1; i < n; i++) offset[i] = offset[i-1] + deg[i-1];
	for (i = 0; i < n; i++) deg[i] = 0;
	for (i = 0; i < nedges; i++) {
		list[offset[v1[i]] + deg[v1[i]]] = v2[i];
		w[offset[v1[i]] + deg[v1[i]]] = weight[i];
		deg[(int)v1[i]]++;
	}
	for (i = 0; i < n; i++) {
		list[offset[i] + deg[i]] = m + i;
		w[offset[i] + deg[i]] = 0;
		deg[i]++;
	}
	for (i = 0; i < n; i++) {
		l1[i] = 0;
		for (j = 0; j < deg[i]; j++) {
			if (w[offset[i]+j] > l1[i]) l1[i] = w[offset[i] + j];
		}
	}
	for (i = 0; i < n; i++) {
		match1[i] = -1;
	}
	for (i = 0; i < n + m; i++) {
		l2[i] = 0;
		match2[i] = -1;
	}
	for (i = 0; i < n; i++) {
		for(j = 0; j < n + m; j++) t[j] = -1;
		s[p = q = 0] = i;
		for(; p <= q; p++) {
			if (match1[i] >= 0) break;
			k = s[p];
			for (r = 0; r < deg[k]; r++) {
				if (match1[i] >= 0) break;
				j = list[offset[k] + r];
				if (w[offset[k] + r] < l1[k] + l2[j] - 1e-8) continue;
				if (t[j] < 0) {
					s[++q] = match2[j];
					t[j] = k;
					if (match2[j] < 0) {
						for(; j>=0 ;) {
							k = match2[j] = t[j];
							p = match1[k];
							match1[k] = j;
							j = p;
						}
					}
				}
			}
		}
		if (match1[i] < 0) {
			al = 1e20;
			for (j = 0; j < p; j++) {
				t1 = s[j];
				for (k = 0; k < deg[t1]; k++) {
					t2 = list[offset[t1] + k];
					if (t[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al) {
						al = l1[t1] + l2[t2] - w[offset[t1] + k];
					}
				}
			}
			for (j = 0; j < p; j++) l1[s[j]] -= al;
			for (j = 0; j < n + m; j++) if (t[j] >= 0) l2[j] += al;
			i--;
			continue;
		}
	}
	ret = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < deg[i]; j++) {
			if (list[offset[i] + j] == match1[i]) {
				ret += w[offset[i] + j];
			}
		}
	}
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) (*noutedges)++;
    }
    *noutedges = 0;
    for (i = 0; i < n; i++) {
        if (match1[i] < m) {
            out1[*noutedges] = i;
            out2[*noutedges] = match1[i];
            (*noutedges)++;
        }
    }

    delete [] l1;
    delete [] l2;
    delete [] v1;
    delete [] v2;
    delete [] s;
    delete [] t;
	delete [] match1;
	delete [] match2;
    delete [] offset;
    delete [] deg;
    delete [] list;
    delete [] w;	
    
	return ret;
}

double ProdTensor::norm(double *v, unsigned long n) { 
	unsigned register long int i;  
	register double sqsum = 0;
	
	for (i=0; i < n; i++) {
		sqsum += (v[i]*v[i]);
	}
    
    return sqrt(sqsum);
}

/* v_norm = v / norm(v, 2);
* Input) v: a vector of size n, n: len of vector v
* Output) v_norm
*/
double ProdTensor::normalize(double *v_norm, double *v, unsigned long n) { 
    double Norm =  norm(v, n);
    if(Norm < 1e-52) 
		return 0;		
    double scale = 1 / Norm;    

	for (register unsigned int i = 0; i < n; i++) {
		v_norm[i] = v[i]*scale;
	}
	
	return Norm;
}




int ProdTensor::initMatchingDatasets() {
	rows = (double *)calloc(n, sizeof(double));
	if(rows == NULL) {
		fprintf(stderr, "ProdTensor::initMatchingDatasets Can't allocate memory for 'rows'\n");
		return -3;
	}
	
	
	cols = (double *)calloc(n, sizeof(double));
	if(cols == NULL) {
		fprintf(stderr, "ProdTensor::initMatchingDatasets Can't allocate memory for 'cols'\n");
		return -4;
	}

	edge_weights = (double *)calloc(n, sizeof(double));
	if(edge_weights == NULL) {
		fprintf(stderr, "ProdTensor::initMatchingDatasets Can't allocate memory for 'edge_weights'\n");
		return -5;
	}
		
    
	mi = (double *)calloc(min(n1, n2), sizeof(double));
	if(mi == NULL) {
		fprintf(stderr, "ProdTensor::initMatchingDatasets Can't allocate memory for 'mi'\n");
		return -6;
	}
	
	mj = (double *)calloc(min(n1, n2), sizeof(double));
	if(mj == NULL) {
		fprintf(stderr, "ProdTensor::initMatchingDatasets Can't allocate memory for 'mj'\n");
		return -7;
	}	
	
	return 0;
}


ProdTensor::ProdTensor(Graph *G, TriangleCounts *T_G, Graph *H, TriangleCounts *T_H, double *w, Sparsity_Type sparsity_type, char *output_path, char *prefix) {
	//tensor parameters
	this->m = 3;
	this->G = G;
	this->H = H;
	this->T_G = T_G;
	this->T_H = T_H;

	this->n1 = this->G->n;
	this->n2 = this->H->n;
	this->n = n1*n2; /* number of edges in L */

	memcpy(this->prefix, prefix, 1024);	
	memcpy(this->output_path, output_path, 1024);	

	if(initMatchingDatasets() != 0) {
		fprintf(stderr, "ProdTensor:: Can't initialize matching datasets.\n");
		exit(-1);
	}

	this->sparsity_type = sparsity_type;
	printf("Sparsity flag = %d\n", sparsity_type);
	unsigned register int i, j;

	// WARNING! This has to be set even in sparse case since we need it currently. This feature is developed for fixing "known" alignments, but needs revision
	unmatched_rows.resize(n1);
	for(i = 0; i < n1; i++) {
		unmatched_rows[i] = i;
	}
	unmatched_cols.resize(n2);
	for(i = 0; i < n2; i++) {
		unmatched_cols[i] = i;
	}


	if( this->sparsity_type != NoSparsity ) {
		int total_matches = 0;
		for(i = 0; i < n1; i++) {
			for(j = 0; j < n2; j++) {
				if( w[i*n2 + j] != 0 ) {
					total_matches ++;
					unmatched_pairs_row.push_back(i);
					unmatched_pairs_col.push_back(j);
				}
			}
		}
		printf("Total matches in L = %d, size unmatched_pairs = %d\n", total_matches, (int)unmatched_pairs_row.size());
	}
}


void ProdTensor::impTTV(double *new_x, double *x) {

	if( this->sparsity_type == NoSparsity ) {
		unsigned int *Ni1, *Ni2, *p, *q;
		register unsigned int i1, j1, k1, i2, j2, k2, idx, cur_row, cur_col;

		// Ensure fixed submatrices are untouched.
		memcpy(new_x, x, n*sizeof(double));


		for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
			i1 = unmatched_rows[cur_row];
			Ni1 = T_G->getAdjTriangles(i1);
			for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {
				i2 = unmatched_cols[cur_col];
				Ni2 = T_H->getAdjTriangles(i2);
				idx = n2*i1 + i2; // X[i1][i2];
				new_x[idx] = 0; // resets new_x[idx]
				for(p = Ni1+1; p <= Ni1 + *Ni1; ) {
					j1 = *p++; k1 = *p++;
					for(q = Ni2+1; q <= Ni2 + *Ni2; ) {
						j2 = *q++; k2 = *q++;
						new_x[idx] += (x[n2*j1 + j2] * x[n2*k1 + k2]) + (x[n2*j1 + k2] * x[n2*k1 + j2]); // (m-1)!*( (m-1)! - 1 ) non-redundant combinations
					}
				}
				new_x[idx] *= 2; // To account for the symmetry in T, i.e., X[j1][j2]*X[k1][k2] vs X[k1][k2]*X[j1][j2], (m-1)! redundancy in permutations.
			}
		}
	}
	else {
		bzero(new_x, this->n*sizeof(double));
		unsigned int *Ni1, *Ni2, *p, *q;

		register unsigned int cur_pair, i1, i2, j1, j2, k1, k2, idx;
		for(cur_pair = 0; cur_pair < unmatched_pairs_row.size(); cur_pair++) {
			i1 = unmatched_pairs_row[cur_pair];
			i2 = unmatched_pairs_col[cur_pair];
			idx = n2*i1 + i2; // X[i1][i2];

			Ni1 = T_G->getAdjTriangles(i1);
			Ni2 = T_H->getAdjTriangles(i2);

			for(p = Ni1+1; p <= Ni1 + *Ni1; ) {
				j1 = *p++; k1 = *p++;
				for(q = Ni2+1; q <= Ni2 + *Ni2; ) {
					j2 = *q++; k2 = *q++;
					new_x[idx] += (x[n2*j1 + j2] * x[n2*k1 + k2]) + (x[n2*j1 + k2] * x[n2*k1 + j2]); // (m-1)!*( (m-1)! - 1 ) non-redundant combinations
				}
			}
			new_x[idx] *= 2; // To account for the symmetry in T, i.e., X[j1][j2]*X[k1][k2] vs X[k1][k2]*X[j1][j2], (m-1)! redundancy in permutations.
		}
	}
}

// x is row-order linearized version of matrix X
stats ProdTensor::evalX(double *x) {
	stats myStats;
	int noutedges;

	unsigned register int j, k, l;


	double v1 = 0, v2 = 0, v4 = 0, temp;
	for(j = 0; j < (unsigned int)n; j++) {
		temp = x[j];
		v1 += temp;
		temp *= temp;
		v2 += temp;
		temp *= temp;
		v4 += temp;
	}
	myStats.locality = round(v2*v2/v4);
	

	setupBipartiteGraph(x);
	myStats.matching_score = match(n1, n2, unmatched_rows.size()*unmatched_cols.size(), rows, cols, edge_weights, mi, mj, &noutedges);	
	myStats.num_matches = noutedges;

	myStats.correctness = 0;
	int common_subgraph_size = 2000;
	for(j = 0; j < (unsigned)noutedges; j++) {
		if(mi[j] < common_subgraph_size && mj[j] < common_subgraph_size && mi[j] == mj[j])
			myStats.correctness ++;
	}
	myStats.correctness = 100 * myStats.correctness / common_subgraph_size;

	wall_clock timer;
	timer.tic();

	// Construct alignment graph
	bool **alignGraph = (bool **) calloc(noutedges, sizeof(bool *));
	if(alignGraph == NULL) {
		fprintf(stderr, "\t\tissHOPM:: Can't allocate memory for alignGraph.\n");
		exit(-1);
	}
	for(j = 0; j < (unsigned)noutedges; j++) {
		alignGraph[j] = (bool *) calloc(noutedges, sizeof(bool));
		if(alignGraph[j] == NULL) {
			fprintf(stderr, "\t\tissHOPM:: Can't allocate memory for alignGraph[%d].\n", j);
			exit(-1);
		}
	}

	register unsigned int i1, i2, j1, j2;
	for(j = 0; j < (unsigned)noutedges; j++) {
		i1 = mi[j];
		i2 = mj[j];		
		for(k = 0; k < (unsigned)noutedges; k++) {
			j1 = mi[k];
			j2 = mj[k];		
			if(G->getEdge(i1, j1) && H->getEdge(i2, j2)) {
				alignGraph[j][k] = true;
			}
		}		
	}

	timer.tic();
	// Count conserved triangles/edges
	myStats.T = myStats.M = 0;
	for(j = 0; j < (unsigned)noutedges; j++) {
		for(k = j+1; k < (unsigned)noutedges; k++) {
			if( !alignGraph[j][k] )
				continue;	
			myStats.M++;				
			for(l = k+1; l < (unsigned)noutedges; l++) {
				if( alignGraph[j][l] && alignGraph[k][l] ) {
						myStats.T++;
				}										
			}
		}
	}							

	return myStats;	
}

//void ProdTensor::Shift(double *x_hat, double *x, double *w, double alpha, double beta, int shift_type) {
void ProdTensor::Shift(int shift_type) {

/*
	unsigned register int j;
	switch(shift_type) {
		case(VARIABLE_SHIFT):
			// Shift using x (inertia idea)
			for(j = 0; j < n; j++) {
				x_hat[j] += alpha*x[j];
			}
			break;
		case(AFFINE_SHIFT):
			// Shift using SeqSim
			for(j = 0; j < n; j++) {
				x_hat[j] += (beta*this->rho*w[j]);
			}
			break;
		case(NORM_AFFINE_SHIFT):
			// Normalize first:: x_hat <-- ( x_hat / norm(x_hat, 2) )
			if(normalize(x_hat, x_hat, n)  == 0)  {
			  fprintf(stderr, "issHOPM:: x_hat vector has norm 0 after first norm.\n");
			  res->flag = -1;
			  return;
			}
			// Then shift using SeqSim
			for(j = 0; j < n; j++) {
				x_hat[j] = beta*x_hat[j] + (1-beta)*w[j];
			}
			break;
		case(COMBINED_SHIFT):
			for(j = 0; j < n; j++) {
				x_hat[j] += ((1 - beta)*arma::norm(x_hat_vec)*w[j] + alpha*x[j]);
			}
			break;
	}
	*/
	printf("\t\t\t\tShift:: Shift type = %d, alpha = %.2e, beta=%.2e, sigma = %.2e, rho =%.2f\n", shift_type, alpha, beta, sigma, rho);
	printf("\t\t\t\tShift:: ||x_hat_vec|| = %.2f, Effective x_vec = %.2e, Effective weight w = %.2e\n", arma::norm(x_hat_vec), this->alpha / arma::norm(x_hat_vec), this->sigma / arma::norm(x_hat_vec));

/*	register int i;
	vec x_norm = x_hat_vec;
	double x_hat_sum = 0, w_sum = 0;
	for (i = 0; i < 2000; i++) {
		x_hat_sum += x_norm(i*n2+i);
		w_sum += w_vec(i*n2+i);
	}
	printf("\t\t\t\tSum(diag(w)) = %e, sum(diag(x_norm)) = %e\n", x_hat_sum, w_sum);
*/

	switch(shift_type) {
		case(VARIABLE_SHIFT):
			x_hat_vec = x_vec + this->alpha*x_vec;
			break;
		case(AFFINE_SHIFT):
			// Shift using SeqSim
			x_hat_vec = x_hat_vec + this->sigma*this->rho*w_vec;
			break;
		case(NORM_AFFINE_SHIFT):
			x_hat_vec = (x_hat_vec / arma::norm(x_hat_vec)) + this->sigma*w_vec;
			break;
		case(COMBINED_SHIFT):
			if(this->beta == 0) // Topological weight = 0 -> Only SeqSim
				x_hat_vec = w_vec + this->alpha*x_vec;
			else if(this->beta == 1) // Topological weight = 1 -> SS-HOPM
				x_hat_vec = x_hat_vec + this->alpha*x_vec;
			else {
				x_hat_vec = x_hat_vec + this->alpha*x_vec + this->sigma*w_vec;//this->sigma*arma::norm(x_hat_vec)*w_vec;
				/*
				x_hat_vec = x_hat_vec + this->alpha*x_vec;
				x_hat_vec = x_hat_vec % (this->sigma*w_vec);
				*/
			}
			break;
	}
}


void ProdTensor::addPref(alignment* align, double *weights, int B) {
	register unsigned int i, j;
	printf("Running bMatching to identify Pref sets\n");

	/************************
	 *    Allocate memory
	 ***********************/
	printf("\tConstructing initial graph for bMatching\n"); fflush(stdout);
    CSR G;
    G.DenseVec2CSR(weights, this->n1, this->n2);

    int* nlocks=(int*)malloc(G.nVer*sizeof(int));    //required for all schemes
    int* start=(int*)malloc(G.nVer*sizeof(int));     //required for sorted and part. sorted
    int* end=(int*)malloc(G.nVer*sizeof(int));       //required for part sorted

    int *b=(int*)malloc(G.nVer*sizeof(int));
    for(i = 0; i < (unsigned int)G.nVer; i++) {
		int card = 0;
		for(j = G.verPtr[i]; j < (unsigned int)G.verPtr[i+1]; j++)
			if(0 < G.verInd[j].weight)
				card++;

		if(B < card)
			b[i]=B;
		else b[i]=max(1, card);
    }

    Node* S=(Node*)malloc(G.nVer*sizeof(Node));      //Heap data structure
	for(i = 0; i < (unsigned int)G.nVer; i++) {
		S[i].heap=(Info*)malloc(max(1, b[i])*sizeof(Info));      //Each heap of size b
		S[i].curSize = 0;
	}


	/***********************
	 *    Perform bMatching
	 ***********************/
	localDom(&G, S, b, start, end, nlocks, 3, 2, 0);


	/***********************************
	 *    Add bMatches to the Pref Sets
	 ***********************************/
	switchCandidate cand;

	printf("\tAdding bMatches to H preferred set\n"); fflush(stdout);
	for (i = 0; i < this->n1; i++) {
/*		printf("\t\t%s:\n", this->G->getVertexName(i));*/
		for(j = 0; j < (unsigned int) S[i].curSize; j++) {
			cand.vertex_id = (S[i].heap[j].id) - this->n1;
/*			printf("\t\t\t-> %s\n", this->H->getVertexName(cand.vertex_id));*/

			align->PrefH[i].push_back( cand );
		}
	}


	printf("\tAdding bMatches to G preferred set\n"); fflush(stdout);
	for (i = 0; i < this->n2; i++) {
/*		printf("\t\t%s:\n", this->H->getVertexName(i));*/

		for(j = 0; j < (unsigned int) S[i+this->n1].curSize; j++) {
			cand.vertex_id = S[i+this->n1].heap[j].id;
/*			printf("\t\t\t-> %s\n", this->G->getVertexName(cand.vertex_id));*/

			align->PrefG[i].push_back( cand );
		}
	}




/*
	for(i = 0; i < align->match_no; i++) {
		printf("Mi[%d] =  %d (%s) in G\n", i, align->left_match[i], this->G->getVertexName(align->left_match[i]));
		for(j = 0; j < align->PrefH[align->left_match[i]].size(); j++) {
			printf("\tVertex %d (%s) in H\n", align->PrefH[align->left_match[i]][j].vertex_id, this->H->getVertexName(align->PrefH[align->left_match[i]][j].vertex_id));
		}
	}


	for(i = 0; i < align->match_no; i++) {
		printf("Mj[%d] =  %d (%s) in H\n", i, align->right_match[i], this->H->getVertexName(align->right_match[i]));
		for(j = 0; j < align->PrefG[align->right_match[i]].size(); j++) {
			printf("\tVertex %d (%s) in G\n", align->PrefG[align->right_match[i]][j].vertex_id, this->G->getVertexName(align->PrefG[align->right_match[i]][j].vertex_id));
		}
	}
*/

	/**************************
	 *    Deallocate memory
	 **************************/
    free(nlocks);
    free(start);
    free(end);
    free(b);
	for(i = 0; i < (unsigned int)G.nVer; i++) {
		free(S[i].heap);
	}
    free(S);
}

alignment* ProdTensor::postprocess(double *x_final, int max_iter, int topoDeg, int seqDeg, double seqSim_threshold) {
	printf("Post-processing ... \n");

	register unsigned int i, i_prime, j, j_prime, k, l, s, it;
	alignment *align = new alignment;


	/*******************************************************
	 *    Initialize alignment with MWM over TAME matrix
	 *******************************************************/
	pruned_w = new double[this->n];
	memcpy(pruned_w, w_vec.memptr(), this->n*sizeof(double));
	for(i = 0; i < this->n; i++) {
		if(pruned_w[i] < seqSim_threshold)
			pruned_w[i] = 0;
	}

	double *temp_x = new double[this->n];
	memcpy(temp_x, x_final, this->n*sizeof(double));
	stats myStats = evalX(temp_x);

	align->conserved_edges = myStats.M;
	align->conserved_triangles = myStats.T;
	align->match_no = myStats.num_matches;


	align->left_match.resize(align->match_no);
	align->right_match.resize(align->match_no);

	align->PrefH = new vector<switchCandidate>[this->n1];
	align->right_project = new int[this->n1];
	for(i = 0; i < this->n1; i++) {
		align->right_project[i] = -1;
	}

	align->PrefG = new vector<switchCandidate>[this->n2];
	align->left_project = new int[this->n2];
	for(i = 0; i < this->n2; i++) {
		align->left_project[i] = -1;
	}

	switchCandidate cand;

	for (i = 0; i < align->match_no; i++) {
		align->left_match[i] = (int)mi[i];
		align->right_match[i] = (int)mj[i];


		cand.vertex_id = align->left_match[i];
		align->PrefG[align->right_match[i]].push_back(cand);


		cand.vertex_id = align->right_match[i];
		align->PrefH[align->left_match[i]].push_back(cand);


		align->left_project[align->right_match[i]] = align->left_match[i];
		align->right_project[align->left_match[i]] = align->right_match[i];

		temp_x[align->left_match[i]*this->n2+align->right_match[i]] = 0;
	}

	align->total_seqSim = align->ortho_count = 0;
	for (l = 0; l < align->match_no; l++) {
		align->total_seqSim += pruned_w[this->n2 * align->left_match[l] + align->right_match[l]];
		align->ortho_count += seqSim_threshold <= pruned_w[this->n2 * align->left_match[l] + align->right_match[l]]? 1:0;
	}
	printf("\tInitial matches = %ld, edges = %ld, triangles = %ld, total SeqSim = %f (ortho count = %d)\n", align->match_no, align->conserved_edges, align->conserved_triangles, align->total_seqSim, align->ortho_count);


	/****************************************
	 *    Add Candidates to the Pref Sets
	 ****************************************/
	// Add (topoDeg-1)-matches from Topo
	addPref(align, temp_x, topoDeg-1);

	// Add seqDeg-matches from pruned seqSim
	addPref(align, pruned_w, seqDeg);



	/******************************************
	 *   Process each candidate individually
	******************************************/
	Move new_move, best_move;

	unsigned int max_length = 2; // We only consider moves of size up to 2 (size 4 paths)
	new_move.m_id.resize(max_length);
	new_move.e.resize(max_length);
	for(i = 0; i < max_length; i++) {
		new_move.e[i].resize(2); // reserve size for end points of matching edges
	}

	best_move.m_id.resize(max_length);
	best_move.e.resize(max_length);
	for(i = 0; i < max_length; i++) {
		best_move.e[i].resize(2); // reserve size for end points of matching edges
	}

	double total_improvement = 0;
	long Delta_tri = 0, Delta_edge = 0, Delta_seqsim = 0;

	vector<Match_deg> match_degree;
	match_degree.resize(align->match_no);
	Delta curr_match_delta;
	double row_sum, col_sum;
	vector<int> right_unmatched, left_unmatched, matched_candidates;
	for(it = 0; it < (unsigned int)max_iter; it++) {
		printf("Iteration %d \n", it);
		total_improvement = 0;
		Delta_tri = Delta_edge = Delta_seqsim = 0;

		/*********************************************
		 *   Decide in which order to process matches
		*********************************************/
		for(k = 0; k < align->match_no; k++) {
			match_degree[k].match_id = k;
			curr_match_delta = Delta_removeMatch(align->left_match, align->right_match, k);
			match_degree[k].edge_deg = curr_match_delta.edge;
			match_degree[k].tri_deg = curr_match_delta.triangle;
			row_sum = col_sum = 0;

			for(l = 0; l < this->n1; l++) {
				col_sum += x_final[this->n2*l + align->right_match[k]];
			}
			for(l = 0; l < this->n2; l++) {
				row_sum += x_final[this->n2*align->left_match[k] + l];
			}

			match_degree[k].score_deg = (row_sum/this->n2 + col_sum/this->n1);
		}
/*		stable_sort(match_degree.begin(), match_degree.end(), triDeg_cmp);*/
/*		stable_sort(match_degree.begin(), match_degree.end(), edgeDeg_cmp);*/
		stable_sort(match_degree.begin(), match_degree.end(), scoreDeg_cmp);


		/*********************************************
		 *   Process each match in the given order
		*********************************************/
		for(unsigned int d = 0; d < align->match_no; d++) {
			k = match_degree[d].match_id;
			best_move.move_no = 0;
			bzero(&(best_move.delta), sizeof(Delta));
			i = align->left_match[k]; i_prime = align->right_match[k];

			left_unmatched.clear(); // Nodes on the left side (G) that are in Pref set of i' and are either matched/unmatched
 			for(l = 0; l < align->PrefG[i_prime].size(); l++) {
				j = align->PrefG[i_prime][l].vertex_id;
				if(align->right_project[j] == -1) {
					left_unmatched.push_back(j);
				}
			}

			right_unmatched.clear();
			for(l = 0; l < align->PrefH[i].size(); l++) {
				j_prime = align->PrefH[i][l].vertex_id;
				if(align->left_project[j_prime] == -1) {
					right_unmatched.push_back(j_prime);
				}
			}

			new_move.move_no = 1;
			new_move.m_id[0] = k;

			new_move.e[0][1] = i_prime;
			for(l = 0; l < left_unmatched.size(); l++) {
				j = left_unmatched[l];
				new_move.e[0][0] = j;

				//printf("Changing %d <-> %d to %d <-> %d\n", i, i_prime, new_move.e[0][0], new_move.e[0][1]); fflush(stdout);
				evaluateMove(new_move, align);
				if(best_move.delta.score < new_move.delta.score ) {
					copyMove(best_move, new_move);
				}
			}

			new_move.e[0][0] = i;
			for(l = 0; l < right_unmatched.size(); l++) {
				j_prime = right_unmatched[l];
				new_move.e[0][1] = j_prime;

				//printf("Changing %d <-> %d to %d <-> %d\n", i, i_prime, new_move.e[0][0], new_move.e[0][1]); fflush(stdout);
				evaluateMove(new_move, align);
				if(best_move.delta.score < new_move.delta.score ) {
					copyMove(best_move, new_move);
				}
			}

			// Now look at already matched edges and evaluate them as swap candidates
			matched_candidates.clear();
			new_move.move_no = 2;
			for(l = 0; l < align->match_no; l++) {
				if(l == k)
					continue;
				j = align->left_match[l];
				j_prime = align->right_match[l];
				new_move.m_id[1] = l;

				new_move.e[0][0] = i;
				new_move.e[0][1] = j_prime;
				new_move.e[1][0] = j;
				new_move.e[1][1] = i_prime;

				// Prune based on preferred sets
				// 1) Make sure j' \in PrefH(i)
				for(s = 0; s < align->PrefH[i].size(); s++) {
					if(align->PrefH[i][s].vertex_id == (int)j_prime)
						break;
				}
				if(s == align->PrefH[i].size())
					continue;

				// 1) Make sure j \in PrefG(i')
				for(s = 0; s < align->PrefG[i_prime].size(); s++) {
					if(align->PrefG[i_prime][s].vertex_id == (int)j)
						break;
				}
				if(s == align->PrefG[i].size())
					continue;

				// prune based on seqsim, if needed

				// Prune based on edges

				//printf("Swap %d <-> %d and %d <-> %d with %d <-> %d and %d <-> %d\n", i, i_prime, j, j_prime, new_move.e[0][0], new_move.e[0][1], new_move.e[1][0], new_move.e[1][1]); fflush(stdout);
				evaluateMove(new_move, align);
				if(best_move.delta.score < new_move.delta.score ) {
					copyMove(best_move, new_move);
				}
			}
			if(0 < best_move.delta.score) {
				applyMove(best_move, align);
				total_improvement += best_move.delta.score;
				Delta_edge += best_move.delta.edge;
				Delta_tri += best_move.delta.triangle;
				Delta_seqsim += best_move.delta.seqsim;
			}
		}
		if(0 < total_improvement) {
			printf("Total improvement = %lf\n", total_improvement);
			printf("\tDelta SeqSim = %ld\n", Delta_seqsim);
			printf("\tDelta Edges = %ld\n", Delta_edge);
			printf("\tDelta Triangles = %ld\n", Delta_tri);

			printf("\tTotal Edges = %ld\n", countEdgesUnderAlignment(align->left_match, align->right_match));
			printf("\tTotal Triangles = %ld\n", countTrianglesUnderAlignment(align->left_match, align->right_match));
			int ortho_count = 0;
			for(l = 0; l < align->match_no; l++) {
				if(pruned_w[this->n2*align->left_match[l] + align->right_match[l]] > seqSim_threshold) {
					ortho_count ++;
				}
			}
			printf("\tTotal orthologs = %d\n", ortho_count);

		}
		else {
			printf("No more improvement during a whole iteration. Exiting now.");
			break;
		}

	}
	align->conserved_edges = countEdgesUnderAlignment(align->left_match, align->right_match);
	align->conserved_triangles = countTrianglesUnderAlignment(align->left_match, align->right_match);
	align->total_seqSim = 0;
	for (l = 0; l < align->match_no; l++) {
		align->total_seqSim += pruned_w[this->n2 * align->left_match[l] + align->right_match[l]];
	}

	return align;
}

double ProdTensor::evaluateMove(Move &new_move, alignment* align) {
/*	int i = align->left_match[new_move.m_id[0]];
	int j = align->left_match[new_move.m_id[1]];
	int i_prime = align->right_match[new_move.m_id[0]];
	int j_prime = align->right_match[new_move.m_id[1]];
	if(new_move.move_no == 1) {
		printf("\tEvaluating %d <-> %d versus %d <-> %d ... \n", i, i_prime, new_move.e[0][0], new_move.e[0][1]); fflush(stdout);
	}
	else if(new_move.move_no == 2) {
		printf("\tEvaluating %d <-> %d and %d <-> %d versus %d <-> %d and %d <-> %d ... \n", i, i_prime, j, j_prime, new_move.e[0][0], new_move.e[0][1], new_move.e[1][0], new_move.e[1][1]); fflush(stdout);
	}
*/

	vector<int> new_mi = align->left_match;
	vector<int> new_mj = align->right_match;
	Delta delta_loss, delta_gain;
	bzero(&(new_move.delta), sizeof(Delta));
	unsigned register int k;

	if(new_move.move_no == 1) {
		delta_loss = Delta_removeMatch(align->left_match, align->right_match, new_move.m_id[0]);
		delta_gain = Delta_addMatch(align->left_match, align->right_match, new_move.m_id[0], new_move.e[0]);

		new_move.delta.triangle = delta_gain.triangle - delta_loss.triangle;
		new_move.delta.edge = delta_gain.edge - delta_loss.edge;
		new_move.delta.seqsim = pruned_w[this->n2 * new_move.e[0][0] + new_move.e[0][1]] -
				pruned_w[this->n2 * align->left_match[new_move.m_id[0]] + align->right_match[new_move.m_id[0]]];

		new_move.delta.ortho_count = (pruned_w[this->n2 * new_move.e[0][0] + new_move.e[0][1]] > 0? 1:0) -
				(pruned_w[this->n2 * align->left_match[new_move.m_id[0]] + align->right_match[new_move.m_id[0]]] > 0? 1:0);
	}
	else  {
		// Aggregated removal cost and update alignment with the new move

		for(k = 0; k < new_move.move_no; k++) {
			delta_loss = Delta_removeMatch(align->left_match, align->right_match, new_move.m_id[k]);

			new_move.delta.triangle -= delta_loss.triangle;
			new_move.delta.edge -= delta_loss.edge;
			new_move.delta.seqsim -= pruned_w[this->n2 * align->left_match[new_move.m_id[k]] + align->right_match[new_move.m_id[k]]];
			new_move.delta.ortho_count -= (pruned_w[this->n2 * align->left_match[new_move.m_id[k]] + align->right_match[new_move.m_id[k]]] > 0? 1:0);

			new_mi[new_move.m_id[k]] = new_move.e[k][0];
			new_mj[new_move.m_id[k]] = new_move.e[k][1];
		}
/*		printf("\tStep 1 = %lf\n", new_move.new_move.delta.triangle);*/


		unsigned register int i = new_move.m_id[0], j = new_move.m_id[1];
		// 1) triangles that we removed twice: add them back in!
		if(G->getEdge(align->left_match[i], align->left_match[j]) && H->getEdge(align->right_match[i], align->right_match[j]) ) { // Removed (alignment) nodes were connected, so they could have been part of shared triangles
			new_move.delta.edge ++;
			for(k = 0; k < align->match_no; k++) {
				if(k == i || k == j)
					continue;

				if( (G->getEdge(align->left_match[i], align->left_match[k]) && H->getEdge(align->right_match[i], align->right_match[k])) ) {
					if ((G->getEdge(align->left_match[j], align->left_match[k]) && H->getEdge(align->right_match[j], align->right_match[k])) ) {
						new_move.delta.triangle ++;
					}
				}
			}
		}
/*		printf("\tStep 2 = %lf\n", new_move.new_move.delta.triangle);*/

		// Compute the gain from added edges
		for(k = 0; k < new_move.move_no; k++) {
			delta_gain =  Delta_addMatch(new_mi, new_mj, new_move.m_id[k], new_move.e[k]);

			new_move.delta.triangle += delta_gain.triangle;
			new_move.delta.edge += delta_gain.edge;
			new_move.delta.seqsim += pruned_w[this->n2 * new_move.e[k][0] + new_move.e[k][1]];
			new_move.delta.ortho_count += (pruned_w[this->n2 * new_move.e[k][0] + new_move.e[k][1]] > 0? 1:0);
		}
/*		printf("\tStep 3 = %lf\n", new_move.new_move.delta.triangle);*/


		// 2) triangles that we added twice: remove them!
		if(G->getEdge(new_mi[new_move.m_id[0]], new_mi[new_move.m_id[1]]) && H->getEdge(new_mj[new_move.m_id[0]], new_mj[new_move.m_id[1]])) { // Added (matching) nodes are connected (in the alignment graph), so they can be part of shared triangles
			new_move.delta.edge --;
			for(k = 0; k < align->match_no; k++) {
				if((int)k == new_move.m_id[0] || (int)k == new_move.m_id[1])
					continue;

				if( (G->getEdge(new_mi[new_move.m_id[0]], align->left_match[k]) && H->getEdge(new_mj[new_move.m_id[0]], align->right_match[k])) ) {
					if( (G->getEdge(new_mi[new_move.m_id[1]], align->left_match[k]) && H->getEdge(new_mj[new_move.m_id[1]], align->right_match[k])) ) {
						new_move.delta.triangle --;
					}
				}
			}
		}
/*		printf("\tStep 4 = %lf\n", new_move.new_move.delta.triangle);*/

/*		printf("Actual edges = %ld, triangles = %ld\n", countEdgesUnderAlignment(new_mi, new_mj) - countEdgesUnderAlignment(align->left_match, align->right_match), countTrianglesUnderAlignment(new_mi, new_mj) - countTrianglesUnderAlignment(align->left_match, align->right_match));
		printf("Estimated edges = %ld, triangle = %ld\n", new_move.delta.edge, new_move.delta.triangle);*/
	}



/*	new_move.delta.score = ((double)new_move.delta.edge / align->conserved_edges) + ((double)new_move.delta.triangle / align->conserved_triangles) + (round(new_move.delta.seqsim) / align->total_seqSim);*/
	alpha = 0;

/*	new_move.delta.score = new_move.delta.triangle; //alpha*round( align->match_no*(double)new_move.delta.triangle / align->conserved_triangles ) + (1-alpha)*round( align->match_no*new_move.delta.seqsim / align->total_seqSim );*/




/*	new_move.delta.score = alpha*round( (double)align->match_no*new_move.delta.triangle / align->conserved_triangles ) + (1-alpha)*round( (double)align->match_no*new_move.delta.ortho_count / align->ortho_count );*/
	new_move.delta.score = alpha*round( (double)align->match_no*new_move.delta.triangle / align->conserved_triangles ) + (1-alpha)*round( (double)align->match_no*new_move.delta.seqsim / align->total_seqSim );

	return new_move.delta.score;
}

void ProdTensor::copyMove(Move &dst, Move& src) {
	register unsigned int i;

	memcpy(&(dst.delta), &(src.delta), sizeof(Delta));

	dst.move_no = src.move_no;
	for(i = 0; i < src.move_no; i++) {
		dst.m_id[i] = src.m_id[i];
		dst.e[i][0] = src.e[i][0];
		dst.e[i][1] = src.e[i][1];
	}
}

void ProdTensor::applyMove(Move &best_move, alignment* align) {
	register unsigned int i;

	for(i = 0; i < best_move.move_no; i++) {
		align->left_match[best_move.m_id[i]] = best_move.e[i][0];
		align->right_match[best_move.m_id[i]] = best_move.e[i][1];
	}
}

eigen *ProdTensor::issHOPM(int max_it, double shift_param, double weight_param, double epsilon, double *w, double *x0, int init_type) {
	
	register unsigned int i, j;
	char x_path[1024];
	
	this->best_x = vec(x0, this->n);
	this->x_vec = vec(x0, this->n); //.set_size(this->n, 1);
	this->x_hat_vec.set_size(this->n, 1);
	w_vec = vec(w, this->n);
	
	this->X.set_size(this->n1, this->n2);


	this->res = new eigen;
	this->res->flag = -2;

	// Misuse of notation in the paper. Let's change the order to avoid confusion
	this->alpha = shift_param;

	printf("ISS-HOPM::Max Iter = %d, Alpha (shift) = %e, Beta (weight) = %e, Epsilon = %e\n", max_it, this->alpha, this->beta, epsilon); 	fflush(stdout);


	stats myStats, bestStats; bestStats.T = bestStats.M = -1;
	wall_clock timer;

	
	double old_lambda = -datum::inf, delta_lambda;

	for(i = 0; i < (unsigned int)max_it; i++) {
		timer.tic();
		printf("\t\tIteration %d ... \n", i+1); fflush(stdout);


		/*************************************
		 * Eigen-pair computations starts here
		 *************************************/
		impTTV(x_hat_vec.memptr(), x_vec.memptr()); // x_hat_vec(i) <-- T*x_(i-1)^2
		printf("\t\t\tdt impTTV = %f (norm x_hat = %.2f)\n", timer.toc(), arma::norm(x_hat_vec));

		// Compute Lambda first: x_(i-1) * x_hat_(i)
		timer.tic();
		res->lambda = dot(x_hat_vec, x_vec);
		printf("\t\t\tdt dot = %f\n", timer.toc());

		// Then shift x_hat_vec(i) ...
		timer.tic();
		printf("\t\t\tShifting now:\n");
		//Shift(x_hat_vec.memptr(), x_vec.memptr(), w, alpha, beta, shift_type);
		//Shift(shift_type);
		x_hat_vec += this->alpha*x_vec;			
		printf("\t\t\tdt Shift = %f (norm x_hat_vec = %e)\n", timer.toc(), arma::norm(x_hat_vec));



		// And finally normalize x_hat_vec(i) wrt previous vectors to have ||x_vec(i)|| = 1
		timer.tic();
		x_vec = x_hat_vec / arma::norm(x_hat_vec);
		printf("\t\t\tdt normalization = %f (norm(x_vec) = %.2f)\n", timer.toc(), arma::norm(x_vec));
		fflush(stdout);



		/*
		printf("\t\tMixing x_vec with SeqSim\n");
		timer.tic();

		double topo_median = arma::median(nonzeros(x_vec));
		printf("\t\t\tdt median = %f, median = %e\n", timer.toc(), topo_median);

		timer.tic();
		for(ii = 0; ii < this->n; ii++) {
			if(0 < (x_vec.memptr())[ii])
				topo_absolute_diff[ii] = abs((x_vec.memptr())[ii] - topo_median);
		}
		printf("\t\t\tdt abs diff = %f\n", timer.toc());

		timer.tic();
		double topo_mad = arma::median(nonzeros(topo_absolute_diff));
		printf("\t\t\tdt MAD = %f, MAD = %e\n", timer.toc(), topo_mad);

		timer.tic();
		for(ii = 0; ii < this->n; ii++) {
			if(0 < (x_vec.memptr())[ii]) {
				double topo_vec_z_neg = 0.6745*(topo_median - (x_vec.memptr())[ii]) / topo_mad;
				(topo_sigmoid.memptr())[ii] = 1 / (1 + exp(topo_vec_z_neg));
			}
		}
		printf("\t\t\tdt computing sigmoid %f\n", timer.toc());
		mixed_sim = this->beta*topo_sigmoid + (1-this->beta)*seqsim_sigmoid;
*/


		/*********************************************************
		 * Done with computation, time to evaluate and export now!
		 ********************************************************/
		// Evaluate current results
		timer.tic();
		delta_lambda = abs(res->lambda-old_lambda);
		//myStats = evalX(x_vec.memptr());
		myStats = evalX(x_vec.memptr()); // also sets mi, mj because it runs matching
		printf("\t\t\tmatchin_score = %e, matched edges = %ld, matched tries = %ld, locality = %.2e, dt = %f\n", myStats.matching_score, myStats.M, myStats.T, myStats.locality, timer.toc());
		fflush(stdout);

		//x_vec = (this->beta/myStats.matching_score)*x_vec + (1-this->beta)*w_vec;
		

		/*
		// Run 2nd round of matching after scaling
		timer.tic();		
		mixed_x = (this->beta/myStats.matching_score)*arma::log10(1+100*x_vec) + (1-this->beta)*w_vec;
		myStats = evalX(mixed_x.memptr()); // also sets mi, mj because it runs matching
		printf("\t\t\tmatchin_score (2nd round) = %e, matched edges = %ld, matched tries = %ld, dt = %f\n", myStats.matching_score, myStats.M, myStats.T, timer.toc());
		fflush(stdout);
		*/

		// Export binarized matrix X
		timer.tic();
		sprintf(x_path, "%s/%s_sparsity=%d_alpha=%.2e_it%d_X.smat", output_path, prefix, sparsity_type, alpha, i+1);
		FILE *fd = fopen(x_path, "w");
		fprintf(fd, "%d\t%d\t%d\n", this->n1, this->n2, myStats.num_matches);
		for(j = 0; j < myStats.num_matches; j++) {
			fprintf(fd, "%d\t%d\t1\n", (int)mi[j], (int)mj[j]);
		}
		fclose(fd);
		printf("\t\t\tdt binary export = %f\n", timer.toc());


		if(bestStats.T < myStats.T) {
			bestStats.matching_score = myStats.matching_score;
			bestStats.M = myStats.M;
			bestStats.T = myStats.T;
			this->best_x = x_vec;
		}			
					



		if(delta_lambda < epsilon) {
			break;
		}
		else {
			old_lambda = res->lambda;
		}
	}


	/**********************************************
	 * 			Clean output
	 * *******************************************/	 
	// Export full matrix X
/*
	timer.tic();
	X = (reshape(best_x, n2, n1)).t();
	sprintf(x_path, "%s/%s_alpha=%.2e_beta=%.2e_X.mat", output_path, prefix, alpha, beta);
	X.save(x_path, raw_ascii);
	printf("\t\t\tdt full export = %f\n", timer.toc());
*/

	alignment* result = postprocess(best_x.memptr(), 10, 200, 10, 5);
	printf("After post processing:: Triangles = %ld, edges = %ld\n ", result->conserved_triangles, result->conserved_edges);
	

	timer.tic();
	sprintf(x_path, "%s/%s_sparsity=%d_alpha=%.2e_X.smat", output_path, prefix, sparsity_type, alpha);
	printf("Exporting final solution to %s ... ", x_path);
	FILE *fd = fopen(x_path, "w");
	fprintf(fd, "%d\t%d\t%d\n", this->n1, this->n2, (int)result->match_no);
	for(j = 0; j < (unsigned int)result->match_no; j++) {
		fprintf(fd, "%d\t%d\t1\n", result->left_match[j], result->right_match[j]);
	}
	fclose(fd);
	printf("Done!\n");


/*	// Export binarized matrix X
	myStats = evalX(best_x.memptr()); // also sets mi, mj because it runs matching
	printf("\t\t\tFINAL matchin_score = %e, matched edges = %ld, matched tries = %ld, dt = %f\n", myStats.matching_score, myStats.M, myStats.T, timer.toc());
	fflush(stdout);

	timer.tic();
	sprintf(x_path, "%s/%s_alpha=%.2e_beta=%.2e_X.smat", output_path, prefix,  alpha, beta);
	FILE *fd = fopen(x_path, "w");
	fprintf(fd, "%d\t%d\t%d\n", this->n1, this->n2, myStats.num_matches);
	for(j = 0; j < myStats.num_matches; j++) {
		fprintf(fd, "%d\t%d\t1\n", (int)mi[j], (int)mj[j]);
	}
	fclose(fd);
*/


/*	vector<int> MI, MJ;
	for(i = 0; i < myStats.num_matches; i++) {
		MI.push_back(mi[i]);
		MJ.push_back(mj[i]);
	}
	printf("Tris = %ld\n", countTrianglesUnderAlignment(MI, MJ));

	long Tsum = 0;
	for(i = 0; i < myStats.num_matches; i++) {
		long curr_sum = Delta_removeMatch(MI, MJ, i);
		printf("%d- %ld\n", i, curr_sum);
		Tsum += curr_sum;
	}
	printf("Total sum = %ld\n", Tsum);
*/

	res->flag = 0;
	return res;
}


Delta ProdTensor::Delta_removeMatch(vector<int> mi, vector<int> mj, unsigned int i) {
	register unsigned int j, k;

	Delta delta;
	delta.triangle = 0;
	delta.edge = 0;
	delta.TriSim =  w_vec[this->n2*mi[i] + mj[i]];

	for(j = 0; j < mi.size(); j++) {
		if(i == j)
			continue;

		if( !G->getEdge(mi[i], mi[j]) || !H->getEdge(mj[i], mj[j]) )
			continue;
		else {
			delta.edge++;
			delta.TriSim += w_vec[this->n2*mi[j] + mj[j]];
		}

		for(k = j+1; k < mi.size(); k++) {
			if(i == k)
				continue;
			if( (G->getEdge(mi[i], mi[k]) && H->getEdge(mj[i], mj[k])) ) {
				if ( (G->getEdge(mi[j], mi[k]) && H->getEdge(mj[j], mj[k])) ) {
					delta.triangle ++;
				}
			}
		}
	}

	return delta;
}


Delta ProdTensor::Delta_addMatch(vector<int> mi, vector<int> mj, unsigned int i, vector<int> e) {
	register unsigned int j, k;


	Delta delta;
	delta.triangle = 0;
	delta.edge = 0;
	delta.TriSim =  w_vec[this->n2*e[0] + e[1]];

	for(j = 0; j < mi.size(); j++) {
		if(i == j)
			continue;

		if( !G->getEdge(e[0], mi[j]) || !H->getEdge(e[1], mj[j]) )
			continue;
		else {
			delta.edge++;
			delta.TriSim += w_vec[this->n2*mi[j] + mj[j]];
		}

		for(k = j+1; k < mi.size(); k++) {
			if(i == k)
				continue;
			if( (G->getEdge(e[0], mi[k]) && H->getEdge(e[1], mj[k])) ) {
				if ( (G->getEdge(mi[j], mi[k]) && H->getEdge(mj[j], mj[k])) ) {
					delta.triangle ++;
				}
			}
		}
	}

	return delta;
}


long ProdTensor::countTrianglesUnderAlignment(vector<int> mi, vector<int> mj) {
	register unsigned int i, j, k;
	 
	 long tri_count = 0;

	 // First count the triangles that fall exclusively with the current block 
     for(i = 0; i < mi.size(); i++) {
		for(j = i+1; j < mi.size(); j++) {
			if(!G->getEdge(mi[i], mi[j]) || !H->getEdge(mj[i], mj[j]) )
				continue;				
			for(k = j+1; k < mi.size(); k++) {				
				if( (G->getEdge(mi[i], mi[k]) && H->getEdge(mj[i], mj[k])) &&
					(G->getEdge(mi[j], mi[k]) && H->getEdge(mj[j], mj[k])) ) {
					tri_count ++;
				}
			}
		}
	}	
	
	return tri_count;
}

long ProdTensor::countEdgesUnderAlignment(vector<int> mi, vector<int> mj) {
	register unsigned int i, j;

	 long edge_count = 0;

     for(i = 0; i < mi.size(); i++) {
		for(j = i+1; j < mi.size(); j++) {
			if(G->getEdge(mi[i], mi[j]) && H->getEdge(mj[i], mj[j]) )
				edge_count++;
		}
	}

	return edge_count;
}

int ProdTensor::InitX(double *x, int init_type, double *w) {
	register unsigned int i, j, cur_row, cur_col;
	printf("InitX:: Initializing x0 with type %d\n", init_type);
	
	bzero(x, n*sizeof(double));	
	printf("\tInitX::Sparsity flag = %d\n", this->sparsity_type);
	switch(init_type) {
		case(Uniform):	{
			double uniform_val = 1 / sqrt(n);
			for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
				i = unmatched_rows[cur_row];
				for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {		
					j = unmatched_cols[cur_col];
					x[i * n2 + j] = uniform_val;
				}
			}
			break;		
		}	
		case(Random): {  // Warning:: Using random initialization with sparse formulation is highly discouraged.
			fprintf(stderr, "\tInitX: *Warning* Using random initialization with sparse formulation is highly discouraged.\n");
			srand(time(NULL));
			for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
				i = unmatched_rows[cur_row];
				for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {		
					j = unmatched_cols[cur_col];
					x[i * n2 + j] = (rand() % n);
				}
			}		
			break;		
		}	
		case(SeqSim): {
			for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
				i = unmatched_rows[cur_row];
				for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {		
					j = unmatched_cols[cur_col];
					x[i * n2 + j] = w[i*n2 + j];
				}
			}		
			break;		
		}
		default:
			fprintf(stderr, "InitX:: Unknown initialization type: %d\n", init_type);
			return -1;
	}
	normalize(x, x, n);

		
	return 0;
}

// Initializes the bipartite graph between vertices in G and H over the set of unmatched rows/cols
unsigned int ProdTensor::setupBipartiteGraph(double *x_in) {
	register unsigned int i, j, cur_row, cur_col;
	unsigned int num_edges = 0;
	if( this->sparsity_type == NoSparsity ) {
		num_edges = unmatched_rows.size() * unmatched_cols.size();

		// <i, j>: original vertex<row, col> index in G and H, <cur_row, cur_col>: Index in the unmatched submatrix
		for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
			i = unmatched_rows[cur_row];
			for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {
				j = unmatched_cols[cur_col];
				rows[cur_row*unmatched_cols.size()+ cur_col] = i;
				cols[cur_row*unmatched_cols.size()+ cur_col] = j;
				edge_weights[cur_row*unmatched_cols.size()+ cur_col] = x_in[i*n2+j];
			}
		}
	}
	else {
		register unsigned int i;
		num_edges = unmatched_pairs_row.size();
		for (i = 0; i < num_edges; i++) {
			rows[i] = unmatched_pairs_row[i];
			cols[i] = unmatched_pairs_col[i];
			edge_weights[i] = x_in[unmatched_pairs_row[i]*n2 + unmatched_pairs_col[i]];
		}
	}

	return num_edges;
}

ProdTensor::~ProdTensor() {   
	
	free(G_isMatched);
	free(H_isMatched);		
	G_isMatched = H_isMatched = NULL;
				
	free(rows);
	free(cols);
	free(edge_weights);
	free(mi);
	free(mj);
	rows = cols = edge_weights = mi = mj = NULL;
}
					
