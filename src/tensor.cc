#include <tensor.h>
//#define VERBOSE_DEBUG

extern char *FileType_names[];
extern char *ShiftType_names[];
extern char *InitType_names[];

char cmd[2048];			

bool *G_isMatched, *H_isMatched;
double *rows, *cols, *edge_weights, *mi, *mj; // matching aux. vectors



extern FILE* log_fd;

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
 int counter = 0;
 
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


double ProdTensor::unmatched_norm(double *v, unsigned long n) { 

	register double sqsum = 0;
	register unsigned int cur_row, cur_col;
	unsigned register long int i, j;  

	for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
		i = unmatched_rows[cur_row];
		for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {		
			j = unmatched_cols[cur_col];
			sqsum += (v[i*n2+j]*v[i*n2+j]);
		}
	}
    
    return sqrt(sqsum);
}


// Normalizes v to have norm = weight over the unmatched subset vertices
double ProdTensor::unmatched_normalize(double *v_norm, double *v, unsigned long n, double weight) {
    
    double Norm =  unmatched_norm(v, n);
    if(Norm < 1e-52) 
		return 0;		
    double scale = weight / Norm;    



	register unsigned int cur_row, cur_col;
	unsigned register long int i, j;  
	for(cur_row = 0; cur_row < unmatched_rows.size(); cur_row++) {
		i = unmatched_rows[cur_row];
		for(cur_col = 0; cur_col < unmatched_cols.size(); cur_col++) {		
			j = unmatched_cols[cur_col];
			v_norm[i*n2+j] = scale*v[i*n2+j];
		}
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


void ProdTensor::scaled_impTTV(double *new_x, double *x) {
	unsigned int *Ni1, *Ni2, *p, *q;
	register unsigned int i1, j1, k1, i2, j2, k2, idx, cur_row, cur_col;
	
	counter_type **D_G = T_G->get_triEdgeDegs(); // Sym: D_G(i, j) = D_G(j, i)
	counter_type **D_H = T_H->get_triEdgeDegs(); // Sym: D_H(i', j') = D_G(j', i')
	

	// Ensure fixed submatrices are untouched.
	memcpy(new_x, x, n*sizeof(double));
	
	double scale_factor;
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
					scale_factor = (D_G[j1][k1]*D_H[j2][k2]);
					if(scale_factor == 0) {
						printf("Fatal error in scaled_impTTV: scale factor is 0 ( <%d, %d, %d> vas <%d, %d, %d> )\n", i1, j1, k1, i2, j2, k2);
						exit(255);
					}
					new_x[idx] += ( (x[n2*j1 + j2] * x[n2*k1 + k2]) + (x[n2*j1 + k2] * x[n2*k1 + j2]) ) / scale_factor; 
				}
			}
			new_x[idx] *= 2; // To account for the symmetry in T, i.e., X[j1][j2]*X[k1][k2] vs X[k1][k2]*X[j1][j2], (m-1)! redundancy in permutations.
		}		
	}		
}

eigen* ProdTensor::MLPR(double alpha, double *w, double *x0, int max_it, double epsilon, bool verbose) {
	register unsigned int i, j, k, l;
	char x_path[1024];
	FILE *x_fd;
	
	printf("\t\tMulti-Linear PageRank:: Alpha = %e, Max Iter = %d, Epsilon = %e\n", alpha, max_it, epsilon); fflush(stdout);
		
	
	double delta = 0, matching_score;//, residual_weight;	  
    int noutedges;
    
	// x0 and w has been previously norm-2 normalized. Re-normalize to ensure norm(res->x, 1) == 1
	memcpy(res->x, x0, n*sizeof(double));  	
	double sum_counter = 0, sum_counter2 = 0;
	for(i = 0; i < n; i++) {
		sum_counter += res->x[i];
		sum_counter2 += w[i];
	}
	for(i = 0; i < n; i++) {
		res->x[i] /= sum_counter;
		w[i] /= sum_counter2;
	}	

	double *y = new (nothrow) double [n];
	if(y == 0) {
		fprintf(stderr, "issHOPM:: Can't allocate memory for 'y'\n");
		exit(-1);
	}
		
	res->flag = -2;
	
	
	for(i = 0; i < (unsigned int)max_it; i++) {	
		clock_t startTime = clock();		
					
		res->its = i+1;
		scaled_impTTV(y, res->x); // y <-- T_tilde*res->x ^2
//		printVec(stdout, y, n1, n2);
	
		sum_counter = 0;
		for(j = 0; j < n; j++) {
			sum_counter += y[j];
		}	

		for(j = 0; j < n; j++) {
			y[j] = alpha * y[j] / sum_counter + (1-alpha) * w[j];
		}		
		
		delta = 0; 
		double newx_norm = 0;
		for(j = 0; j < n; j++) {
			newx_norm += y[j];
			delta += fabs(res->x[j] - y[j]);
		}				
		memcpy(res->x, y, n*sizeof(double));		
		printf("delta = %e, ||x||_1 = %e\n", delta, newx_norm);

		 
/*		sum_counter = 0;
		for(j = 0; j < n; j++) {
			sum_counter += y[j];
			y[j] *= alpha;
		}	

		residual_weight = 1 - (alpha * sum_counter);
//		printf("alpha = %e, sum = %e, Residual weight = %e\n", alpha, sum_counter, residual_weight);
//		printVec(stdout, y, n1, n2);

		for(j = 0; j < n; j++) {
			y[j] += residual_weight * w[j];
		}		
		
		delta = 0; 
		double newx_norm = 0;
		for(j = 0; j < n; j++) {
			newx_norm += y[j];
			delta += fabs(res->x[j] - y[j]);
		}				
		memcpy(res->x, y, n*sizeof(double));		
//		printf("delta = %e, ||x||_1 = %e\n", delta, newx_norm);
//		printVec(stdout, y, n1, n2);
*/


		// Export X per iteration
		sprintf(x_path, "%s/%s_MLPR_%.2e_round%03d_x%d.mat", output_path, prefix, alpha, Round, i+1);
		x_fd = fopen(x_path, "w");
		printVec(x_fd, res->x, n1, n2);
		fclose(x_fd);
		sprintf(cmd, "gzip %s -f", x_path);
		if(system(cmd) == -1) {
			fprintf(stderr, "issHOPM:: Error compressing %s", x_path);
		}



		/*******************************************************
		 * 
		 * Keep user entertained by showing stats per iteration
		 * 
		*******************************************************/		 
		// Compute locality of current solution
		double v1 = 0, v2 = 0, v4 = 0, temp, locality;
		for(j = 0; j < (unsigned int)n; j++) {
			temp = res->x[j];
			v1 += temp;
			temp *= temp;
			v2 += temp;
			temp *= temp;
			v4 += temp;
		}

		locality = round(v2*v2/v4);
		//locality = round(v1*v1/v2);

		if(verbose) {			
			double num_edges = setupBipartiteGraph(res->x);
			matching_score = match(n1, n2, num_edges, rows, cols, edge_weights, mi, mj, &noutedges);	

			sprintf(x_path, "%s/%s_MLPR_%.2e_round%03d_x%d.smat", output_path, prefix, alpha, Round, i+1);
			x_fd = fopen(x_path, "w");
			fprintf(x_fd, "%d\t%d\t%d\n", n1, n2, noutedges);
			for(k = 0; k < (unsigned)noutedges; k++) {
				fprintf(x_fd, "%d\t%d\t%d\n", (int)(mi[k]), (int)(mj[k]), 1);
			}		
			fclose(x_fd);
			
			
			bool **alignGraph = (bool **) calloc(noutedges, sizeof(bool *));
			if(alignGraph == NULL) {
				fprintf(stderr, "\t\tissHOPM:: Can't allocate memory for alignGraph.\n");
				exit(-1);
			}
			for(j = 0; j < (unsigned)noutedges; j++) {
				alignGraph[j] = (bool *) calloc(noutedges, sizeof(bool));
				if(alignGraph[j] == NULL) {
					fprintf(stderr, "\t\tissHOPM:: Can't allocate memory for alignGraph[%d].\n", i);
					exit(-1);
				}
			}

			// Construct alignment graph
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
			

			long T = 0, M = 0;
			for(j = 0; j < (unsigned)noutedges; j++) {
				for(k = j+1; k < (unsigned)noutedges; k++) {
					if( !alignGraph[j][k] )
						continue;	
					M++;				
					for(l = k+1; l < (unsigned)noutedges; l++) {
						if( alignGraph[j][l] && alignGraph[k][l] ) {
								T++;
						}										
					}
				}
			}							

	
			clock_t endTime = clock();
			double time_per_iteration = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;
						 
			printf("\t\tIteration %d\tDelta =%e\tMatching score = %f\tPR = %d \tT = %ld\tM = %ld\tdt = %lf\n", i+1, delta, matching_score, (int)locality, T, M, time_per_iteration);
			fprintf(log_fd, "\t\tIteration %d\tDelta =%e\tMatching score = %f\tPR = %d \tT = %ld\tM = %ld\tdt = %lf\n", i+1, delta, matching_score, (int)locality, T, M, time_per_iteration);
			fflush(log_fd);
			fflush(stdout);
		}
		else {
			clock_t endTime = clock();
			double time_per_iteration = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;			
			printf("\t\tIteration %d\tDelta =%e\tdt = %lf s\n", i+1, delta, time_per_iteration);		
			fprintf(log_fd, "Iteration %d\tDelta =%e\tdt = %lf s\n", i+1, delta, time_per_iteration);		
			fflush(log_fd);
			fflush(stdout);			
		}
		
		if(fabs(delta) < epsilon) {
			break;
		}		
		res->trace[i] = delta;
	}
		
	delete [] y;

	res->flag = 0;
	return res;
	
}


vec ProdTensor::Combine_basis(mat B) {
	stats myStats;
	vec v_out = vec(size(B, 1));	
	mat W = zeros(size(B, 1));

	double *temp_x = new  (nothrow) double[n];
	if(temp_x == 0) {
		fprintf(stderr, "Combine_basis:: Can't allocate memory for 'temp_x'\n");
		exit(-1);
	}
		
	for(unsigned int j = 0; j < this->subspace_dim; j++) {
		v_out = B.col(j);
		
		memcpy(temp_x, v_out.memptr(), n*sizeof(double) );					
		myStats = evalX(temp_x);
		W(j, j) = myStats.T;			
	}

	W = W / sum(W.diag());	
	printf("Sum of W = %f\n", (float)sum(sum(W)));	
	W.print();	
	v_out = sum(B*W, 1);	
	
	
	delete [] temp_x;
	return v_out;
}


// x is row-order linearized version of matrix X
stats ProdTensor::evalX(double *x) {
	stats myStats;
	int noutedges;

	unsigned register int j, k, l;


	/*double v1 = 0, v2 = 0, v4 = 0, temp;
	for(j = 0; j < (unsigned int)n; j++) {
		temp = x[j];
		v1 += temp;
		temp *= temp;
		v2 += temp;
		temp *= temp;
		v4 += temp;
	}
	myStats.locality = round(v2*v2/v4);
	*/
	
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
void ProdTensor::ScaleX(double *x, int n) {
	double EPSILON = numeric_limits<double>::epsilon();
	double MIN_DOUBLE = numeric_limits<double>::min();
	double MAX_DOUBLE = numeric_limits<double>::max();
	printf("Min = %e, Max = %e, Epsilon = %e\n", MIN_DOUBLE, MAX_DOUBLE, EPSILON);

	double Min = MAX_DOUBLE, Max = MIN_DOUBLE;
	register int i;

	for(i = 0; i < n; i++) {
		if(x[i] > EPSILON) {
			x[i] = log(x[i]);
			if(x[i] < Min)
				Min = x[i];
			else if(x[i] > Max)
				Max = x[i];
		}
		else { // get rid = negatives, if exist
			x[i] = 0;
		}
	}
	printf("Min of vec = %e, Max of vec = %e\n", Min , Max);

	double Range = Max - Min;
	for(i = 0; i < n; i++) {
		if(abs(x[i]) > EPSILON) { // if it's non-zero
			x[i] = (x[i] - Min) / Range;
		}
	}
}


void ProdTensor::Orthogonalize(mat M) {
	qr_econ(Q, R, M);
	S = diagmat(sign(R.diag()));
	B = Q*S;
}

vec ProdTensor::combineB(mat B) {
	vec x_vec(this->n), z = zeros(this->n);
	/*
			// TODO: How to optimally integrate?
			for(j = 0; j < this->subspace_dim; j++) {
				x_vec = B.col(j);

				memcpy(res->x, x_vec.memptr(), n*sizeof(double) );
				myStats = evalX(res->x);
				S(j, j) = myStats.T;
			}
			S = S / sum(S.diag());
			x_vec = sum(B*S, 1);
		*/

	x_vec = arma::max(arma::max(B, 1), z);
	x_vec = x_vec / sum(x_vec);

	return x_vec;
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

void ProdTensor::Project(int k) {

	/*
	 double r;
	 // Orthogonalization
	 for(register int j = 0; j < k; j++) {
	 	 r = dot(B.col(j), x_vec);
	 	 x_vec = x_vec - (r*B.col(j));
	 }

	 // Normalization
	 x_vec = x_vec / arma::norm(x_vec);
	 */
	qr_econ(Q, R, B.cols(0, k));
	S = diagmat(sign(R.diag()));
	B.cols(0, k) = Q*S;
}

void ProdTensor::CombineB(mat B, double *w, int k) {
	vec z = zeros(this->n), current_col(this->n);

/*

	x_vec = arma::max(arma::max(B.cols(0, k), 1)+vec(w, this->n), z);
	x_vec = x_vec / arma::norm(x_vec);
*/
	x_vec = z;
	for(register int i = 0; i <=k; i++) {
		current_col = arma::max(B.col(i), z); // remove negative elements
		current_col = current_col / arma::norm(current_col); // re-normalize
		x_vec = x_vec + current_col;
	}
	x_vec = x_vec / arma::norm(x_vec);

	x_vec = x_vec + vec(w, this->n);
	x_vec = x_vec / arma::norm(x_vec);
}

eigen *ProdTensor::issHOPM(int subspace_dim, int max_it, double alpha, double beta, double epsilon, double *w, double *x0, int init_type, int shift_type) {
	
	register unsigned int i, j, k;
	char x_path[1024], log_str[1024];

	vec tri_count = zeros(max_it, 1);

	this->subspace_dim = subspace_dim;

	this->x_vec = vec(x0, this->n); //.set_size(this->n, 1);
	this->x_hat_vec.set_size(this->n, 1);
	this->w_vec.set_size(this->n, 1);

	this->X.set_size(this->n1, this->n2);

	this->Q.set_size(n, this->subspace_dim);
	this->R.set_size(this->subspace_dim, this->subspace_dim);
	this->S.set_size(this->subspace_dim, this->subspace_dim);
	this->M.set_size(this->n, this->subspace_dim);
	this->B.set_size(this->n, this->subspace_dim);

	this->res = new eigen;
	this->res->flag = -2;

	this->alpha = alpha;
	this->beta = beta;
	
	if(beta != 0) {
		this->sigma = (1-beta)/beta;
	}
	else {
		printf("beta = 0, there is nothing to be done!\n");
		return NULL;
	}

	printf("ISS-HOPM:: Dim = %d, Max Iter = %d, Alpha = %e, Beta = %e, Epsilon = %e, Shift_type = %s\n", this->subspace_dim, max_it, alpha, beta, epsilon, ShiftType_names[shift_type]); 	fflush(stdout);


	stats myStats;
	wall_clock timer;
	/*
	if(w != NULL) {
		if(normalize(w, w, n)  == 0)  {
		  fprintf(stderr, "tame:: w vector has norm 0.\n");
		}
	}
	*/


	/*int noutedges = 0;
	setupBipartiteGraph(w);
	myStats.matching_score = match(n1, n2, n1*n2, rows, cols, edge_weights, mi, mj, &noutedges);
	for (i = 0; i < (unsigned int)noutedges; i++) {
		this->w_vec(mi[i]*n2+mj[i]) = edge_weights[i];
	}*/
	int noutedges;
	setupBipartiteGraph(w);
	double matching_score = match(n1, n2, unmatched_rows.size()*unmatched_cols.size(), rows, cols, edge_weights, mi, mj, &noutedges);	
	this->w_vec = vec(w, this->n) / matching_score; // % vec(w, this->n);
	printf("\tNorm-1 of w matching = %e\n", matching_score);	
	

/*	// Estimate relative weight of seqsim/tri match (rho)
	if(w != NULL) {
		if(normalize(w, w, n)  == 0)  {
		  fprintf(stderr, "tame:: w vector has norm 0.\n");
		}

		// Computing relative weight of triangular matching and node similarity scores
		double mean_seqsim = 0, mean_uniTris = 0, epsilon = 1.0 / n;
		unsigned int *Ni1, *Ni2;

		unsigned int G_trie_deg, H_trie_deg;
		for(i = 0; i < n1; i++) {
			Ni1 = T_G->getAdjTriangles(i);
			G_trie_deg = Ni1[0];
			for(j = 0; j < n2; j++) {
				Ni2 = T_H->getAdjTriangles(j);
				H_trie_deg = Ni2[0];
				mean_uniTris += 12*epsilon*G_trie_deg*H_trie_deg; // Each triangle contributes 3epsilon, there are 2 permutations (jj', kk' vs jk', kj'), and there is a factor of 2 for symmetry (jj', kk' vs kk', jj')
				mean_seqsim += w[n2*i+j];
			}
		}
		this->rho = mean_uniTris / mean_seqsim; // Weight of SeqSim
		printf("\t\tSeqSim weight = %e\tTri weight = %e\tSeqSim Normalization factor = %e\n", mean_seqsim / n, mean_uniTris / n, rho);
	}

	else {
		this->rho = 0;
	}*/
	
	
	double old_lambda = -datum::inf, delta_lambda;
	long old_tries = 0, delta_tries;
	//vec x_subset_vec; x_subset_vec.set_size(this->n);


	/*register unsigned int ii;
	vec topo_sigmoid = zeros(this->n, 1), topo_absolute_diff = zeros(this->n, 1), seqsim_sigmoid = zeros(this->n, 1), seqsim_absolute_diff = zeros(this->n, 1), mixed_sim;
	printf("\t\tNormalizing w_vec\n");
	timer.tic();
	double seqsim_median = arma::median(nonzeros(w_vec));
	printf("\t\t\tdt median = %f, median = %e\n", timer.toc(), seqsim_median);

	timer.tic();
	for(ii = 0; ii < this->n; ii++) {
		if(0 < (w_vec.memptr())[ii])
			seqsim_absolute_diff[ii] = abs((w_vec.memptr())[ii] - seqsim_median);
	}
	printf("\t\t\tdt abs diff = %f\n", timer.toc());

	timer.tic();
	double seqsim_mad = arma::median(nonzeros(seqsim_absolute_diff));
	printf("\t\t\tdt MAD = %f, MAD = %e\n", timer.toc(), seqsim_mad);

	timer.tic();
	for(ii = 0; ii < this->n; ii++) {
		if(0 < (w_vec.memptr())[ii]) {
			double seqsim_vec_z_neg = 0.6745*(seqsim_median - (w_vec.memptr())[ii]) / seqsim_mad;
			(seqsim_sigmoid.memptr())[ii] = 1 / (1 + exp(seqsim_vec_z_neg));
		}
	}
	printf("\t\t\tdt computing sigmoid %f\n", timer.toc());
*/


	for(k = 0; k < this->subspace_dim; k++) {
		timer.tic();
		printf("\tComponent %d ...\n", k+1); fflush(stdout);

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
			printf("\t\t\t\tNorm of x_hat(k) = %.2f\t Norm of x(k-1) = %.2f\n", arma::norm(x_hat_vec), arma::norm(x_vec));
			//Shift(x_hat_vec.memptr(), x_vec.memptr(), w, alpha, beta, shift_type);
			//Shift(shift_type);
			x_hat_vec += alpha*x_vec;			
			printf("\t\t\tdt Shift = %f (norm x_hat_vec = %e)\n", timer.toc(), arma::norm(x_hat_vec));



			// And finally ortho-normalize x_hat_vec(i) wrt previous vectors to have ||x_vec(i)|| = 1
			timer.tic();
			/*
			B.col(k) = x_hat_vec;
			Project(k); // Projects x_hat_vec(i) to the orthogonal subspace of the span of previous vector
			x_vec = B.col(k);
			*/
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
			printf("\t\t\tdt matching (1st round) matchin_score = %e, dt = %f\n", myStats.matching_score, timer.toc());
			fflush(stdout);

			// Run 2nd round of matching after scaling
			timer.tic();
			x_vec = (this->beta/myStats.matching_score)*x_vec + (1-this->beta)*w_vec;
			myStats = evalX(x_vec.memptr()); // also sets mi, mj because it runs matching
			printf("\t\t\tdt matching (2nd round) matchin_score = %e, dt = %f\n", myStats.matching_score, timer.toc());
			fflush(stdout);


			// Export binarized matrix X
			timer.tic();
			sprintf(x_path, "%s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_it%d_X.smat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1, i+1);
			FILE *fd = fopen(x_path, "w");
			fprintf(fd, "%d\t%d\t%d\n", this->n1, this->n2, myStats.num_matches);
			for(j = 0; j < myStats.num_matches; j++) {
				fprintf(fd, "%d\t%d\t1\n", (int)mi[j], (int)mj[j]);
			}
			fclose(fd);
			printf("\t\t\tdt binary export = %f\n", timer.toc());


			
			/*
			timer.tic();
			X = (reshape(x_vec, n2, n1)).t();
			sprintf(x_path, "%s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_it%d_X.mat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1, i+1);
			X.save(x_path, raw_ascii);
			printf("\t\t\tdt full export = %f\n", timer.toc());
			*/


			// Report stats
			sprintf(log_str, "\t\tEnd of iteration %d\tLambda = %.6e\tDelta Lambda=%.2e\tMatching score = %.2f\tT = %ld\tM = %ld\t%% correct = %.2f\n", i+1, res->lambda, delta_lambda, myStats.matching_score, myStats.T, myStats.M, myStats.correctness);
			printf("%s", log_str); fflush(stdout);
			fprintf(log_fd, "%s", log_str); fflush(log_fd);

			tri_count(i) = myStats.T;


			if(delta_lambda < epsilon) {
				break;
			}
			else {
				old_lambda = res->lambda;
			}
		}

		// Export full matrix (X) after convergence ...
/*
		timer.tic();
		X = (reshape(x_vec, n2, n1)).t();
		sprintf(x_path, "%s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_it%d_X.mat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1, i);
		X.save(x_path, raw_ascii);
		printf("\t\tdt full export = %f\n", timer.toc());

		// Now combine all x vectors and export the combined matrix
		CombineB(B, w, k); // x_vec = max(max(B(:, 1:k), 1), 0)

		myStats = evalX(x_vec.memptr());
		delta_tries = myStats.T - old_tries;
		sprintf(log_str, "\tEnd of Component %d\tMatching score = %.2f\tT = %ld\tM = %ld\n", k+1, myStats.matching_score, myStats.T, myStats.M);
		printf("%s", log_str); fflush(stdout);
		fprintf(log_fd, "%s", log_str); fflush(log_fd);

		// Export binarized version of the combined matrix
		timer.tic();
		sprintf(x_path, "%s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_X.smat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1);
		FILE *fd = fopen(x_path, "w");
		fprintf(fd, "%d\t%d\t%d\n", this->n1, this->n2, myStats.num_matches);
		for(j = 0; j < myStats.num_matches; j++) {
			fprintf(fd, "%d\t%d\t1\n", (int)mi[j], (int)mj[j]);
		}
		fclose(fd);


		// Export full version of the combined matrix

		X = (reshape(x_vec, n2, n1)).t();
		sprintf(x_path, "%s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_X.mat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1);
		X.save(x_path, raw_ascii);
		printf("\t\t\tdt export = %f\n", timer.toc());/
		*/
		int max_it, max_val = 0;
		for(j = 0; j < i; j++) {
			if(max_val < tri_count(j) ) {
				max_val = tri_count(j);
				max_it = j;
			}
		}
		printf("\t\tMax # of triangles obtained is %d @ it%d\n", max_val, max_it); fflush(stdout);
		delta_tries = max_val - old_tries;

		sprintf(x_path, "%s/max.txt", output_path);
		printf("\tWriting max value to %s\n", x_path);
		FILE *max_file = fopen(x_path, "w");
		fprintf(max_file, "%d", max_val);
		fclose(max_file);

		char cpcmd[1024];
		sprintf(cpcmd, "cp %s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_it%d_X.smat %s/%s%s_alpha=%.2e_beta=%.2e_comp%03d_X.smat", output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1, max_val, output_path, prefix, ShiftType_names[shift_type], alpha, beta, k+1);
		int ret = system(cpcmd);
		if(ret) {
			fprintf(stderr, "issHOPM:: can't copy best it X to final X\n");
			exit(255);
		}

		if (delta_tries < 1) { // negligible triangle increase per component
			break;
		}
		else {
			old_tries = myStats.T;
		}
	}


	res->flag = 0;
	return res;
}


double ProdTensor::opt_tri_stat(double *x, double *mi, double *mj, vector<int> match_perm) {
	printf("\t\t\tRunning opt_tri_stat ....\n");
	long double tri_count = 0, penalty = 1;
	register unsigned int i, j, k;
	register unsigned int i1, i2, j1, j2, k1, k2, incident_vertices = 0;
	bool isolated = true;

	double tri_stat, top_tri_stat = 0; // Best triangle aggregation statistic found so far
	int top_tri_index = -1; // Index of cut in the sorted list which optimizes the triangle aggregation statistic
	for(i = 0; i < match_perm.size(); i++) { // vertex being added ...
		i1 = mi[match_perm[i]];
		i2 = mj[match_perm[i]];

		isolated = true;
		for(j = 0; j < i; j++) {
			j1 = mi[match_perm[j]];
			j2 = mj[match_perm[j]];			
			if(!G->getEdge(i1, j1) || !H->getEdge(i2, j2) )
				continue;				
			for(k = 0; k < j; k++) {
				k1 = mi[match_perm[k]];
				k2 = mj[match_perm[k]];			
				if( (G->getEdge(i1, k1) && H->getEdge(i2, k2)) &&
					(G->getEdge(j1, k1) && H->getEdge(j2, k2)) ) {
					tri_count ++;
					isolated = false;
				}				
			}
		}
		if(!isolated)
			incident_vertices++;
			
		penalty = pow(incident_vertices, 1.5);
		tri_stat = tri_count / penalty;
		if(top_tri_stat < tri_stat) {
			top_tri_stat = tri_stat;
			top_tri_index = i;
		}
		//printf("\t\t\t%d- val = %e, tri count = %Le, penalty = %Le, current stat = %e, Top stat = %e, Top index = %d\n", i, val, tri_count, penalty, tri_stat, top_tri_stat, top_tri_index);
	}	
	
	if(top_tri_index == -1) {
		printf("opt_tri_stat:: Block with zero enriched triangles!\n");
		return 0;
	}
	else {
		return top_tri_index+1;
	}
}

double ProdTensor::nonzeros_estimate(vector<double> v) {
	printf("\t\t\tRunning  nonzeros_estimate (bin_type = %d) ....\n", bin_type);
	register int i;
	unsigned int nnz = 0;
	
	int len = v.size();
	double v1 = 0, v2 = 0, v4 = 0, temp;
	
	for(i = 0; i < len; i++) {
		temp = abs(v[i]);
		v1 += temp;
		temp *= temp;
		v2 += temp;
		temp *= temp;
		v4 += temp;
	}
	
	switch(bin_type) {
		case(Norm12_bin): 	// sum(v_i^2)^2 / sum (v_i^4)
			nnz = round(v2*v2/v4);
			break;
		case(Norm24_bin): // sum(v_i)^2 / sum (v_i^2)
			nnz = round(v1*v1/v2);
			break;
			
		case(Z_bin): // positive Z-score -- above average
			v1 /= len; // compute mean;
			for(i = 0; i < len; i++) {
				if(v1 <= v[i])
					nnz++;
			}
			break;
			
		default:
			fprintf(stderr, "nonzeros_estimate: Binarization code %d unknown.\n", bin_type);						
	}
	
	return nnz;
}


typedef pair<double, int> Pair;
typedef vector<Pair>::const_iterator I;

struct CmpPair
{
    bool operator()(const Pair& a, const Pair& b)
    { return b.first < a.first; }
};

void sortingPermutation(const vector<double>& values, vector<double>& sorted_values, vector<int>& permutation) {
    vector<Pair> pairs;
    for (register int i = 0; i < (int)values.size(); i++)
        pairs.push_back(Pair(values[i], i));

    sort(pairs.begin(), pairs.end(), CmpPair());

    for (I p = pairs.begin(); p != pairs.end(); ++p) {		
		sorted_values.push_back(p->first);
        permutation.push_back(p->second);
	}
        
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
	
	// Then count triangles with two end points in the block, and one among already matched pairs in previous blocks


	// Finally count triangles with only one end points in the current block
	
	return tri_count;
}


 // Computes the weight of each Block. rel_weight is the overl weight given to the previous blocks. If -1, the weight of the remaining block is set based on the average of previous blocks
 block_stats ProdTensor::computeWeights(double rel_weight) {

	register unsigned int i, j, k;
	block_stats res;
	res.w.resize(Blocks.size() + 1);
	res.tri_contrib.resize(Blocks.size() + 1); // +1 to reserve room for the remaining unmatched vertices
	
	vector<int> accumulated_mi, accumulated_mj;
	vector<int> Block_id;
	for(k = 0; k < Blocks.size(); k++) {
		for(i = 0; i < Blocks[k].size; i++) {
			accumulated_mi.push_back(Blocks[k].rows[i]);
			accumulated_mj.push_back(Blocks[k].cols[i]);
			Block_id.push_back(k);
		}
	}

	unsigned int size = accumulated_mi.size(); // Size (#nodes) of the alignment graph		
	bool **alignGraph = (bool **) calloc(size, sizeof(bool *));
	if(alignGraph == NULL) {
		fprintf(stderr, "\t\t\tcomputeWeights:: Can't allocate memory for alignGraph.\n");
		exit(-1);
	}
	for(i = 0; i < size; i++) {
		alignGraph[i] = (bool *) calloc(size, sizeof(bool));
		if(alignGraph[i] == NULL) {
			fprintf(stderr, "\t\t\tcomputeWeights:: Can't allocate memory for alignGraph[%d].\n", i);
			exit(-1);
		}
	}

	// Construct alignment graph
	register unsigned int i1, i2, j1, j2;
	for(i = 0; i < size; i++) {
		i1 = accumulated_mi[i];
		i2 = accumulated_mj[i];		
		for(j = 0; j < size; j++) {
			j1 = accumulated_mi[j];
			j2 = accumulated_mj[j];		
			if(G->getEdge(i1, j1) && H->getEdge(i2, j2)) {
				alignGraph[i][j] = true;
			}
		}		
	}
	

	for(i = 0; i <= Blocks.size(); i++) {
		res.tri_contrib[i] = 0;
	}
	
	printf("\t\tCounting triangles ... "); fflush(stdout);
	long T = 0;
	double end_point_contrib = 1.0/3;
	for(i = 0; i < size; i++) {
		for(j = i+1; j < size; j++) {
			if( !alignGraph[i][j] )
				continue;	
				
			for(k = j+1; k < size; k++) {
				if( alignGraph[i][k] && alignGraph[j][k] ) {
						T++;
						res.tri_contrib[Block_id[i]] += end_point_contrib;
						res.tri_contrib[Block_id[j]] += end_point_contrib;
						res.tri_contrib[Block_id[k]] += end_point_contrib;
				}										
			}
		}
	}	
	printf("Total triangles conserved so far = %ld\n", T);


	printf("\t\tComputing weights (rel_weight = %e) ...\n", rel_weight);
	int num_unmatched_pairs = min(unmatched_rows.size(), unmatched_cols.size());				
	if(rel_weight == 0) { // Full ortho
		for(i = 0; i < Blocks.size(); i++) {
			res.w[i] = 0;
			printf("\t\t\tres.w[%d] = %f (tries contrib = %f)\n", i, res.w[i], res.tri_contrib[i]);
		}
		res.w[i] = 1;
		printf("\t\t\tres.w[%d] = %f (tries contrib = %f)\n", i, res.w[i], res.tri_contrib[i]);
		return res;
	}
	else {
		if(rel_weight == -1) { // Adaptive weighting
			res.tri_contrib[Blocks.size()] = num_unmatched_pairs? round((double)T / Blocks.size()):0;	// If there are no more unmatched pairs, set to 0, O.W., to the average of previous blocks
			printf("\t\tEstimated # of unmatched triangles = %.2f\n", res.tri_contrib[Blocks.size()]);
			T += res.tri_contrib[Blocks.size()];
		}
		else { // Fixed weighting
			double beta = (1 - rel_weight) / rel_weight;
			res.tri_contrib[Blocks.size()] = beta*T;	
			T += res.tri_contrib[Blocks.size()];		
		}
		
		for(i = 0; i <= Blocks.size(); i++) {
			res.w[i] = res.tri_contrib[i]/T;
			printf("\t\t\tres.w[%d] = %f (tries contrib = %f)\n", i, res.w[i], res.tri_contrib[i]);
		}
	}	
	
	for(i = 0; i < size; i++) {
		free(alignGraph[i]);
	}
	free(alignGraph);	
	
	fflush(stdout);	
	return res;
}


block ProdTensor::Binarize(double *x, int method) {
	printf("\tBinarizing solution (method = %d):\n", method);
    int num_matches;
	register unsigned int i, j, k;
	block res_block;
	res_block.size = 0; 
	
	if(norm(x, n) == 0) {
		printf("Binarize:: vector has norm zero.\n");
		return res_block;
	}
	
	// Ensure one-to-one matching ...
	printf("\t\tBi-partite matching: "); fflush(stdout);
	double num_edges = setupBipartiteGraph(x); // Updates rows, cols, and edge_weights based on unmatched pairs in X
	res_block.match_score = match(n1, n2, num_edges, rows, cols, edge_weights, mi, mj, &num_matches);		
	printf("Match score = %f\n", res_block.match_score);
	
	
	// Find number of non-zeros in the matching vector ...
	printf("\t\tSelecting \"non-zeros\" from the match weights: \n"); fflush(stdout);
	vector<int> match_perm, rows, cols;
	vector<double> match_weights(num_matches), sorted_match_weights(num_matches);
	for(k = 0; k < (unsigned int)num_matches; k++) { 
		match_weights[k] = x[(int)(mi[k]) * n2 + (int)(mj[k])];
	}				
    sortingPermutation(match_weights, sorted_match_weights, match_perm);

	switch(method) {
		case 0: // Triangle stat optimization method
			res_block.size = opt_tri_stat(x, mi, mj, match_perm);	
			break;
		case -1: // No binarization-- use all remaining unmatched vertices
			res_block.size = num_matches;
			break;
		default: // non-zero estimate method
			res_block.size = nonzeros_estimate(match_weights);
			break;
	}
	if(res_block.size <= 0) {
		printf("Binarize:: None of matches were selected.\n");
		return res_block;		
	}
	res_block.locality = 1 - ((double)res_block.size / num_matches);	
	printf("\t\t# pairs = %d, locality (perc. not selected-- zeros) = %.2f %%\n", (int)res_block.size, res_block.locality*100); fflush(stdout);
	
	
	// Extract significant matches ...
    //double val;
	for(k = 0; k < res_block.size; k++) { 
		rows.push_back(mi[match_perm[k]]);
		cols.push_back(mj[match_perm[k]]);
	//	val = x[rows[k] * n2 + cols[k]];
	}				

	
	/***********************************************************
	 * 	Remove matched pairs from the list of unmatched pairs
	 **********************************************************/
    //unmatched_rows = unmatched_rows \ rows;
	res_block.rows = rows;
    sort(rows.begin(), rows.end());    
	unsigned int idx1 = 0, idx2 = 0;
    for( i = 0; i < unmatched_rows.size(); i++) {
		if(unmatched_rows[i] != rows[idx2]) {
			unmatched_rows[idx1++] = unmatched_rows[i];
		}
		else {
			idx2++;
			if(idx2 == rows.size()) {
				i++;
				break;
			}
		}	
	}
	for(; i < unmatched_rows.size(); i++) {
		unmatched_rows[idx1++] = unmatched_rows[i];		
	}
	unmatched_rows.resize(idx1);
	res_block.unmatched_rows = unmatched_rows;
	
	
	
    //unmatched_cols = unmatched_cols \ cols;
    res_block.cols = cols;
    sort(cols.begin(), cols.end());    
	idx1 = 0, idx2 = 0;
    for( i = 0; i < unmatched_cols.size(); i++) {
		if(unmatched_cols[i] != cols[idx2]) {
			unmatched_cols[idx1++] = unmatched_cols[i];
		}
		else {
			idx2++;
			if(idx2 == cols.size()) {
				i++;
				break;
			}
		}	
	}
	for(; i < unmatched_cols.size(); i++) {
		unmatched_cols[idx1++] = unmatched_cols[i];		
	}
	unmatched_cols.resize(idx1);
	res_block.unmatched_cols = unmatched_cols;


	/***********************************************
	 * 			Strore B, H, and V blocks
	 **********************************************/
	// Diagonal B block
	register int row, col;
	res_block.B = (double *)malloc(res_block.size*res_block.size*sizeof(double));
	if(res_block.B == NULL) {
		fprintf(stderr,"Binarize: unable to allocate memory for B\n");
		exit(-1);
	}
	for(i = 0; i < res_block.size; i++) {
		for(j = 0; j < res_block.size; j++) {
			row = res_block.rows[i];
			col = res_block.cols[j];
			res_block.B[i*res_block.size + j] = x[row * n2 + col];
		}		
	}

    // In order, horizontal H block
	res_block.H = (double *)malloc(res_block.size*unmatched_cols.size()*sizeof(double));
	if(res_block.H == NULL) {
		fprintf(stderr,"Binarize: unable to allocate memory for H\n");
		exit(-1);
	}
	for(i = 0; i < res_block.rows.size(); i++) {
		for(j = 0; j < unmatched_cols.size(); j++) {
			row = res_block.rows[i];
			col = unmatched_cols[j];
			res_block.H[i*unmatched_cols.size() + j] = x[row * n2 + col];
		}		
	}


    // In order, vertical v block
	res_block.V = (double *)malloc(unmatched_rows.size()*res_block.size*sizeof(double));
	if(res_block.V == NULL) {
		fprintf(stderr,"Binarize: unable to allocate memory for V\n");
		exit(-1);
	}
	for(i = 0; i < unmatched_rows.size(); i++) {
		for(j = 0; j < res_block.cols.size(); j++) {
			row = unmatched_rows[i];
			col = res_block.cols[j];
			res_block.V[i*res_block.size + j] = x[row * n2 + col];
		}		
	}
    
    
    /****************************************
     * 		Count conserved triangles 
     ***************************************/	
	res_block.tries = countTrianglesUnderAlignment(res_block.rows, res_block.cols);
	printf("\t\tConserved internal triangles = %ld\n", res_block.tries); fflush(stdout);


	return res_block;	
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

int ProdTensor::TAME(double alpha, double beta, double *w, int max_it, double epsilon, int shift_type, int init_type, int num_rounds, double locality_threshold) {
	char path[1024];
	register unsigned int i, j, k;
	block res_block;
	FILE *fd;
	clock_t startTime, endTime;
	double dt;
	
	Round = 1;
	if(num_rounds < 0 || (int)min(n1, n2) < num_rounds) {
		num_rounds = min(n1, n2); // In worst case, we add a single pair in each round
	}	


	
		
	// Remove known alignment(s) from further consideration by storing it as the first block
	if(0 < this->prior_alignments[0].size()) {
		printf("\tProcessing known alignment\n");				
		block prior_block;
		prior_block.size = this->prior_alignments[0].size(); 

		
		// Ensure that the block has norm 1
		prior_block.B = (double *)calloc(prior_block.size*prior_block.size, sizeof(double));
		if(prior_block.B == NULL) {
			fprintf(stderr,"TAME: unable to allocate memory for B\n");
			exit(-1);
		}
		double partial_val = sqrt(1/prior_block.size);
		for(i = 0; i < prior_block.size; i++) {
			prior_block.rows.push_back(prior_alignments[0][i]);
			prior_block.cols.push_back(prior_alignments[1][i]);			
			prior_block.B[i*(prior_block.size + 1)] = partial_val;
			
			printf("\t\tAligning %s <->	 %s with weight %e\n", this->G->getVertexName(prior_alignments[0][i]), this->H->getVertexName(prior_alignments[1][i]), partial_val);					
		}
		printf("\tNorm of prior block = %e\n", norm(prior_block.B, prior_block.rows.size()*prior_block.cols.size()));

		// Update unmatched vertices
		sort(this->prior_alignments[0].begin(), this->prior_alignments[0].end()); // Aligned vertices of the first graph
		unsigned int idx1 = 0, idx2 = 0;
		for( i = 0; i < unmatched_rows.size(); i++) {
			if(unmatched_rows[i] != this->prior_alignments[0][idx2]) {
				unmatched_rows[idx1++] = unmatched_rows[i];
			}
			else {
				printf("Removing %s from unmatched vertices of G\n", this->G->getVertexName(unmatched_rows[i]));
				idx2++;
				if(idx2 == this->prior_alignments[0].size()) {
					i++;
					break;
				}
			}
		}
		for(; i < unmatched_rows.size(); i++) {
			unmatched_rows[idx1++] = unmatched_rows[i];
		}
		unmatched_rows.resize(idx1);
		prior_block.unmatched_rows = unmatched_rows;


		// Now remove aligned vertices of the second graph from the unmatched_cols
		sort(this->prior_alignments[1].begin(), this->prior_alignments[1].end()); 
		idx1 = 0; idx2 = 0;
		for( i = 0; i < unmatched_cols.size(); i++) {
			if(unmatched_cols[i] != this->prior_alignments[1][idx2]) {
				unmatched_cols[idx1++] = unmatched_cols[i];
			}
			else {
				printf("Removing %s from unmatched vertices of H\n", this->H->getVertexName(unmatched_cols[i]));
				idx2++;
				if(idx2 == this->prior_alignments[1].size()) {
					i++;
					break;
				}
			}
		}
		for(; i < unmatched_cols.size(); i++) {
			unmatched_cols[idx1++] = unmatched_cols[i];
		}
		unmatched_cols.resize(idx1);		 
		prior_block.unmatched_cols = unmatched_cols;


		// Reset H& V in the prior block
		prior_block.H = (double *)calloc(prior_block.size*unmatched_cols.size(), sizeof(double));
		if(prior_block.H == NULL) {
			fprintf(stderr,"TAME: unable to allocate memory for H\n");
			exit(-1);
		}
		
		prior_block.V = (double *)calloc(unmatched_rows.size()*prior_block.size, sizeof(double));
		if(prior_block.V == NULL) {
			fprintf(stderr,"TAME: unable to allocate memory for V\n");
			exit(-1);
		}
		Blocks.push_back(prior_block);
	}
	

	this->res = new eigen;	
	res->x = (double *) calloc(this->n*this->subspace_dim, sizeof(double)) ;
	if(res->x == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for res->x\n");
		exit(-1);
	}		
	
	res->trace = (double *) calloc(max_it, sizeof(double)); 
	if(res->trace == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for res->trace\n");
		exit(-1);
	}		
	

	// Initialize X0
	double *x0 = (double *)calloc(n, sizeof(double));
	if(x0 == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for 'x0'\n");
		return -1;
	}
	InitX(x0, init_type, w);		
	
	// Estimate relative weight of seqsim/tri match (rho)
	if(w != NULL) {
		if(normalize(w, w, n)  == 0)  { 
		  fprintf(stderr, "tame:: w vector has norm 0.\n");	  
		}
		
		// Computing relative weight of triangular matching and node similarity scores
		double mean_seqsim = 0, mean_uniTris = 0, epsilon = 1.0 / n;
		unsigned int *Ni1, *Ni2;
		
		unsigned int G_trie_deg, H_trie_deg;
		for(i = 0; i < n1; i++) {
			Ni1 = T_G->getAdjTriangles(i);
			G_trie_deg = Ni1[0];
			for(j = 0; j < n2; j++) {
				Ni2 = T_H->getAdjTriangles(j);
				H_trie_deg = Ni2[0];
				mean_uniTris += 12*epsilon*G_trie_deg*H_trie_deg; // Each triangle contributes 3epsilon, there are 2 permutations (jj', kk' vs jk', kj'), and there is a factor of 2 for symmetry (jj', kk' vs kk', jj')
				mean_seqsim += w[n2*i+j];
			}
		}
		this->rho = mean_uniTris / mean_seqsim; // Weight of SeqSim
		printf("SeqSim weight = %e\tTri weight = %e\tSeqSim Normalization factor = %e\n", mean_seqsim / n, mean_uniTris / n, rho);
	}
	else {
		this->rho = 0;
	}




	printf("Round 1 (initial phase):\n\t");
	sprintf(path, "%s/%s%s_%.2e_round001.mat", output_path, prefix, ShiftType_names[shift_type], alpha);		
	if(exist(path)) { // If we have pre-computed results of SS-HOPM and could successfully load it
		if(readVec(path, res->x, n1, n2) != 0) {
			fprintf(stderr, "Error reading file %s\n", path);
			return -1;
		}		
	}
	else { // Rerun SS-HOPM
		printf("File %s does not exist. Running SS-HOPM ...\n", path);
		startTime = clock();	
		this->unmatched_weight = 1;	
		if(shift_type != MULTI_LINEAR) {
			res = issHOPM(subspace_dim, max_it, alpha, beta, epsilon, w, x0, init_type, shift_type);
		}
		else {
			res = MLPR(alpha, w,  x0,  max_it,  epsilon, true);			
		}
		endTime = clock();	
		dt = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;	
		printf("\tElapsed time %lf (s)\n", dt); 		

		// Exort original matrix& trace
		printf("\tExporting initial phase results ... ");
		sprintf(path, "%s/%s%s_%.2e_round%03d.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
		fd = fopen(path, "w");
		printVec(fd, res->x, n1, n2);
		fclose(fd);			

					
		sprintf(path, "%s/%s%s_%.2e_round%03d_trace.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
		fd = fopen(path, "w");
		printVec(fd, res->trace, res->its, 1);
		fclose(fd);					
		printf("done\n");
	}

	/*******************************************
	 * 			Process eigenvector
	 ******************************************/
	res_block = Binarize(res->x, this->bin_type);
	if(res_block.size != 0) {
		Blocks.push_back(res_block); // Store B_k, V_k, H_k, P_k (rows), and Q_k (cols)

		sprintf(path, "%s/%s%s_%.2e_round%03d_B.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
		fd = fopen(path, "w");
		printVec(fd, res_block.B, res_block.size, res_block.size);
		fclose(fd);					

		sprintf(path, "%s/%s%s_%.2e_round%03d_H.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
		fd = fopen(path, "w");
		printVec(fd, res_block.H, res_block.size, unmatched_cols.size());
		fclose(fd);					

		sprintf(path, "%s/%s%s_%.2e_round%03d_V.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
		fd = fopen(path, "w");
		printVec(fd, res_block.V, unmatched_rows.size(), res_block.size);
		fclose(fd);		
	}




	Round ++;	
	block_stats Stats;
	double block_norm, B_norm, V_norm, H_norm;
	register int row, col;
	for (; Round <= (unsigned int) num_rounds; Round++) {		
		/*******************************************
		 *   Decide if we should run more rounds
		 ******************************************/
		if( unmatched_rows.size() == 0 || //**** TODO: DECIDE how to stop rounds with a few unmatched pairs left ... : If SeqSim tries were greater
			unmatched_cols.size() == 0) {
				printf("=> Terminating: All pairs are matched\n");
				Round --;
				break;
		} 
/*		else if(Blocks[Blocks.size()-1].locality < locality_threshold) { // Stop if eigenvector de-localizes
				printf("=> Terminating: Eigen-value locality = %.2f %% < %.2f %%\n", Blocks[Blocks.size()-1].locality, locality_threshold);
				Round --;
				break;		
		}
	*/	
		//else if{} // Additional stopping criterias
		else {
			printf("Round %d ... \n", Round);

			/************************************************
			* 				Update x0
			***********************************************/	
			printf("\tInitialize x0 ... \n"); fflush(stdout);			


			
			printf("\t\tComputing weight of each block ... \n");
			/* 1 - this->relative_weight would be the weight of unmatched block
			Stats.[0..Blocks.size()-1]: weight of each block, Stats.w[Blocks.size()]: weight of remaining unmatched block, sum(Stats.w) = 1; */
			Stats = computeWeights(this->relative_weight); 

			// First we have to reset x0 to prior vector over unmatched pairs and scale ...
			InitX(x0, init_type, w); // x0 has norm 1 here. -- For elements matched in previous blocks, x0 will be 0 here.
			double weight = sqrt(Stats.w[Blocks.size()]); 
			printf("\t\tRescaling  x0 with %.2f estimated triangle contribution (weight = %e, scale = %e) ... \n", Stats.tri_contrib[Blocks.size()], Stats.w[Blocks.size()], weight);				
			for(i = 0; i < n; i++) {
				x0[i] *= weight;
			}
			printf("\t\t\tnorm of initial x0 = %e\n", norm(x0, n));
			
			
			// Now merge previous blocks into x0
			printf("\t\tMerging previous blocks into x0:\n");
			for(k = 0; k < Blocks.size(); k++) {
				B_norm = norm(Blocks[k].B, Blocks[k].rows.size()*Blocks[k].cols.size());
				H_norm = Blocks[k].unmatched_cols.size() == 0? 0:norm(Blocks[k].H, Blocks[k].rows.size()*Blocks[k].unmatched_cols.size());
				V_norm = Blocks[k].unmatched_rows.size() == 0? 0:norm(Blocks[k].V, Blocks[k].unmatched_rows.size()*Blocks[k].cols.size());												
				block_norm =  sqrt(B_norm*B_norm + H_norm*H_norm + V_norm*V_norm); //||y||_2				
				weight = sqrt(Stats.w[k]) / block_norm;

				printf("\t\t\tMerging Block %d with %.2f triangle contribution (weight = %e, scale = %e)\n", k, Stats.tri_contrib[k], Stats.w[k], weight);								
				for(i = 0; i < Blocks[k].rows.size(); i++) { // B_k
					row = Blocks[k].rows[i];
					for(j = 0; j < Blocks[k].cols.size(); j++) {
						col = Blocks[k].cols[j];
						x0[row * n2 + col] = weight*Blocks[k].B[i*Blocks[k].cols.size() + j];
					}
				}
								
				for(i = 0; i < Blocks[k].rows.size(); i++) { // H_k
					row = Blocks[k].rows[i];
					for(j = 0; j < Blocks[k].unmatched_cols.size(); j++) {
						col = Blocks[k].unmatched_cols[j];						
						x0[row * n2 + col] = weight*Blocks[k].H[i*Blocks[k].unmatched_cols.size() + j];
					}
				}				

				for(i = 0; i < Blocks[k].unmatched_rows.size(); i++) { // V_k
					row = Blocks[k].unmatched_rows[i];
					for(j = 0; j < Blocks[k].cols.size(); j++) {
						col = Blocks[k].cols[j];						
						x0[row * n2 + col] = weight*Blocks[k].V[i*Blocks[k].cols.size() + j];
					}
				}	
				printf("\t\t\tNorm of current x0 = %f\n", norm(x0, n));							
			}
			printf("\tFinal norm of x0 = %f\n", norm(x0, n)); 


			startTime = clock();
			printf("\tRunning SS-HOPM ...\n"); fflush(stdout);			
			this->unmatched_weight = sqrt(Stats.w[Blocks.size()]); // Ensures that the relative contribution of unmatched pairs to the norm of X remains the same
			if(shift_type != MULTI_LINEAR) {
				res = issHOPM(subspace_dim, max_it, alpha, beta, epsilon, w, x0, init_type, shift_type);
			}
			else {
				res = MLPR(alpha, w,  x0,  max_it,  epsilon, true);			
			}
					
			endTime = clock();	
			dt = difftime(endTime, startTime)/(double)CLOCKS_PER_SEC;	
			printf("\tElapsed time %lf (s)\n", dt); 							


			// Exort original matrix& trace
			printf("\tExporting X of round %d ... ", Round);
			sprintf(path, "%s/%s%s_%.2e_round%03d.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
			fd = fopen(path, "w");
			printVec(fd, res->x, n1, n2);
			fclose(fd);			
			sprintf(cmd, "gzip %s -f", path);
			if(system(cmd) == -1) {
				fprintf(stderr, "issHOPM:: Error compressing %s", path);
			}
			
			sprintf(path, "%s/%s%s_%.2e_round%03d_trace.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
			fd = fopen(path, "w");
			printVec(fd, res->trace, res->its, 1);
			fclose(fd);					
			printf("done\n");			
			
			/*******************************************
			 * 			Process eigenvector
			 ******************************************/
			res_block = Binarize(res->x, this->bin_type);
			if(res_block.size != 0) {
				Blocks.push_back(res_block); // Store B_k, V_k, H_k, P_k (rows), and Q_k (cols)

				sprintf(path, "%s/%s%s_%.2e_round%03d_B.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
				fd = fopen(path, "w");
				printVec(fd, res_block.B, res_block.size, res_block.size);
				fclose(fd);					

				sprintf(path, "%s/%s%s_%.2e_round%03d_H.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
				fd = fopen(path, "w");
				printVec(fd, res_block.H, res_block.size, unmatched_cols.size());
				fclose(fd);					

				sprintf(path, "%s/%s%s_%.2e_round%03d_V.mat", output_path, prefix, ShiftType_names[shift_type], alpha, Round);
				fd = fopen(path, "w");
				printVec(fd, res_block.V, unmatched_rows.size(), res_block.size);
				fclose(fd);					
			}
			else {
				printf("=> Terminating: Empty block after binarization\n");
				Round --;
				break;		
			}
		}
	}	
	fflush(stdout);
	
	unsigned int *P = (unsigned int *) calloc(n1, sizeof(unsigned int)); // Row permutations
	if(P == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for P\n");
		exit(-1);		
	}
	unsigned int *Q = (unsigned int *) calloc(n2, sizeof(unsigned int)); // Column permutations
	if(Q == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for Q\n");
		exit(-1);		
	}

	printf("\tEstimate final weight of each block\n");
	res_block = Binarize(res->x, -1); // Extract sorted block for all remaining unmatched vertices
	Blocks.push_back(res_block);	
	Stats = computeWeights(1.0);
					
					

	// Now merge previous blocks into x0
	printf("\tMerging blocks into result vector:\n");
	bzero(x0, n*sizeof(double));			

	int match_counter = 0;	
	for(k = 0; k < Blocks.size(); k++) {		
		B_norm = norm(Blocks[k].B, Blocks[k].rows.size()*Blocks[k].cols.size());
		/*
		H_norm = Blocks[k].unmatched_cols.size() == 0? 0:norm(Blocks[k].H, Blocks[k].rows.size()*Blocks[k].unmatched_cols.size());
		V_norm = Blocks[k].unmatched_rows.size() == 0? 0:norm(Blocks[k].V, Blocks[k].unmatched_rows.size()*Blocks[k].cols.size());
		block_norm =  sqrt(B_norm*B_norm + H_norm*H_norm + V_norm*V_norm); //||y||_2
		*/
		block_norm =  sqrt(B_norm*B_norm);
		double weight = sqrt(Stats.w[k]) / block_norm;
		
		printf("\t\tMerging Block %d with %.2f triangle contribution (weight = %e, scale = %e)\n", k, Stats.tri_contrib[k], Stats.w[k], weight);				
		for(i = 0; i < Blocks[k].rows.size(); i++) { // B_k
			P[match_counter] = Blocks[k].rows[i];
			Q[match_counter] = Blocks[k].cols[i];
			match_counter++;					
			row = Blocks[k].rows[i];
			for(j = 0; j < Blocks[k].cols.size(); j++) {
				col = Blocks[k].cols[j];
				x0[row * n2 + col] = weight*Blocks[k].B[i*Blocks[k].cols.size() + j];
			}
		}
						
/*		for(i = 0; i < Blocks[k].rows.size(); i++) { // H_k
			row = Blocks[k].rows[i];
			for(j = 0; j < Blocks[k].unmatched_cols.size(); j++) {
				col = Blocks[k].unmatched_cols[j];						
				x0[row * n2 + col] = weight*Blocks[k].H[i*Blocks[k].unmatched_cols.size() + j];
			}
		}				

		for(i = 0; i < Blocks[k].unmatched_rows.size(); i++) { // V_k
			row = Blocks[k].unmatched_rows[i];
			for(j = 0; j < Blocks[k].cols.size(); j++) {
				col = Blocks[k].cols[j];						
				x0[row * n2 + col] = weight*Blocks[k].V[i*Blocks[k].cols.size() + j];
			}
		}*/
		printf("\t\t\tNorm of current result vector = %f\n", norm(x0, n));
	}

	printf("\tFinal norm of result vector = %f\n", norm(x0, n));

						
	
	sprintf(path, "%s/%s%s_%.2e_X_final.mat", output_path, prefix, ShiftType_names[shift_type], alpha);	
	printf("\tWriting final X to %s\n", path);
	fd = fopen(path, "w");
	if(fd == NULL) {
		fprintf(stderr, "\t\tTAME:: Can't open %s\n", path);
		exit(-1);
	}	
	printVec(fd, x0, n1, n2);
	fclose(fd);		
	
	stats myStats = evalX(x0);
	
	printf("Final conserved stats: edges = %ld, triangles = %ld\n", myStats.M, myStats.T);



	sprintf(path, "%s/%s%s_%.2e_P_final.mat", output_path, prefix, ShiftType_names[shift_type], alpha);
	printf("\tWriting final row permutation to %s\n", path);
	fd = fopen(path, "w");
	if(fd == NULL) {
		fprintf(stderr, "\t\tTAME:: Can't open %s\n", path);
		exit(-1);
	}	
	for(i = 0; i < n1; i++) {
		fprintf(fd, "%d\n", P[i]);
	}
	fclose(fd);		
	free(P);


	sprintf(path, "%s/%s%s_%.2e_Q_final.mat", output_path, prefix, ShiftType_names[shift_type], alpha);
	printf("\tWriting final column permutation to %s\n", path);
	fd = fopen(path, "w");
	if(fd == NULL) {
		fprintf(stderr, "\t\tTAME:: Can't open %s\n", path);
		exit(-1);
	}
	for(i = 0; i < n2; i++) {
		fprintf(fd, "%d\n", Q[i]);
	}
	fclose(fd);			
	free(Q);



	free(x0);	
	free(res->x);	
	free(res->trace);
	delete this->res;	
	
	return Round;
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
					
