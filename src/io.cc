#include <io.h>

/* Reads network file in IsoRank format
 * Input) path: Path to the input file in IsoRank format, g: Pointer to an unitialized graph data structure
 * Returns) It initializes g, and returns 0 on success of -1 on error
*/ 
int readISO(char *path, Graph *&g) {
	int res;
	
	FILE *fd = fopen(path, "r");	
	if(fd == NULL) {
		fprintf(stderr, "readISO:: Can't open %s\n", path);
		return -1;
	}
	printf("Reading IsoRank network from %s ...\n", path);

	unsigned int n, ne;	
	unsigned short int *adj;	
	vector<string> id2gene;	
	map <string, int> gene2id;
			
	register int vertex_id = 0;
	char src[1024], dst[1024];
	unsigned long line_no = 0;

	
	//Skip header line
	res = fscanf(fd, "%*s %*s");
	
	
	// Read edges, and construct a unified list of vertex id(s)
	while(!feof(fd)) {
		res = fscanf(fd, "%1024s %1024s", src, dst);
		
		// No more edges
		if(feof(fd))
			break;
		
		if(res != 2) {
			fprintf(stderr, "Syntax errror on line %ld\n", line_no);
			return -1;
		}
		
		if(gene2id.find(src) == gene2id.end()) {
			gene2id[src] = vertex_id++;
			id2gene.push_back(src);
		}
		if(gene2id.find(dst) == gene2id.end()) {
			gene2id[dst] = vertex_id++;			
			id2gene.push_back(dst);
		}
		line_no++;
	}
	n = gene2id.size();
	

	long mem_size = n*n*sizeof(unsigned short int);
	adj = (unsigned short int*)calloc(1, mem_size);
	if(adj == NULL) {
		fprintf(stderr, "readISO:: Error allocating memory for adjacency matrix.\n");
		return -1;
	}
	
	fseek(fd, 0, SEEK_SET);
	res = fscanf(fd, "%*s %*s");
	
	for(register unsigned int i = 0; i < line_no; i++) {
		res = fscanf(fd, "%1024s %1024s", src, dst);
		
		if(!strcmp(src, dst)) // Remove self-loops
			continue;
		
		// Ensure bi-directionality of edges
		adj[gene2id[src]*n + gene2id[dst]] = 1;
		adj[gene2id[dst]*n + gene2id[src]] = 1;
	}
	
	ne = 0;
	for(register unsigned int i = 0; i < n; i++) {
		for(register unsigned int j = i+1; j < n; j++) {
			if(adj[i*n + j] == 1)
				ne++;
		}		
	}
	
	printf("\tNumber of lines == %ld\n", line_no);
	printf("\tNumber of nodes = %d\n", n);
	printf("\tNumber of unique edges = %d\n", ne);
		
	printf("Done\n");
	fclose(fd);
	
	g = new Graph(n, ne, adj, id2gene);	
	free(adj);	
	return 0;
}


/* Reads network file as an edge-list in Tab-separated format
 * Input) path: Path to the input file in IsoRank format, g: Pointer to an unitialized graph data structure
 * Returns) It initializes g, and returns 0 on success of -1 on error
*/ 
int readTab(char *path, Graph *&g) {
	int res;
	
	FILE *fd = fopen(path, "r");	
	if(fd == NULL) {
		fprintf(stderr, "readTab:: Can't open %s\n", path);
		return -1;
	}
	printf("Reading Tab-separated network from %s ...\n", path);

	unsigned int n, ne;	
	unsigned short int *adj;	
	vector<string> id2gene;	
	map <string, int> gene2id;
			
	register int vertex_id = 0;
	char src[1024], dst[1024];
	unsigned long line_no = 0;

	
	// Read edges, and construct a unified list of vertex id(s)
	while(!feof(fd)) {
		res = fscanf(fd, "%1024s %1024s", src, dst);
		
		// No more edges
		if(feof(fd))
			break;
		
		if(res != 2) {
			fprintf(stderr, "Syntax errror on line %ld\n", line_no);
			return -1;
		}
		
		if(gene2id.find(src) == gene2id.end()) {
			gene2id[src] = vertex_id++;
			id2gene.push_back(src);
		}
		if(gene2id.find(dst) == gene2id.end()) {
			gene2id[dst] = vertex_id++;			
			id2gene.push_back(dst);
		}
		line_no++;
	}
	n = gene2id.size();
	

	long mem_size = n*n*sizeof(unsigned short int);
	adj = (unsigned short int*)calloc(1, mem_size);
	if(adj == NULL) {
		fprintf(stderr, "readTAB:: Error allocating memory for adjacency matrix.\n");
		return -1;
	}
	
	fseek(fd, 0, SEEK_SET);	
	for(register unsigned int i = 0; i < line_no; i++) {
		res = fscanf(fd, "%1024s %1024s", src, dst);
		
		if(!strcmp(src, dst)) // Remove self-loops
			continue;
		
		// Ensure bi-directionality of edges
		adj[gene2id[src]*n + gene2id[dst]] = 1;
		adj[gene2id[dst]*n + gene2id[src]] = 1;
	}
	
	ne = 0;
	for(register unsigned int i = 0; i < n; i++) {
		for(register unsigned int j = i+1; j < n; j++) {
			if(adj[i*n + j] == 1)
				ne++;
		}		
	}
	
	printf("\tNumber of lines == %ld\n", line_no);
	printf("\tNumber of nodes = %d\n", n);
	printf("\tNumber of unique edges = %d\n", ne);
		
	printf("Done\n");
	fclose(fd);
	
	g = new Graph(n, ne, adj, id2gene);	
	free(adj);	
	return 0;
}







/* Reads network file sparse matrix (SMAT) format
 * Input) path: Path to the input file in IsoRank format, g: Pointer to an unitialized graph data structure
 * Returns) It initializes g, and returns 0 on success of -1 on error
*/ 
int readSMAT(char *path, Graph *&g) {
	int res;
	register unsigned int i, j;
	
	FILE *fd = fopen(path, "r");	
	if(fd == NULL) {
		fprintf(stderr, "readTab:: Can't open %s\n", path);
		return -1;
	}
	printf("Reading SMAT network from %s ...\n", path);

	unsigned int n, ne;	
	res = fscanf(fd,"%d %*d %d", &n, &ne);
	printf("# vertices = %d\t# edges = %d\n", n, ne);


	long mem_size = n*n*sizeof(unsigned short int);
	unsigned short int *adj = (unsigned short int*)calloc(1, mem_size);
	if(adj == NULL) {
		fprintf(stderr, "readTAB:: Error allocating memory for adjacency matrix.\n");
		return -1;
	}
	
	char label[1024];
	vector<string> id2gene;		
	for(i = 0; i < n; i++) {
		sprintf(label, "%d", i);
		id2gene.push_back(label);
	}	
	
	
	// Read non-zeros
	unsigned long line_no = 0;	
	while(!feof(fd)) {
		res = fscanf(fd, "%d %d %*f", &i, &j);
		
		// No more edges
		if(feof(fd))
			break;
		
		if(res != 2) {
			fprintf(stderr, "Syntax errror on line %ld\n", line_no);
			return -1;
		}
		
		if(i == j) // remove self-loops
			continue;
			
		adj[i*n + j] = adj[j*n + i] = 1;
		line_no++;
	}
	
		
	printf("Done\n");
	fclose(fd);
	
	g = new Graph(n, ne, adj, id2gene);	
	free(adj);	
	return 0;
}










/* Reads sequence similarity scores for protein pairs in input networks
 * Input) w: output vector to be initialized, path: Path to the input file in eval format, g and h: Initialized input graphs
 * Returns) It initializes w (un-normalized), and returns 0 on success of -1 on error
*/ 
double *readSeqSim(Graph *&g, Graph *&h, char *path, char *out_path) {
	int res;
	double *w = NULL;
	if(strcmp(path, "") == 0)
		return NULL;
		
	FILE *fd = fopen(path, "r");	
	if(fd == NULL) {
		fprintf(stderr, "readSeqSim:: Can't open %s\n", path);
		return NULL;
	}
	printf("Reading sequence similarities from %s ...\n", path);

	int mem_size = g->n*h->n*sizeof(double);
	w = (double *)calloc(1, mem_size);
	if(w == NULL) {
		fprintf(stderr, "readSeqSim:: Error allocating memory for w.\n");
		return NULL;
	}
	
	
	double val;
	register int row = -1, col = -1;
	char src[1024], dst[1024];
	long line_no=0, matched = 0;
	while(!feof(fd)) {
		res = fscanf(fd, "%s %s %lf", src, dst, &val);
		if(res == EOF)
			break;
		
		line_no++;
		if( (row = g->getVertexID(src)) != -1) {
			if( (col = h->getVertexID(dst)) != -1) {
				matched++;
				w[h->n*row + col] = val;
			}						
		}		
	}

	fclose(fd);		
	printf("\t%ld lines parsed successfully (matched lines = %ld)\n", line_no, matched);
	

	printf("\tExporting results to %s ...\n", out_path);
	FILE *out_fd = fopen(out_path, "w");
	if(out_fd == NULL) {
		fprintf(stderr, "readSeqSim:: Can't open %s for exporting the results\n", out_path);
		return NULL;
	}
	
	fprintf(out_fd, "%d\t%d\t%ld\n", g->n, h->n, matched);
	for(register unsigned int i = 0; i < g->n; i++)
		for(register unsigned int  j = 0; j < h->n; j++)
			if(w[h->n*i + j] > 0)
				fprintf(out_fd, "%d\t%d\t%f\n", i, j, w[h->n*i + j]);			

	fclose(out_fd);

	printf("Done.\n");
	return w;
}



double *readSeqSim_SMAT(char *path) {
	int res;
	double *w = NULL;
	if(strcmp(path, "") == 0)
		return NULL;
		
	FILE *fd = fopen(path, "r");	
	if(fd == NULL) {
		fprintf(stderr, "readSeqSim_SMAT:: Can't open %s\n", path);
		return NULL;
	}
	printf("Reading sequence similarities from SMAT file %s ...\n", path);

	unsigned int n, m, ne;	
	res = fscanf(fd,"%d %d %d", &m, &n, &ne);
	printf("# m = %d\tn=%d\tnnz=%d\n", m, n, ne);
	
	int mem_size = m*n*sizeof(double);
	w = (double *)calloc(1, mem_size);
	if(w == NULL) {
		fprintf(stderr, "readSeqSim:: Error allocating memory for w.\n");
		return NULL;
	}
	
	
	double val;
	register int row = -1, col = -1;
	
	long line_no=0;
	while(!feof(fd)) {
		res = fscanf(fd, "%d %d %lf", &row, &col, &val);
		if(res == EOF)
			break;
		
		line_no++;
		w[n*row + col] = val;
	}

	fclose(fd);		
	printf("\t%ld lines parsed successfully\n", line_no);

	printf("Done.\n");
	return w;
}

/* Exports newtork to SMAT format
 * Input) path: Output path, name: Output filename, g: Pointer to an initialized graph data structure
 * Returns) Writes adjacency matrix of graph "g" to "path/name.smat" and its node labelings to "path/name.attr". Returns 0 on success of -1 on error
*/ 
int writeSMAT(char *path, char *name, Graph *g) {
	if(g == NULL) {
		fprintf(stderr, "writeSMAT:: Uninitialized graph\n");
		return -1;
	}
	
	// Exporting adjacency matrix
	char smat_path[1024];
	sprintf(smat_path, "%s/%s.smat", path, name);
	FILE *fd = fopen(smat_path, "w");	
	if(fd == NULL) {
		fprintf(stderr, "writeSMAT:: Can't open %s\n", smat_path);
		return -1;
	}
	printf("Exporting adjacency matrix to %s...\n", smat_path);
	
	fprintf(fd, "%d\t%d\t%d\n", g->n, g->n, g->ne);	
	for(register unsigned int i = 0; i < g->n; i++) {
		for(register unsigned int j = i; j < g->n; j++) {
			if(g->getEdge(i, j) == 1) {
				fprintf(fd, "%d\t%d\t1\n", i, j);
			}
		}
	}
	fclose(fd);

	// Exporting node labelings
	char attr_path[1024];
	sprintf(attr_path, "%s/%s.attr", path, name);
	fd = fopen(attr_path, "w");	
	if(fd == NULL) {
		fprintf(stderr, "writeSMAT:: Can't open %s\n", attr_path);
		return -1;
	}
	printf("Exporting node attributes to %s...\n", attr_path);

	for(register unsigned int i = 0; i < g->n; i++) {
		fprintf(fd, "%s\n", g->getVertexName(i));
	}
	fclose(fd);
		
	printf("Done\n");
	return 0;
}


/* Exports list of triangles to SSTEN format
 * Input) path: Output path, name: Output filename, g: Pointer to an initialized graph data structure, t: Pointer to a preprocessed triangle count data structure
 * Returns) Writes list of triangles to "path/name.sten". Returns 0 on success of -1 on error
*/ 
int writeSSTEN(char *path, char *name, Graph *g, TriangleCounts *t) {
	if(g == NULL) {
		fprintf(stderr, "writeSSTEN:: Uninitialized graph\n");
		return -1;
	}

	if(t == NULL) {
		fprintf(stderr, "writeSSTEN:: Uninitialized TriangleCount\n");
		return -1;
	}
		
	// Exporting adjacency matrix
	char sten_path[1024];
	sprintf(sten_path, "%s/%s.sten", path, name);
	FILE *fd = fopen(sten_path, "w");	
	if(fd == NULL) {
		fprintf(stderr, "writeSSTEN:: Can't open %s\n", sten_path);
		return -1;
	}
	printf("Exporting list of triangles to %s...\n", sten_path);

	fprintf(fd, "%d\t%d\n", 3, g->n);
	vector< vector<unsigned int> > myTriangles = t->getTriangles();
	for(register unsigned int i = 0; i < myTriangles.size(); i++) {				
		fprintf(fd, "%d\t%d\t%d\t1\n", myTriangles[i][0], myTriangles[i][1], myTriangles[i][2]);
	}	
	
	fclose(fd);
	return 0;
}


void printTriangles(Graph *G, TriangleCounts &T_G) {
	vector< vector<unsigned int> > myTriangles = T_G.getTriangles();
	for(register unsigned int i = 0; i < T_G.size(); i++) {				
		printf("\t\t\t%d\t%d\t%d (<%s, %s, %s>)\n", myTriangles[i][0], myTriangles[i][1], myTriangles[i][2],
			G->getVertexName(myTriangles[i][0]), G->getVertexName(myTriangles[i][1]), G->getVertexName(myTriangles[i][2]));
	}
}

void printVec(FILE *fd, double *v, int n1, int n2) {
	register int i, j;  
	if(n1 == 0 || n2 == 0)
		return;
		
	for (i = 0; i < n1; i++) {
		fprintf(fd, "%le", v[i*n2]);
		for(j = 1; j < n2; j++) {	
			fprintf(fd, "\t%le", v[i*n2 + j]);
		}
		fprintf(fd, "\n");
	}
}

int readVec(const char *path, double *v, int n1, int n2) {
	FILE *fd = fopen(path, "r");
	if(fd == NULL) {
		fprintf(stderr, "readISO:: Can't open %s\n", path);
		return -1;
	}
	printf("Reading vector from %s ...\n", path);

	register int i, j;  

	for (i = 0; i < n1; i++) {
		for(j = 0; j < n2; j++) {	
			if(fscanf(fd, "%le", v + i*n2 + j) == EOF)
				return -1;
		}
	}
	
	fclose(fd);
	return 0;
}

int exist(char *name) {
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}

float getTotalSystemMemory() { // Return memory in MB	
	unsigned long long pages, page_size;
	pages = sysconf(_SC_AVPHYS_PAGES);	
	page_size = sysconf(_SC_PAGESIZE);
		
/*	pages = sysconf(_SC_PHYS_PAGES);
    page_size = sysconf(_SC_PAGE_SIZE);*/
    
    return ((float)(pages * page_size)) / (1<<20);
}
