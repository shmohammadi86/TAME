#include <generics.h>
#include <graph.h>
#include <tensor.h>
#include <io.h>

//#define DEBUG
//#define VERBOSE_DEBUG

char *FileType_names[]  = {"IsoRank", "Tab", "SMAT"};
char *ShiftType_names[] = {"Affine-SSHOPM", "NormalizedAffine-SSHOPM", "SSHOPM", "CombinedShift-SSHOPM", "MultiLinear"};
char *InitType_names[]  = {"Uniform", "Random", "SeqSim"};

void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: tri-match -G <1st_net> -H <2nd_net> [optional args]\n ");
    fprintf (stream,
		 "\t-h	--help\t\t\t\tDisplay this usage information. \n"
		 "\t-t	--type <tab|isorank|smat>\tType of input files (default=tab)\n"
		 "\t-G	--first <1st_net>\t\tPath to the first input network.\n"
		 "\t-H	--second <1st_net>\t\tPath to the second input network.\n"
		 "\t-S	--seq <seq_file>\t\tPath to the SeqSim file(default="").\n" 
		 "\t-o	--output <output_folder>\tOutput folder (default=output/).\n"

		 "\t-Y	--sparsity <level>\t\tSparsity level (0, 1, and 2).\n"

    	 "\t-C	--continue <X.txt>\t\tContinues iterations starting from matrix stored in file.\n"
		 "\t-x 	--x0 <uniform|random|seqsim>\tInitialization type (default=uniform).\n" 
		 "\t-i 	--iter <its>\t\t\tNumber of iterations (default=100).\n"

		 "\t-a 	--alpha <weight>\t\trelative weight of SeqSim/Triangles (default = 1.0).\n"
		 "\t-b 	--beta <shift_val>\t\tShift value (default = 0.0).\n" 
		 "\t-e 	--epsilon <eps>\t\t\tepsilon threshold (default=1e-16).\n" 		
		 
		 "\t-p 	--b_seq <bSeq>\t\t\tb-matching degree for sequences (default=50).\n" 		
		 "\t-q 	--b_topo <bSeq>\t\t\tb-matching degree for topolotgy (default=200).\n" 		
		 "\t-I 	--post_iter <iter_count>\tnumber of iterations for post-processing step (default=10).\n" 		
		  		 
		 
	);	     
    exit (exit_code);
}

int main(int argc, char **argv) {
	
	

  int next_option;
  const char *const short_options = "hG:H:t:S:o:x:i:a:b:e:Y:C:I:p:q:";
  const struct option long_options[] = {
		{"help",   0, NULL, 'h'},
		{"first",  1, NULL, 'G'},
		{"second", 1, NULL, 'H'},
		{"type", 1, NULL, 	't'},
		{"seq",   1, NULL, 	'S'},		
		{"output",	1, NULL,'o'},
		{"x0",   1, NULL, 'x'},		
		{"iter",   1, NULL, 'i'},		
		{"alpha",   1, NULL, 'a'},
		{"beta",	1, NULL, 'b'},
		{"epsilon",   1, NULL, 'e'},		
		{"sparsity", 1, NULL, 'Y'},
		{"continue", 1, NULL, 'C'},
		{"b_seq", 1, NULL, 'p'},
		{"b_topo", 1, NULL, 'q'},
		{"post_iter", 1, NULL, 'I'},
		{NULL,     0, NULL,  0 }		
	};

	char first[1024] = "", second[1024] = "", output_path[1024] = "output", prior_output[1024] = "";
	FileType file_type = Tab;
	InitType init_type = SeqSim;
	Sparsity_Type sparsity_type = NoSparsity;
	
	double alpha = 0, beta = 1;
	int max_it = 100;
	double epsilon = 1e-06;
	double val;
	char full_path1[1024] = "", full_path2[1024] = "", seq_path[1024] = "";	
	
	int b_seq=50, b_topo=200, post_iter = 10;
	
	if(argc < 3)
   		print_usage (stdout, -1);
	
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'G':
				strcpy(full_path1, optarg);
	    		break;

			case 'H':
				strcpy(full_path2, optarg);
	    		break;

			case 'o':
				strcpy(output_path, optarg);
	    		break;

			case 'C':
				strcpy(prior_output, optarg);
	    		break;

			case 't':
				if(!strcmp(optarg, "isorank")) {
					file_type = IsoRank;
				}
				else if(!strcmp(optarg, "tab")) {
					file_type = Tab;
				}
				else if(!strcmp(optarg, "smat")) {
					file_type = SMAT;
				}
				else {
					fprintf(stderr, "Main:: Unknown input filetype: %s\n", optarg);
					exit(-1);
				}
	    		break;

			case 'x':
				if( !strcasecmp(optarg, "random") ) {
					init_type = Random;
				}
				else if( !strcasecmp(optarg, "uniform") ) {
					init_type = Uniform;
				}
				else if( !strcasecmp(optarg, "seqsim") ) {
					init_type = SeqSim;
				}
				else {
					fprintf(stderr, "Main:: Unknown initialization type: %s\n", optarg);
					exit(-1);					
				}
				break;

			case 'i':
				max_it = atoi(optarg);
				break;
				
			case 'b':
				beta = atof(optarg);
				break;
				
			case 'a':
				val = strtod(optarg, NULL);
				/* Check for various possible errors */
				if ( (errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) ) {
					fprintf(stderr, "Main:: Can't convert alpha: %s\n", optarg);
					exit(EXIT_FAILURE);					
				}
				alpha = val;
	    		break;

			case 'e':
				val = strtod(optarg, NULL);
				/* Check for various possible errors */
				if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
					   || (errno != 0 && val == 0)) {
					fprintf(stderr, "Main:: Can't convert epsilon: %s\n", optarg);
					exit(EXIT_FAILURE);
				}
				epsilon = val;
	    		break;
				
			case 'S':
				strcpy(seq_path, optarg);
				break;


			case 'Y':
				sparsity_type = (Sparsity_Type)atoi(optarg);
				break;
				
			case 'p':
				b_seq = atoi(optarg);
				break;
				
			case 'q':
				b_topo = atoi(optarg);
				break;
				
			case 'I':
				post_iter = atoi(optarg);
				break;
				
			case '?':
	    		print_usage (stdout, 1);
				
			case -1:		/* Done with options. */
			    break;
			
			default:		/* Something else: unexpected. */
				printf("Unknown option %c\n", next_option);
                print_usage (stderr, -1);
		}
    } while (next_option != -1);

	
	char *str = full_path1;
	int len = strlen(str);
	int status = 0;
	register int i, j;
	for(i = len-1; i >= 0; i--) {
		if(str[i] == '.') {
			status = 1;
			break;
		}
		else if(str[i] == '/') {
			status = 2;
			break;
		}		
	}
	if(status == 0)
		strcpy(first, str);
	else if(status == 2)
		strcpy(first, str+i+1);
	else {
		for(j = i; j >=0; j--) {
			if(str[j] == '/') {
				break;
			}	
		}
		strncpy(first, str+j+1, i-j-1);
		
	}


	str = full_path2;
	len = strlen(str);
	status = 0;
	for(i = len-1; i >= 0; i--) {
		if(str[i] == '.') {
			status = 1;
			break;
		}
		else if(str[i] == '/') {
			status = 2;
			break;
		}		
	}
	if(status == 0)
		strcpy(second, str);
	else if(status == 2)
		strcpy(second, str+i+1);
	else {
		for(j = i; j >=0; j--) {
			if(str[j] == '/') {
				break;
			}	
		}
		strncpy(second, str+j+1, i-j-1);		
	}
	
	char *cptr = strrchr(full_path1, '/');	
	strncpy(first, cptr+1, strlen(cptr)-5-(file_type==SMAT?1:0));
	cptr = strrchr(full_path2, '/');
	strncpy(second, cptr+1, strlen(cptr)-5-(file_type==SMAT?1:0));

	printf("Initializng output folder %s ... ", output_path);
	if(!exist(output_path)) {
		char mkcmd[1024];
		sprintf(mkcmd, "mkdir -p %s", output_path);
		int ret = system(mkcmd);
		if(ret) {
			fprintf(stderr, "can't create output dir %s\n", output_path);
			exit(255);
		}
	    //mkdir(output_path, 0700);
	    printf("done\n");
	}
	else {
		printf("already exists\n");
	}
	


	printf("Path to first graph: %s\n"
	"Path to second graph: %s\n"
	"Path to sequence similarity eval file: %s\n"
	"Input file type: %s\n"
	"Output path: %s\n"
	"Initialization type: %s\n"
	"Alpha: %e\n"
	"Beta: %e\n"
	"Number of iterations: %d\n"
	"Epsilon: %e\n"
	"Prior output: %s\n"
	,full_path1, full_path2, seq_path, FileType_names[file_type], output_path, InitType_names[init_type], alpha, beta, max_it, epsilon, prior_output);


	char prefix[1024];
	sprintf(prefix, "%s_vs_%s_%s", first, second, InitType_names[init_type]);


	
	/*******************************************************************************************************************************	
	 *****    Reading input graphs 
	 *******************************************************************************************************************************/
	Graph *G=NULL, *H=NULL; 		
	switch(file_type) {
		case IsoRank:
			if(readISO(full_path1, G) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path1);
				exit(-1);
			}
			
			if(readISO(full_path2, H) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path2);
				exit(-1);
			}
			break;
		
		case Tab:
			if(readTab(full_path1, G) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path1);
				exit(-1);
			}
			
			if(readTab(full_path2, H) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path2);
				exit(-1);
			}
			break;				
			
		case SMAT:
			if(readSMAT(full_path1, G) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path1);
				exit(-1);
			}
			
			if(readSMAT(full_path2, H) == -1) {
				fprintf(stderr, "Error reading network file %s\n", full_path2);
				exit(-1);
			}			
			break;
		default:
			fprintf(stderr, "Main:: Filetype %s is not currently supported.\n", FileType_names[file_type]);
	}
 	writeSMAT(output_path, first, G);
	writeSMAT(output_path, second, H);



	
	/*******************************************************************************************************************************	
	 *****    Initialing input vectors
	 *******************************************************************************************************************************/	
	double *w = NULL;
	if(strcmp(seq_path, "")) {
		char out_path[1024];
		sprintf(out_path, "%s/%s_vs_%s.smat", output_path, first, second);

		if(file_type == SMAT) {
			w = readSeqSim_SMAT(seq_path); // w is not normalized, yet. It will be normalized in ProdTensor::TAME
		}
		else {
			w = readSeqSim(G, H, seq_path, out_path); // w is not normalized, yet. It will be normalized in ProdTensor::TAME
		}
	}
	if( w == NULL) {
		fprintf(stderr, "Main:: Error reading sequence similarity file. Resetting to uniform prior\n");
		long N = (G->n*H->n);
		int mem_size = N*sizeof(double);
		w = (double *)calloc(1, mem_size);
		if(w == NULL) {
			fprintf(stderr, "main:: Error allocating memory for w.\n");
			exit(255);
		}
		double uni_val = 1.0 / N;
		for(register int i = 0; i < N; i++) {
			w[i] = uni_val;
		}
	}
	
	/*******************************************************************************************************************************	
	 *****    Removing "orphan" vertices from the input graphs (nodes with no matching vertices (based on w) in the other graph
	 *******************************************************************************************************************************/
	if(sparsity_type == FullSparsity) {
		register unsigned int i, j;
		int *num_matches_in_H = (int *)calloc(G->n, sizeof(int));
		int *num_matches_in_G = (int *)calloc(H->n, sizeof(int));
		for(i = 0; i < G->n; i++) {
			for(j = 0; j < H->n; j++) {
				if(w[H->n*i + j] != 0) {
					num_matches_in_H[i]++;
					num_matches_in_G[j]++;
				}
			}
		}
		vector <int> G_orphan, H_orphan;

		for(i = 0; i < G->n; i++) {
			if(num_matches_in_H[i] == 0) {
				G_orphan.push_back(i);
				//printf("vertex %d is orphan in G!\n", i+1);
			}
		}
		printf("Number of orphans in G = %d\n", (int)G_orphan.size());
		G->pruneVertices(G_orphan);

		for(i = 0; i < H->n; i++) {
			if(num_matches_in_G[i] == 0) {
				H_orphan.push_back(i);
				//printf("vertex %d is orphan in H!\n", i+1);
			}
		}
		printf("Number of orphans in H = %d\n", (int)H_orphan.size());
		H->pruneVertices(H_orphan);
	}
	/*******************************************************************************************************************************
	 *****    Counting triangles in input graphs 
	 *******************************************************************************************************************************/
	
	printf("Preprocessing triangles ... \n\tCounting ...\n");
	TriangleCounts *T_G = new TriangleCounts(G);
	TriangleCounts *T_H = new TriangleCounts(H);
	printf("\t\t# of Triangles in G: %d\n", T_G->size());
	printf("\t\t# of Triangles in H: %d\n", T_H->size());
	writeSSTEN(output_path, first, G, T_G);
	writeSSTEN(output_path, second, H, T_H);


/*	counter_type **D_H = T_H->get_triEdgeDegs();
	for(int ii = 0; ii < H->n; ii++) {
		for(int jj = 0; jj < H->n; jj++) {
			printf("%d	", D_H[ii][jj]);
		}
		printf("\n");
	}*/
	/*******************************************************************************************************************************	
	 *****    Running TAME ...
	 *******************************************************************************************************************************/

	ProdTensor pt(G, T_G, H, T_H, w, sparsity_type, output_path, prefix); // Initialize knocker product tensor: it will copy pointers  to triangle tensor of each graph
	double *x0 = (double *)calloc(pt.n, sizeof(double));
	if(x0 == NULL) {
		fprintf(stderr, "TAME:: Can't allocate memory for 'x0'\n");
		return -1;
	}
	if(strcmp(prior_output, "")) {
		printf("Reading last saved output from %s\n", prior_output);
		mat X;
		X.load(prior_output, raw_ascii);
		vec X_vec = reshape(X.t(), pt.n, 1);
		memcpy(x0, X_vec.memptr(), pt.n*sizeof(double));
	}
	else {
		pt.InitX(x0, init_type, w);
	}

	pt.issHOPM(max_it, alpha, beta, epsilon, w, x0, init_type, b_seq, b_topo, post_iter);


	delete T_G;
	delete T_H;
	return 0;
}
