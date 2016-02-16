#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <armadillo>


#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_64BIT_WORD


using namespace std;
using namespace arma;



enum FileType {IsoRank, Tab, SMAT};
enum ShiftType {AFFINE_SHIFT, NORM_AFFINE_SHIFT, VARIABLE_SHIFT, COMBINED_SHIFT, MULTI_LINEAR};
enum InitType {Uniform, Random, SeqSim};
enum BinarizationType {TriStat=0, Norm12_bin, Norm24_bin, Z_bin};
enum Sparsity_Type {NoSparsity = 0, PartialSparsity, FullSparsity};
