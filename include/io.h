#pragma once
#include <generics.h>
#include <graph.h>
#include <triangle.h>
#include <tensor.h>



int readISO(char *path, Graph *&g);
int readTab(char *path, Graph *&g);
int readSMAT(char *path, Graph *&g);
double *readSeqSim(Graph *&g, Graph *&h, char *path, char* out_path);
double *readSeqSim_SMAT(char *path);
int writeSMAT(char *path, char *name, Graph *g);
int writeSSTEN(char *path, char *name, Graph *g, TriangleCounts *t);

void printTriangles(Graph *G, TriangleCounts &T_G);
void printVec(FILE *fd, double *v, int n1, int n2);
int readVec(const char *in_path, double *v, int n1, int n2);
float getTotalSystemMemory();
int exist(char *name);
