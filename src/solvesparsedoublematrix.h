#ifndef SOLVESPARSEDOUBLEMATRIX_H
#define SOLVESPARSEDOUBLEMATRIX_H


#include "mempool.h"
#include "sparsedoublematrix.h"
#include "partition_double.h"
#include <gsl/gsl_vector.h>



// a = lu
void luSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a);
void luPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid);

// return 1 if it can be drop, else 0
int dropLij(const double Lij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *u,const int i,const int j);
int dropUij(const double Uij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *u,const int i,const int j);

// incomplete lu
void iluSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const double tol);
void iluPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid,const double tol);

//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b); 

void triNoPSolveSparseDoubleMatrix(double *x, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b); 

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *A,const double *b);

// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b with pre-defined permutation matrix - p and pTrans
void solveWithPermutationSparseDoubleMatrix(double *x, const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *A,const double *b);


// a = ll^t ... the fake version
void cholSparseDoubleMatrix(SparseDoubleMatrix *l, const SparseDoubleMatrix *a);


// l and u is the ilu factoes
// use ILU
void parallelPCG(const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const SparseDoubleMatrix *a, const double *x_init, const double *b, double *sol,const int threadNum);


#endif
