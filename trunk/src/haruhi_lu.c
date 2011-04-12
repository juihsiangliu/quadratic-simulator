#include "mempool.h"
#include "sparsedoublematrix.h"
#include "solvesparsedoublematrix.h"
#include "parallel_lu_double.h"



double *readB(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	long long nnz = 0;
	long i = 0;
	fscanf(fp,"%lld\n",&nnz);
	double *b = malloc(sizeof(double)*nnz);
	for(i=0;i<nnz;i++)
	{
		fscanf(fp,"%lf\n",&b[i]);
	}
	return b;
}




// argv[1]: a matrix
// argv[2]: b matrix
int main(int argc, char *argv[])
{
	enum OOCFlag oocFlag = ooc;
	time_t t1,t2;
	const int thread = 4;
//	const int goalPartition = 64;
	const int goalPartition = 8;
	// mempool init
	int list[17] = {512,2048,8,8,8,8,8,8,8,8,8,8,8,8,0,0,0};
	createMempoolSet(16,17,1024*1024*1024,list,goalPartition*2); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"

//	int list[27] = {512,512,8,8,8,8,8,8,8,8,8,8,8,8,0,0,0,0,0,0,0,0,0,0,0,0,0};
//	createMempoolSet(16,27,1.5*1024*1024*1024,list,goalPartition*2); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"
	// read a from file
	SparseDoubleMatrix *a = read_ind_SparseDoubleMatrix(argv[1]);
	fprintf(stderr,"loading complete\n");
	// allocate the aRefine and permutation, l, u matrix
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *p = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(a->totalRow,a->totalRow);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	const int aRow = a->totalRow;
	const int aCol = a->totalCol;
	// construct elimination tree
	ParallelETree *tree = createParallelETree(goalPartition*4);
	time(&t1);
	partitionSparseDoubleMatrix(p,pTrans,tree,aRefine,a,goalPartition);
	time(&t2);
	fprintf(stderr,"reorder time:%g\n",difftime(t2,t1));
	freeSparseDoubleMatrix(a);
	clearPidMempoolSet(0);
	// parallel lu
	time(&t1);
	struct OOCInfo *oocInfoList = parallelLUDouble(l,u,tree,aRefine,p,thread,oocFlag);	
	time(&t2);
	fprintf(stderr,"lu time:%g\n",difftime(t2,t1));
	freeParallelETree(tree);
	// tri solve
	double *b = readB(argv[2]);
	double *x = malloc(sizeof(double)*aRow);
	
	// ic
	if(oocFlag == ic)	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	// ooc
	time(&t1);
	if(oocFlag == ooc)  oocTriSolveSparseDoubleMatrix(x,oocInfoList,p,pTrans,b);
	time(&t2);
	fprintf(stderr,"tri solve time:%g\n",difftime(t2,t1));
	freeOOCInfoList(oocInfoList);

	int i;
	fprintf(stdout,"%d\n",aRow);
	for(i=0;i<aRow;i++)fprintf(stdout,"%lf\n",x[i]);
	// free the arrays
	if(oocFlag != ooc) freeSparseDoubleMatrix(aRefine);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(pTrans);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	free(x);
	free(b);

	// finalize the mempool
	FILE *fp_mem = fopen("mempool.log","w");
	usageMempoolSet(fp_mem);
	fclose(fp_mem);
	freeMempoolSet();
}
