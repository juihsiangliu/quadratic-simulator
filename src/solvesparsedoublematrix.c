#include "solvesparsedoublematrix.h"



// ============================================




// a = plu
void luSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a)
{
	luPidSparseDoubleMatrix(l,u,a,0);
}




void luPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid)
{
	iluPidSparseDoubleMatrix(l,u,a,pid,0.0);
}



// return 1 if it can be drop, else 0
int dropLij(const double Lij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *u,const int i,const int j)
{
	const double ujj = u->rowIndex[j]->rowLink->data;
	if( fabs(Lij*ujj) >= fabs(tol*colNormA[j]) ) return 0;
//	if( fabs(Lij) >= tol ) return 0;
	else return 1;
}


// return 1 if it can be drop, else 0
int dropUij(const double Uij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *u,const int i,const int j)
{
	if( fabs(Uij) >= fabs(tol*colNormA[j])) return 0;
//	if( fabs(Uij) >= tol) return 0;
	else return 1;
}



void iluSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const double tol)
{
	iluPidSparseDoubleMatrix(l,u,a,0,tol);
}



void iluPidSparseDoubleMatrix(SparseDoubleMatrix *l, SparseDoubleMatrix *u,const SparseDoubleMatrix *a,const int pid,const double tol)
{
	// init & alloc
	double scale = 0.0;
	double element = 0.0;
		
	clearPidSparseDoubleMatrix(l,pid);
	clearPidSparseDoubleMatrix(u,pid);

//	preproceesing
	copyPidSparseDoubleMatrix(u,a,pid);
	identityPidSparseDoubleMatrix(l,pid);
	int i,j;

	SparseDoubleElement *uiCache = NULL;	
	SparseDoubleElement *uEachRowCache = NULL;	

	SparseDoubleElement **lRowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*l->totalRow,pid);
	SparseDoubleElement *lColCache = NULL;
	SparseDoubleElement *uRowCache = NULL;
	SparseDoubleElement **uColCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u->totalCol,pid);
	SparseDoubleElement **delRowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u->totalRow,pid);
	SparseDoubleElement *delColCache = NULL;

	for(i=0;i<l->totalRow;i++) lRowCache[i] = NULL;
	for(i=0;i<u->totalCol;i++) uColCache[i] = NULL;
	for(i=0;i<u->totalRow;i++) delRowCache[i] = NULL;

	double *colNormA = NULL;
	if(tol!=0.0)
	{
		colNormA = getMempoolSet(sizeof(double)*a->totalCol);
		for(i=0;i<a->totalCol;i++) 	colNormA[i] = colNormSparseDoubleMatrix(a,i);
	}

	// perform lu
	for(i=0;i<u->totalRow-1;i++) // permutate
	{
		const SparseDoubleElement *uii = u->rowIndex[i]->rowLink;
		const SparseDoubleElement *eachRow = uii->colLink;
		while(eachRow != NULL) // update E
		{
			scale = eachRow->data / uii->data;
			if(scale!=0)
			{
				if(eachRow->row!=i)
				{
//					if( fabs(scale) >= fabs(tol*(colNormA[i]/uii->data)))
					if(tol==0.0 || !dropLij(scale,tol,colNormA,a,u,eachRow->row,i))
						setFastPidSparseDoubleMatrix(l,scale,eachRow->row,i,&lRowCache[eachRow->row],&lColCache,pid);
				}
				else setFastPidSparseDoubleMatrix(l,scale,eachRow->row,i,&lRowCache[eachRow->row],&lColCache,pid);
				SparseDoubleElement *inRow = uii->rowLink;
				while(inRow != NULL)
				{
					element = scale * getFastRowSparseDoubleMatrix(u,i,inRow->col,&uiCache);
					element = getFastRowSparseDoubleMatrix(u,eachRow->row,inRow->col,&uEachRowCache) - element;
					if(eachRow->row!=inRow->col)
					{
						if(tol==0.0 || !dropUij(element,tol,colNormA,a,u,eachRow->row,inRow->col))
							setFastPidSparseDoubleMatrix(u,element,eachRow->row,inRow->col,&uRowCache,&uColCache[inRow->col],pid);
					}
					else setFastPidSparseDoubleMatrix(u,element,eachRow->row,inRow->col,&uRowCache,&uColCache[inRow->col],pid);
					inRow = inRow->rowLink;
				}
			}
			const SparseDoubleElement *next = eachRow->colLink;
			delFastPidSparseDoubleMatrix(u,eachRow->row,i,&delRowCache[eachRow->row],&delColCache,pid);
			eachRow = next;
		}
	}

	if(tol!=0.0)  retMempoolSet(colNormA,sizeof(double)*a->totalCol);
	retPidMempoolSet(lRowCache,sizeof(SparseDoubleElement *)*l->totalRow,pid);
	retPidMempoolSet(uColCache,sizeof(SparseDoubleElement *)*u->totalCol,pid);
	retPidMempoolSet(delRowCache,sizeof(SparseDoubleElement *)*u->totalRow,pid);
}






static double getSumInChol(const SparseDoubleMatrix *l, const int i,const int j)
{
	double sum = 0.0;
	int colI,colJ;
	const SparseDoubleElement *rowIPtr = l->rowIndex[i]->rowLink;
	const SparseDoubleElement *rowJPtr = l->rowIndex[j]->rowLink;
	while(rowIPtr!=NULL && rowJPtr!=NULL)
	{
		colI = rowIPtr->col;
		colJ = rowJPtr->col;
		if(colI>=j || colJ>=j)
		{
			break;
		}
		else
		{
			if(colI == colJ)
			{
				sum += rowIPtr->data * rowJPtr->data;
				rowIPtr = rowIPtr->rowLink;
				rowJPtr = rowJPtr->rowLink;
			}
			else if(colI > colJ)
			{
				rowJPtr = rowJPtr->rowLink;
			}
			else
			{
				rowIPtr = rowIPtr->rowLink;
			}
		}
	}
	return sum;
}


static void cdiv(SparseDoubleMatrix *l,const int j)
{
	SparseDoubleElement *ajjPtr = l->colIndex[j]->colLink;
	ajjPtr->data = sqrt(ajjPtr->data);
	const double ajj = ajjPtr->data;
	SparseDoubleElement *currentPtr = ajjPtr->colLink;
	while(currentPtr!=NULL)
	{
		int row = currentPtr->row;
		const double aij = currentPtr->data;
		currentPtr->data = aij/ajj;
		currentPtr = currentPtr->colLink;
	}
}




static void cmod(SparseDoubleMatrix *l,const int j,const int k,const double ajk)
{
	int i;
	/*
	for(i=j;i<l->totalRow;i++)
	{
		const double aij = getSparseDoubleMatrix(l,i,j);
		const double aik = getSparseDoubleMatrix(l,i,k);
		const double result = aij-aik*ajk;
		if(result!=0) setSparseDoubleMatrix(l,result,i,j);
		else delSparseDoubleMatrix(l,i,j);
	}
*/
	SparseDoubleElement *colKPtr = l->colIndex[k]->colLink;
	SparseDoubleElement *base = NULL;
	while(colKPtr!=NULL)
	{
		i = colKPtr->row;
		if(i < j)
		{
			colKPtr = colKPtr->colLink;
		}
		else
		{
			const double aik = colKPtr->data;
			const double result = aik*ajk;
			if(result!=0)
			{
//				decSparseDoubleMatrix(l,result,i,j);
				base = decFastSparseDoubleMatrix(l,result,i,j,colKPtr,base);
			}
			else
			{
				fprintf(stderr,"damn\n");
				delSparseDoubleMatrix(l,i,j);
			}
			colKPtr = colKPtr->colLink;
		}
	}
}




void cholSparseDoubleMatrix(SparseDoubleMatrix *l, const SparseDoubleMatrix *a)
{
	// copy the lower triangular from a to l
	clearSparseDoubleMatrix(l);
	int i,j,k;
	for(i=0;i<a->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = a->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			const int row = i;
			const int col = eachRow->col;
			if(col > row)
			{
				break;
			}
			else
			{
				setSparseDoubleMatrix(l,eachRow->data,row,col);		
				eachRow = eachRow->rowLink;
			}
		}
	}

/*
//	dumpSparseDoubleMatrix(stderr,l);
	// perfom sparse cholesky , left looking
	for(j=0;j<a->totalRow;j++) // for each row j
	{
		fprintf(stderr,"%d/%d\n",a->totalRow,j);
		SparseDoubleElement *currentPtr = l->rowIndex[j]->rowLink;
		while(currentPtr!=NULL)
		{
			k = currentPtr->col;
			const double ajk = currentPtr->data;
			if(k>=j) break;
			else
			{
				cmod(l,j,k,ajk);
				currentPtr = currentPtr->rowLink;
			}
		}
		cdiv(l,j);
	}
*/

	// cholesky ~ right looking
	for(k=0;k<a->totalCol;k++)
	{
		fprintf(stderr,"%d/%d\n",k,a->totalRow);
		cdiv(l,k);
		SparseDoubleElement *currentPtr = l->colIndex[k]->colLink;
		while(currentPtr!=NULL)
		{
			j = currentPtr->row;
			const double ajk = currentPtr->data;
			if(j>k)
			{
				cmod(l,j,k,ajk);
			}
			currentPtr = currentPtr->colLink;
		}
	}


//	dumpSparseDoubleMatrix(stderr,l);

	/* the dense implementation
	int i;
	int j;
	clearSparseDoubleMatrix(l);
	for(i=0;i<a->totalRow;i++)
	{
		for(j=0;j<=i;j++)
		{
			if(i==j)
			{
				const double aii = getSparseDoubleMatrix(a,i,i);
				const double sum = getSumInChol(l,i,i);
				const double result = sqrt(aii - sum);
				setSparseDoubleMatrix(l,result,i,i);
			}
			else
			{
				const double ljj = getSparseDoubleMatrix(l,j,j);
				const double aij = getSparseDoubleMatrix(a,i,j);
				const double sum = getSumInChol(l,i,j);
				const double result = (aij-sum)/ljj;
				if(result!=0) setSparseDoubleMatrix(l,result,i,j);
			}
		}
	}
	*/


	
}






//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void triSolveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b)
{
	double *y = getMempoolSet(sizeof(double)*u->totalCol);
	double *pb = getMempoolSet(sizeof(double)*p->totalRow);

	double sum = 0.0;
	double element = 0.0;
	
	int i,k;
	SparseDoubleElement *eachRow;
	SparseDoubleElement *inRow;
//	mulVecSparseDoubleMatrix(pb,p,b);
	for(i=0;i<u->totalRow;i++)
	{
		int col = p->rowIndex[i]->rowLink->col;
		pb[i] = b[col];
	}
	
	// solve ly = pb
	for(i=0;i<l->totalRow;i++)
	{
		sum = 0.0;
		eachRow = l->rowIndex[i]->rowLink;
		inRow = eachRow;
		while(inRow->rowLink != NULL)
		{
			const double lij = inRow->data;
			const double yj = y[inRow->col];
			element = lij*yj;
			sum = sum+element;
			inRow = inRow->rowLink;
		}
		element = pb[i];
		element = (element - sum)/inRow->data;
		y[i] = element;
	}
	retMempoolSet(pb,sizeof(double)*p->totalRow);
	pb = NULL;

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		sum = 0;
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const double uij = inRow->data;
			const double xj = x[inRow->col];
			element = uij*xj;
			sum = sum + element;
			inRow = inRow->rowLink;
		}
		element = y[i];
		element = (element-sum) / u->rowIndex[i]->rowLink->data;
		x[i] = element;
	}
	double *xTemp = getMempoolSet(sizeof(double)*u->totalRow);
	for(i=0;i<u->totalRow;i++)
	{
		int col = pTrans->rowIndex[i]->rowLink->col;
		xTemp[i] = x[col];
	}
	memcpy(x,xTemp,sizeof(double)*u->totalRow);
	retMempoolSet(xTemp,sizeof(double)*u->totalRow);
//	mulVecSparseDoubleMatrix(x,pTrans,x);

	retMempoolSet(y,sizeof(double)*u->totalCol);
}





void triNoPSolveSparseDoubleMatrix(double *x, const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const double *b)
{
	double *y = getMempoolSet(sizeof(double)*u->totalCol);

	double sum = 0.0;
	double element = 0.0;
	
	int i,k;
	SparseDoubleElement *eachRow;
	SparseDoubleElement *inRow;
	
	// solve ly = b
	for(i=0;i<l->totalRow;i++)
	{
		sum = 0.0;
		eachRow = l->rowIndex[i]->rowLink;
		inRow = eachRow;
		while(inRow->rowLink != NULL)
		{
			const double lij = inRow->data;
			const double yj = y[inRow->col];
			element = lij*yj;
			sum = sum+element;
			inRow = inRow->rowLink;
		}
		element = b[i];
		element = (element - sum)/inRow->data;
		y[i] = element;
	}

	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		sum = 0;
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const double uij = inRow->data;
			const double xj = x[inRow->col];
			element = uij*xj;
			sum = sum + element;
			inRow = inRow->rowLink;
		}
		element = y[i];
		element = (element-sum) / u->rowIndex[i]->rowLink->data;
		x[i] = element;
	}
	retMempoolSet(y,sizeof(double)*u->totalCol);
}



// directly use luSparseQuadMatrix() and triSolveSparseQuadMatrix() to solve Ax = b
void solveSparseDoubleMatrix(double *x,const SparseDoubleMatrix *a,const double *b)
{
	const int nodeNum = a->totalRow;
	SparseDoubleMatrix *p = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *pTrans = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *l = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(nodeNum,nodeNum);

	amdSparseDoubleMatrix(p,a);
	transSparseDoubleMatrix(pTrans,p);
	mulSparseDoubleMatrix(aRefine,p,a);
	mulSparseDoubleMatrix(aRefine,aRefine,pTrans);

	luSparseDoubleMatrix(l,u,aRefine);
	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	freeSparseDoubleMatrix(p);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(aRefine);
	freeSparseDoubleMatrix(pTrans);
}




void solveWithPermutationSparseDoubleMatrix(double *x,const SparseDoubleMatrix *p, const SparseDoubleMatrix *pTrans, const SparseDoubleMatrix *a, const double *b)
{
	const int nodeNum = a->totalRow;
	SparseDoubleMatrix *l = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *u = createSparseDoubleMatrix(nodeNum,nodeNum);
	SparseDoubleMatrix *aRefine = createSparseDoubleMatrix(nodeNum,nodeNum);

	mulSparseDoubleMatrix(aRefine,p,a);
	mulSparseDoubleMatrix(aRefine,aRefine,pTrans);

	luSparseDoubleMatrix(l,u,aRefine);
	triSolveSparseDoubleMatrix(x,p,pTrans,l,u,b);
	freeSparseDoubleMatrix(l);
	freeSparseDoubleMatrix(u);
	freeSparseDoubleMatrix(aRefine);
}




static double myNorm2(const double *a,const int size)
{
	int i;
	double sum = 0.0;
	for(i=0;i<size;i++) sum += pow(a[i],2);

	return sqrt(sum);
}


static void dumpVector(FILE *fp, const double *vec, const int size)
{
	int i;
	for(i=0;i<size;i++) fprintf(fp,"%lf ",vec[i]);
	fprintf(fp,"\n");
}



// use ILU
void parallelPCG(const SparseDoubleMatrix *l,const SparseDoubleMatrix *u, const SparseDoubleMatrix *a, const double *x_init, const double *b,double *sol,const int threadNum)
{
	const double drop = 1e-5;
	int i;
	int k = 0;
	const int n = a->totalCol;

	double *rBefore = getMempoolSet(sizeof(double)*n);
	double *rAfter = getMempoolSet(sizeof(double)*n);

	double *xBefore = getMempoolSet(sizeof(double)*n);
	double *xAfter = getMempoolSet(sizeof(double)*n);
	
	double *zBefore = getMempoolSet(sizeof(double)*n);
	double *zAfter = getMempoolSet(sizeof(double)*n);
	
	double *pBefore = getMempoolSet(sizeof(double)*n);
	double *pAfter = getMempoolSet(sizeof(double)*n);

	// r0 = b-Ax0
	double *Ax0 = getMempoolSet(sizeof(double)*n);
	parallelMulVecSparseDoubleMatrix(Ax0,a,x_init,threadNum);
	for(i=0;i<n;i++) rBefore[i] = b[i] - Ax0[i];
	retMempoolSet(Ax0,sizeof(double)*n);
	if(myNorm2(rBefore,2) < drop)
	{
		memcpy(sol,x_init,sizeof(double)*n);
		printf("...\n");
		return;
	}
	// z0 = inv(M) * r0
	triNoPSolveSparseDoubleMatrix(zBefore,l,u,rBefore);
	// p0 = z0
	memcpy(pBefore,zBefore,sizeof(double)*n);

	memcpy(xBefore,x_init,sizeof(double)*n);
	double alpha;
	double beta;
	double *Ap = getMempoolSet(sizeof(double)*n);
	for(k=0;k<n;k++)
	{
		double tmp1;
		double tmp2;
		// Ap
		parallelMulVecSparseDoubleMatrix(Ap,a,pBefore,threadNum);
		// alpha
		tmp1 = tmp2 = 0.0;
		for(i=0;i<n;i++) tmp1 += rBefore[i]*zBefore[i];
		for(i=0;i<n;i++) tmp2 += pBefore[i]*Ap[i];
		alpha = tmp1 / tmp2;
		// x_next = x_now + alpha * p
		for(i=0;i<n;i++) xAfter[i] = xBefore[i] + alpha*pBefore[i];
		// r_next = r_current - alpha * Ap
		for(i=0;i<n;i++) rAfter[i] = rBefore[i] - alpha*Ap[i];	
		// if r_next is small enough , break
		if(myNorm2(rAfter,n) < drop) 
		{
			memcpy(sol,xAfter,sizeof(double)*n);
			break;
		}
		else
		{
			printf("k=%d,%g\n",k,myNorm2(rAfter,n));
			sleep(1);
		}
		// z_next = inv(M) * r_next
		// ...
		triNoPSolveSparseDoubleMatrix(zAfter,l,u,rAfter);
		// beta
		tmp2 = 0.0;
		for(i=0;i<n;i++) tmp2 += rAfter[i]*zAfter[i];
		beta = tmp2 / tmp1;
		// p_next = z_next + beta * p_current
		for(i=0;i<n;i++) pAfter[i] = zAfter[i] + beta * pBefore[i];

		memcpy(rBefore,rAfter,sizeof(double)*n);
		memcpy(xBefore,xAfter,sizeof(double)*n);
		memcpy(zBefore,zAfter,sizeof(double)*n);
		memcpy(pBefore,pAfter,sizeof(double)*n);
	}
	retMempoolSet(Ap,sizeof(double)*n);
	retMempoolSet(rBefore,sizeof(double)*n);
	retMempoolSet(rAfter,sizeof(double)*n);
	retMempoolSet(xBefore,sizeof(double)*n);
	retMempoolSet(xAfter,sizeof(double)*n);
	retMempoolSet(zBefore,sizeof(double)*n);
	retMempoolSet(zAfter,sizeof(double)*n);
	retMempoolSet(pBefore,sizeof(double)*n);
	retMempoolSet(pAfter,sizeof(double)*n);

	fprintf(stderr,"converage in iteration %d\n",k);
}
