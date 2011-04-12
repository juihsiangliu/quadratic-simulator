#include "sparsedoublematrix.h"




SparseDoubleElement * createSparseDoubleElement(const double data)
{
	return createPidSparseDoubleElement(data,0);
}


SparseDoubleElement * createPidSparseDoubleElement(const double data,const int pid)
{
	SparseDoubleElement *ptr = getPidMempoolSet(sizeof(SparseDoubleElement),pid);

	ptr->row = ptr->col = -1;
	ptr->rowLink = ptr->colLink = NULL;
	ptr->data = data;
	
	return ptr;
}




void freeSparseDoubleElement(SparseDoubleElement *ptr)
{
	freePidSparseDoubleElement(ptr,0);
}




void freePidSparseDoubleElement(SparseDoubleElement *ptr,const int pid)
{
	retPidMempoolSet(ptr,sizeof(SparseDoubleElement),pid);
}


//========================================================

SparseDoubleMatrix *createSparseDoubleMatrix(const int row,const int col)
{
	return createPidSparseDoubleMatrix(row,col,0);
}




SparseDoubleMatrix *createPidSparseDoubleMatrix(const int row,const int col,const int pid)
{
	SparseDoubleMatrix *ptr = (SparseDoubleMatrix *)getPidMempoolSet(sizeof(SparseDoubleMatrix),pid);
	ptr->totalRow = row;
	ptr->totalCol = col;
	ptr->nnz = 0;

	ptr->rowIndex = (SparseDoubleElement **)getPidMempoolSet(row*sizeof(SparseDoubleElement *),pid);
	ptr->colIndex = (SparseDoubleElement **)getPidMempoolSet(col*sizeof(SparseDoubleElement *),pid);
	int i;
	for(i=0;i<row;i++)	ptr->rowIndex[i] = createPidSparseDoubleElement(0,pid);
	for(i=0;i<col;i++)	ptr->colIndex[i] = createPidSparseDoubleElement(0,pid);

	return ptr;
}





void freeSparseDoubleMatrix(SparseDoubleMatrix *ptr)
{
	freePidSparseDoubleMatrix(ptr,0);
}




void freePidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseDoubleElement *currentEntry;
	SparseDoubleElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry != NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freePidSparseDoubleElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) freePidSparseDoubleElement(ptr->rowIndex[i],pid);
	for(i=0;i<ptr->totalCol;i++) freePidSparseDoubleElement(ptr->colIndex[i],pid);
	
	retPidMempoolSet(ptr->rowIndex,ptr->totalRow*sizeof(SparseDoubleElement *),pid);
	retPidMempoolSet(ptr->colIndex,ptr->totalCol*sizeof(SparseDoubleElement *),pid);
	retPidMempoolSet(ptr,sizeof(SparseDoubleMatrix),pid);
}



// important: element can NOT be null
// set ptr->data[rowIndex][colIndex] = element
// will copy a new element into the ptr->data[rowIndex][colIndex]
// will automatically create the entry if ptr->data[rowIndex][colIndex] is null && element!=null
void setSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
	setPidSparseDoubleMatrix(ptr,element,rowIndex,colIndex,0);
}




void setPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,const int pid)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseDoubleElement *insert = createPidSparseDoubleElement(element,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data = element;
}






void setFastSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element ,const int rowIndex, const int colIndex,SparseDoubleElement **baseRow,SparseDoubleElement **baseCol)
{
	setFastPidSparseDoubleMatrix(ptr,element ,rowIndex,colIndex,baseRow,baseCol,0);
}






void setFastPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element ,const int rowIndex, const int colIndex,SparseDoubleElement **baseRow,SparseDoubleElement **baseCol,const int pid)
{
	SparseDoubleElement *rowTarget = NULL;
	SparseDoubleElement *colTarget = NULL;
	
	if(baseRow == NULL || baseCol == NULL)
	{
		printf("WTF\n");
	}

	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseDoubleElement *insert = createPidSparseDoubleElement(element,pid);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	
		
		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data = element;

	*baseRow = rowTarget;
	*baseCol = colTarget;
}





double getSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex)
{
	SparseDoubleElement *target = ptr->rowIndex[rowIndex]->rowLink;
	// set target
	while(target != NULL)
	{
		if(target->col >= colIndex ) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=colIndex) return 0;
	else return target->data;
}





double getFastRowSparseDoubleMatrix(const SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseRow)
{
	SparseDoubleElement *target = ptr->rowIndex[rowIndex]->rowLink;
	if(*baseRow!=NULL)
	{
		if( (*baseRow)->row == rowIndex && (*baseRow)->col > target->col && (*baseRow)->col <= colIndex )
		{
			target = *baseRow;
		}
	}

	// set target
	while(target != NULL)
	{
		if(target->col >= colIndex ) break;
		else target = target->rowLink;
	}
	if(target == NULL || target->col!=colIndex)
	{
		return 0;
	}
	else
	{
		*baseRow = target;
		return target->data;
	}
}




double getFastColSparseDoubleMatrix(const SparseDoubleMatrix *ptr,const int rowIndex, const int colIndex,SparseDoubleElement **baseCol)
{
	SparseDoubleElement *target = NULL;

	if(*baseCol!=NULL)
	{
		if((*baseCol)->col == colIndex  && (*baseCol)->row <= rowIndex)
		{
			target = *baseCol;
		}
	}
	if(target == NULL) target = ptr->colIndex[colIndex]->colLink;

	// set target
	while(target != NULL)
	{
		if(target->row >= rowIndex) break;
		else target = target->colLink;
	}

	if(target == NULL || target->row != rowIndex)
	{
		return 0;
	}
	else
	{
		*baseCol = target;
		return target->data;
	}
}




void delSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex)
{
	delPidSparseDoubleMatrix(ptr,rowIndex,colIndex,0);
}




void delPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex,const int colIndex,const int pid)
{
	SparseDoubleElement *rowPrev = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colPrev = ptr->colIndex[colIndex];
	SparseDoubleElement *del;
	// set rowTarget
	while(rowPrev->rowLink != NULL)
	{
		if(rowPrev->rowLink->col >= colIndex ) break;
		else rowPrev = rowPrev->rowLink;
	}
	// set colTarget
	while(colPrev->colLink != NULL)
	{
		if(colPrev->colLink->row >= rowIndex ) break;
		else colPrev = colPrev->colLink;
	}
	// return directly if the entry is null
	if(rowPrev->rowLink ==NULL || rowPrev->rowLink->col!=colIndex) return;
	else 
	{
		del = rowPrev->rowLink;
		// update the information of prev
		rowPrev->rowLink = del->rowLink;
		colPrev->colLink = del->colLink;
		freePidSparseDoubleElement(del,pid);
		ptr->nnz--;
	}
}




void delFastPidSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int rowIndex, const int colIndex, SparseDoubleElement **baseRow, SparseDoubleElement **baseCol,const int pid)
{
	SparseDoubleElement *rowTarget = NULL;
	SparseDoubleElement *colTarget = NULL;
	SparseDoubleElement *del;
	
	// check *baseRow is available and legal
	if((*baseRow)!=NULL)
	{
		if((*baseRow)->row == rowIndex && (*baseRow)->col < colIndex)
		{
			rowTarget = *baseRow;
		}
	}
	// check *baseCol is available and legal
	if((*baseCol)!=NULL)
	{
		if((*baseCol)->col == colIndex && (*baseCol)->row < rowIndex)
		{
			colTarget = *baseCol;
		}
	}
	// not set by rowBase 
	if(rowTarget == NULL)  rowTarget = ptr->rowIndex[rowIndex];
	// not set by colBase 
	if(colTarget == NULL)  colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// return directly if the entry is null
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex) 
	{
		// ...	
	}
	else 
	{
		del = rowTarget->rowLink;
		// update the information of prev
		rowTarget->rowLink = del->rowLink;
		colTarget->colLink = del->colLink;
		freePidSparseDoubleElement(del,pid);
		ptr->nnz--;
	}
	*baseRow = rowTarget;
	*baseCol = colTarget;
}





void dumpSparseDoubleMatrix(FILE *fp,const SparseDoubleMatrix *ptr)
{
	int i;
	fprintf(fp,"totalRow:%d ,totalCol:%d\n",ptr->totalRow,ptr->totalCol);
	fprintf(fp,"nnz= %lld\n",ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *target = ptr->rowIndex[i]->rowLink;
		while(target!=NULL)
		{
//			fprintf(fp,"row: %d, col:%d, ",target->row,target->col);
//			fprintf(fp,"node: %p, rowLink: %p, colLink: %p\n",target,target->rowLink,target->colLink);
//			fprintf(fp,"data = %g\n",target->data);
			target = target->rowLink;
		}
	}
}



void plotSparseDoubleMatrix(FILE *fp, const SparseDoubleMatrix *ptr)
{
	int i,j;
	fprintf(fp,"totalRow:%d ,totalCol:%d\n",ptr->totalRow,ptr->totalCol);
	fprintf(fp,"nnz= %ll\n",ptr->nnz);
	char buf[ptr->totalCol];
	memset(buf,'-',ptr->totalCol);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *target = ptr->rowIndex[i]->rowLink;
		while(target!=NULL)
		{
			buf[target->col] = 'x';
			target = target->rowLink;
		}
		for(j=0;j<ptr->totalCol;j++) fprintf(fp,"%c",buf[j]);
		fprintf(fp,"\n");
		memset(buf,'-',ptr->totalCol);
	}
}





void swapRowSparseDoubleMatrix(SparseDoubleMatrix *ptr, const int row1,const int row2)
{
	SparseDoubleMatrix *tmpMatrix = createSparseDoubleMatrix(ptr->totalRow,ptr->totalCol);

	while(ptr->rowIndex[row1]->rowLink != NULL)
	{
		const int rowNew = row2;
		const int col = ptr->rowIndex[row1]->rowLink->col;
		const double tmp = ptr->rowIndex[row1]->rowLink->data;
		setSparseDoubleMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseDoubleMatrix(ptr,row1,col);
	}
	while(ptr->rowIndex[row2]->rowLink != NULL)
	{
		const int rowNew = row1;
		const int col = ptr->rowIndex[row2]->rowLink->col;
		const double tmp = ptr->rowIndex[row2]->rowLink->data;
		setSparseDoubleMatrix(tmpMatrix,tmp,rowNew,col);
		delSparseDoubleMatrix(ptr,row2,col);
	}
	// dump tempMatrix->row1 to ptr->row1
	while(tmpMatrix->rowIndex[row1]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row1]->rowLink->col;
		const double tmp = tmpMatrix->rowIndex[row1]->rowLink->data;
		setSparseDoubleMatrix(ptr,tmp,row1,col);
		delSparseDoubleMatrix(tmpMatrix,row1,col);
	}
	// dump tempMatrix->row2 to ptr->row2
	while(tmpMatrix->rowIndex[row2]->rowLink != NULL)
	{
		const int col = tmpMatrix->rowIndex[row2]->rowLink->col;
		const double tmp = tmpMatrix->rowIndex[row2]->rowLink->data;
		setSparseDoubleMatrix(ptr,tmp,row2,col);
		delSparseDoubleMatrix(tmpMatrix,row2,col);
	}
	freeSparseDoubleMatrix(tmpMatrix);
}




void clearSparseDoubleMatrix(SparseDoubleMatrix *ptr)
{
	clearPidSparseDoubleMatrix(ptr,0);
}




void clearPidSparseDoubleMatrix(SparseDoubleMatrix *ptr,const int pid)
{
	// free the memory row by row
	int i;
	SparseDoubleElement *currentEntry;
	SparseDoubleElement *removeEntry;
	for(i=0;i<ptr->totalRow;i++)
	{
		currentEntry = ptr->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			removeEntry = currentEntry;
			currentEntry = currentEntry->rowLink;
			freePidSparseDoubleElement(removeEntry,pid);
		}
	}

	for(i=0;i<ptr->totalRow;i++) ptr->rowIndex[i]->rowLink = NULL;
	for(i=0;i<ptr->totalCol;i++) ptr->colIndex[i]->colLink = NULL;
	ptr->nnz = 0;
}




void copySparseDoubleMatrix(SparseDoubleMatrix *dest,const SparseDoubleMatrix *src)
{
	copyPidSparseDoubleMatrix(dest,src,0);
}





void copyPidSparseDoubleMatrix(SparseDoubleMatrix *dest,const SparseDoubleMatrix *src,const int pid)
{
	if(dest == src) return;

	int i;
	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*src->totalRow,pid);
	for(i=0;i<src->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*src->totalCol,pid);
	for(i=0;i<src->totalCol;i++) colCache[i] = NULL;

	clearPidSparseDoubleMatrix(dest,pid);
	
	for(i=0;i<src->totalRow;i++)
	{
		const SparseDoubleElement *currentEntry = src->rowIndex[i]->rowLink;
		while(currentEntry!=NULL)
		{
			const int insRow = currentEntry->row;
			const int insCol = currentEntry->col;
			const double data = currentEntry->data;
//			setPidSparseDoubleMatrix(dest,data,insRow,insCol,pid);
			setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[i],&colCache[insCol],pid);
			currentEntry = currentEntry->rowLink;
		}
	}
	
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*src->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*src->totalCol,pid);
}





void addSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a,const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;

	int i;
	SparseDoubleElement **rowCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalRow);
	for(i=0;i<c->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getMempoolSet(sizeof(SparseDoubleElement *)*c->totalCol);
	for(i=0;i<c->totalCol;i++) colCache[i] = NULL;

	for(i=0;i<a->totalRow;i++)
	{
		// add row by row ,
		// increase min(col1,col2) until row a or row b is null
		// then insert the remaining of that row
		entryA = a->rowIndex[i]->rowLink;
		entryB = b->rowIndex[i]->rowLink;
		while(entryA!=NULL && entryB!=NULL)
		{
			if(entryA->col < entryB->col)
			{
				const double targetA = entryA->data;
//				setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
				setFastSparseDoubleMatrix(cResult,targetA,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const double targetB = entryB->data;
//				setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
				setFastSparseDoubleMatrix(cResult,targetB,i,entryB->col,&rowCache[i],&colCache[entryB->col]);
				entryB = entryB->rowLink;
			}
			else
			{
				const double targetA = entryA->data;
				const double targetB = entryB->data;
//				setSparseDoubleMatrix(cResult,targetA+targetB,i,entryA->col);
				setFastSparseDoubleMatrix(cResult,targetA+targetB,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
//			setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
			setFastSparseDoubleMatrix(cResult,targetA,i,entryA->col,&rowCache[i],&colCache[entryA->col]);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const double targetB = entryB->data;
//			setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
			setFastSparseDoubleMatrix(cResult,targetB,i,entryB->col,&rowCache[i],&colCache[entryB->col]);
			entryB = entryB->rowLink;
		}
	}

	retMempoolSet(rowCache,sizeof(SparseDoubleElement *)*c->totalRow);
	retMempoolSet(colCache,sizeof(SparseDoubleElement *)*c->totalCol);

	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}



// c = a - b
void subSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a,const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;
	int i;
	for(i=0;i<a->totalRow;i++)
	{
		// sub row by row ,
		// increase min(col1,col2) until row a or row b is null
		// then insert the remaining of that row
		entryA = a->rowIndex[i]->rowLink;
		entryB = b->rowIndex[i]->rowLink;
		while(entryA!=NULL && entryB!=NULL)
		{
			if(entryA->col < entryB->col)
			{
				const double targetA = entryA->data;
				setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
				entryA = entryA->rowLink;
			}
			else if(entryA->col > entryB->col)
			{
				const double targetB = -1*entryB->data;
				setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
				entryB = entryB->rowLink;
			}
			else
			{
				const double targetA = entryA->data;
				const double targetB = entryB->data;
				setSparseDoubleMatrix(cResult,targetA-targetB,i,entryA->col);
				entryA = entryA->rowLink;
				entryB = entryB->rowLink;
			}
		}
		// insert the remaining of a
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			setSparseDoubleMatrix(cResult,targetA,i,entryA->col);
			entryA = entryA->rowLink;
		}
		// insert the remaining of b
		while(entryB!=NULL)
		{
			const double targetB = -1*entryB->data;
			setSparseDoubleMatrix(cResult,targetB,i,entryB->col);
			entryB = entryB->rowLink;
		}
	}
	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}



// c = a * b
void mulSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a, const SparseDoubleMatrix *b)
{
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(c->totalRow,c->totalCol);
	SparseDoubleElement *entryA;
	SparseDoubleElement *entryB;
	int i,j;
	for(i=0;i<a->totalRow;i++)
	{
		for(j=0;j<b->totalCol;j++)
		{
			int flag = 0;
			double tempSum = 0;
			entryA = a->rowIndex[i]->rowLink;
			entryB = b->colIndex[j]->colLink;
			// calculate the product of row a(i,:) and col b(:,j) 
			// and then set to c(i,j)
			while(entryA!=NULL && entryB!=NULL)
			{
				if(entryA->col < entryB->row)
				{
					entryA = entryA->rowLink;
				}
				else if(entryA->col > entryB->row)
				{
					entryB = entryB->colLink;
				}
				else
				{
					const double targetA = entryA->data;
					const double targetB = entryB->data;
					tempSum = tempSum + targetA * targetB;
					entryA = entryA->rowLink;
					entryB = entryB->colLink;
					flag = 1;
				}
			}
			if(flag==1) setSparseDoubleMatrix(cResult,tempSum,i,j);
		}
	}
	copySparseDoubleMatrix(c,cResult);
	
	freeSparseDoubleMatrix(cResult);
}


// c = k * a , k is a scale
void scaleSparseDoubleMatrix(SparseDoubleMatrix *c,const double k,const SparseDoubleMatrix *a)
{
	int i;
	SparseDoubleElement *entryA;
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);
	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			setSparseDoubleMatrix(cResult,k*targetA,i,entryA->col);
			entryA = entryA->rowLink;
		}
	}
	copySparseDoubleMatrix(c,cResult);
	freeSparseDoubleMatrix(cResult);
}




// c = trans(A)
void transSparseDoubleMatrix(SparseDoubleMatrix *c, const SparseDoubleMatrix *a)
{
	int i;
	SparseDoubleElement *entryA;
	SparseDoubleMatrix *cResult = createSparseDoubleMatrix(a->totalRow,a->totalCol);

	SparseDoubleElement **rowCache = getMempoolSet(sizeof(SparseDoubleElement *)*cResult->totalRow);
	for(i=0;i<cResult->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getMempoolSet(sizeof(SparseDoubleElement *)*cResult->totalCol);
	for(i=0;i<cResult->totalCol;i++) colCache[i] = NULL;

	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
//			setSparseDoubleMatrix(cResult,targetA,entryA->col,i);
			setFastSparseDoubleMatrix(cResult,targetA,entryA->col,i,&rowCache[entryA->col],&colCache[i]);
			entryA = entryA->rowLink;
		}
	}
	
	retMempoolSet(rowCache,sizeof(SparseDoubleElement *)*cResult->totalRow);
	retMempoolSet(colCache,sizeof(SparseDoubleElement *)*cResult->totalCol);

	copySparseDoubleMatrix(c,cResult);
	freeSparseDoubleMatrix(cResult);
}



static void *blockMulVec(void *ptr)
{
	int i;
	struct ParOfMul2 *par = (struct ParOfMul2 *) ptr;
	const SparseDoubleMatrix *a = par->a;
	const double *b = par->b;
	double *c = par->c;

	SparseDoubleElement *entryA;
	double sum = 0;
	for(i = par->rowBegin ; i<par->rowEnd ; i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		sum = 0;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			const double targetB = b[entryA->col];
			sum = sum + targetA * targetB;
			entryA = entryA->rowLink;
		}
		c[i] = sum;
	}
	pthread_exit(0);
}



void parallelMulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b,const int threadNum)
{
	int i,j;
	memset(c,0,sizeof(double)*a->totalCol);

	const int offset = ceil((double)a->totalRow / (double) threadNum);
	int *start = getMempoolSet(sizeof(int)*(threadNum+1));
	for(i=0;i<threadNum;i++) start[i] = i*offset;
	start[i] = a->totalRow;

	struct ParOfMul2 *parList = getMempoolSet(sizeof(struct ParOfMul2)*threadNum);
	for(i=0;i<threadNum;i++)
	{
		parList[i].rowBegin = start[i];
		parList[i].rowEnd = start[i+1];
		parList[i].a = a;
		parList[i].b = b;
		parList[i].c = c;
	}

	pthread_t pid[threadNum];
	for(i=0;i<threadNum;i++) pthread_create(&pid[i],NULL,blockMulVec,&parList[i]);
	for(i=0;i<threadNum;i++) pthread_join(pid[i],NULL);

	retMempoolSet(parList,sizeof(struct ParOfMul2)*threadNum);
	retMempoolSet(start,sizeof(int)*(threadNum+1));
}





// c = a * b, (sparse multiply dense)
// a is a sparse matrix
// b is a n*1 QuadMatrix
// c is a n*1 QuadMatrix
void mulVecSparseDoubleMatrix(double *c,const SparseDoubleMatrix *a, const double *b)
{
	int i,j;
	SparseDoubleElement *entryA;
	double *cResult = getMempoolSet(sizeof(double)*a->totalCol);
	double sum = 0;

	for(i=0;i<a->totalRow;i++)
	{
		entryA = a->rowIndex[i]->rowLink;
		sum = 0;
		while(entryA!=NULL)
		{
			const double targetA = entryA->data;
			const double targetB = b[entryA->col];
			sum = sum + targetA * targetB;
			entryA = entryA->rowLink;
		}
		cResult[i] = sum;
	}
	memcpy(c,cResult,sizeof(double)*a->totalRow);
	
	retMempoolSet(cResult,sizeof(double)*a->totalCol);
}





// set a as identity matrix
void identitySparseDoubleMatrix(SparseDoubleMatrix *a)
{
	identityPidSparseDoubleMatrix(a,0);
}



void identityPidSparseDoubleMatrix(SparseDoubleMatrix *a,const int pid)
{
	clearPidSparseDoubleMatrix(a,pid);
	int i;
	for(i=0;i<a->totalRow;i++) setPidSparseDoubleMatrix(a,1.0,i,i,pid);
}



// a(i,j) = a(i,j) + element
// will allocate the memory if necessary
void incSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr+element,row,col);
	
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseDoubleElement *insert = createSparseDoubleElement(element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data += element;
}


// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
void decSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];
	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseDoubleElement *insert = createSparseDoubleElement(-1*element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
	}
	else rowTarget->rowLink->data -= element;
	
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr-element,row,col);
}



// a(i,j) = a(i,j) - element
// will allocate the memory if necessary
SparseDoubleElement* decFastSparseDoubleMatrix(SparseDoubleMatrix *ptr,const double element,const int rowIndex, const int colIndex,SparseDoubleElement *baseRow,SparseDoubleElement *baseCol)
{
	SparseDoubleElement *rowTarget = ptr->rowIndex[rowIndex];
	SparseDoubleElement *colTarget = ptr->colIndex[colIndex];

	if(baseRow->col < colIndex) rowTarget = baseRow;
//	if(baseCol!=NULL) colTarget = baseCol;

	if(baseCol!=NULL)
	{
		colTarget = baseCol;
		if(baseCol->colLink!= NULL)
		{
			if(baseCol->colLink->row == rowIndex && baseCol->colLink->col ==  rowTarget->col)
			{
				baseCol->colLink->data -= element;
				return baseCol->colLink;
			}
		}
	}


	// set rowTarget
	while(rowTarget->rowLink != NULL)
	{
		if(rowTarget->rowLink->col >= colIndex ) break;
		else rowTarget = rowTarget->rowLink;
	}
	// set colTarget
	while(colTarget->colLink != NULL)
	{
		if(colTarget->colLink->row >= rowIndex ) break;
		else colTarget = colTarget->colLink;
	}
	// allocate memory if necessary
	if(rowTarget->rowLink == NULL || rowTarget->rowLink->col!=colIndex)
	{
		ptr->nnz++;
		SparseDoubleElement *insert = createSparseDoubleElement(-1*element);
		insert->row = rowIndex;
		insert->col = colIndex;
		insert->rowLink = rowTarget->rowLink;
		insert->colLink = colTarget->colLink;	

		// update the information of target
		rowTarget->rowLink = insert;
		colTarget->colLink = insert;
		return insert;
	}
	else
	{
		rowTarget->rowLink->data -= element;
		return rowTarget->rowLink;
	}
	
//	const double ptr = getSparseDoubleMatrix(a,row,col);
//	setSparseDoubleMatrix(a,ptr-element,row,col);
}





void dense2SparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src)
{
	dense2PidSparseDoubleMatrix(dest,src,0);
}






void dense2PidSparseDoubleMatrix(SparseDoubleMatrix *dest,const double *src,const int pid)
{
	clearPidSparseDoubleMatrix(dest,pid);
	const int row = dest->totalRow;
	const int col = dest->totalCol;
	int i,j,k;
	k = 0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			const double tmp = src[k];
			k++;
			if(tmp!=0)
			{
				setPidSparseDoubleMatrix(dest,tmp,i,j,pid);
			}
		}
	}
}





void sparse2DenseDoubleMatrix(double *dest,const SparseDoubleMatrix *src)
{
	memset(dest,0,sizeof(double)*src->totalRow*src->totalCol);
	int i,j;
	for(i=0;i<src->totalRow;i++)
	{
		SparseDoubleElement *ptr = src->rowIndex[i]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
			const double data = ptr->data;
			setMyMatrix(dest,data,src->totalRow,src->totalCol,i,col);
			ptr = ptr->rowLink;
		}
	}
}




void setDenseCol2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src,const int totalRow, const int targetCol)
{
	int i,j;
	for(i=0;i<totalRow;i++)
	{
		const double element = src[i];
		if( element == 0 ) delSparseDoubleMatrix(dest,i,targetCol);
		else setSparseDoubleMatrix(dest,element,i,targetCol);
	}
}



void setDenseColQuick2SparseDoubleMatrix(SparseDoubleMatrix *dest, const double *src,const int totalRow, const int targetCol)
{
	int i,j;
	for(i=0;i<totalRow;i++)
	{
		const double element = src[i];
		if(element==0) continue;
		else setSparseDoubleMatrix(dest,element,i,targetCol);
	}
}





// the boundaries row1,col1,row2,col2 are included
void clearBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const int row1, const int col1, const int row2, const int col2)
{
	int i;
	// del the old data in the block
	for(i=row1;i<=row2;i++)
	{
		SparseDoubleElement *eachRow = dest->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const int delRow = eachRow->row;
				const int delCol = eachRow->col;
				eachRow = eachRow->rowLink;
				delSparseDoubleMatrix(dest,delRow,delCol);
			}
		}
	}
}




// the boundaries row1,col1,row2,col2 are included
void appendBlockSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int row1, const int col1, const int row2, const int col2)
{
	if(dest == src) return;

	int i;


	// insert new data to the block
	for(i=row1;i<=row2;i++)
	{
		SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow != NULL)
		{
			if(eachRow->col < col1) eachRow = eachRow->rowLink;
			else if(eachRow->col > col2) break;
			else
			{
				const double data = eachRow->data;
				const int insRow = eachRow->row;
				const int insCol = eachRow->col;
				setSparseDoubleMatrix(dest,data,insRow,insCol);
				eachRow = eachRow->rowLink;
			}
		}
	}
}




void mergeSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest)
{
	mergePidSparseDoubleMatrix(dest,src,destRow,destCol,ltRowDest,ltColDest,0);
}





void mergePidSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int destRow, const int destCol,const int ltRowDest, const int ltColDest,const int pid)
{
	int i;
	const int baseRow = ltRowDest;
	const int baseCol = ltColDest;
	
	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	for(i=0;i<dest->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalCol,pid);
	for(i=0;i<dest->totalCol;i++) colCache[i] = NULL;

	for(i=0;i<src->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		int insRow = 0;
		if(eachRow!=NULL) insRow = baseRow + eachRow->row; // must put here ... sometimes the whole row may be empty..
		while(eachRow!=NULL)
		{
			const double data = eachRow->data;	
			const int insCol = baseCol + eachRow->col;
			setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[insRow],&colCache[insCol],pid);
			eachRow = eachRow->rowLink;
		}
	}
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*dest->totalCol,pid);
}




void getSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc)
{
	getPidSubSparseDoubleMatrix(dest,src,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,0);
}





void getPidSubSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *src, const int ltRowSrc, const int ltColSrc, const int rbRowSrc, const int rbColSrc,const int pid)
{
	clearPidSparseDoubleMatrix(dest,pid);
	int i;
	const int baseRow = ltRowSrc;
	const int baseCol = ltColSrc;

	SparseDoubleElement **rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	for(i=0;i<dest->totalRow;i++) rowCache[i] = NULL;
	SparseDoubleElement **colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*dest->totalCol,pid);
	for(i=0;i<dest->totalCol;i++) colCache[i] = NULL;

	for(i=ltRowSrc;i<=rbRowSrc;i++)
	{
		const SparseDoubleElement *eachRow = src->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			if(eachRow->col < ltColSrc) eachRow = eachRow->rowLink;
			else if(eachRow->col > rbColSrc) break;
			else
			{
				const double data = eachRow->data;
				const int insRow = eachRow->row - baseRow;
				const int insCol = eachRow->col - baseCol;
//				setPidSparseDoubleMatrix(dest,data,insRow,insCol,pid);
				setFastPidSparseDoubleMatrix(dest,data,insRow,insCol,&rowCache[insRow],&colCache[insCol],pid);
				eachRow = eachRow->rowLink;
			}
		}
	}
	
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*dest->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*dest->totalCol,pid);
}



// permutate 
void permutateSparseDoubleMatrix(SparseDoubleMatrix *dest, const SparseDoubleMatrix *pRow, const SparseDoubleMatrix *pCol, const SparseDoubleMatrix *src)
{
	SparseDoubleMatrix *destTemp = createSparseDoubleMatrix(dest->totalRow,dest->totalCol);
//	SparseDoubleMatrix *destTemp2 = createSparseDoubleMatrix(dest->totalRow,dest->totalCol);
	int i;

	// permutate row
	SparseDoubleElement *rowCache1 = NULL;
	SparseDoubleElement **colCache1 = getMempoolSet(sizeof(SparseDoubleElement *)*destTemp->totalCol);
	for(i=0;i<destTemp->totalCol;i++) colCache1[i] = NULL;
	for(i=0;i<pRow->totalRow;i++)
	{
		const int colForRowI = pRow->rowIndex[i]->rowLink->col;
		const SparseDoubleElement *ptr = src->rowIndex[colForRowI]->rowLink;
		while(ptr!=NULL)
		{
			const int col = ptr->col;
//			setSparseDoubleMatrix(destTemp,ptr->data,i,col);			
			setFastSparseDoubleMatrix(destTemp,ptr->data,i,col,&rowCache1,&colCache1[col]);			
			ptr = ptr->rowLink;
		}
	}
	retMempoolSet(colCache1,sizeof(SparseDoubleElement *)*destTemp->totalCol);

	clearSparseDoubleMatrix(dest);

	// permutate col
	SparseDoubleElement **rowCache2 = getMempoolSet(sizeof(SparseDoubleElement *)*dest->totalRow);
	SparseDoubleElement *colCache2 = NULL;
	for(i=0;i<dest->totalRow;i++) rowCache2[i] = NULL;
	for(i=0;i<pCol->totalCol;i++)
	{
		const int rowForColI = pCol->colIndex[i]->colLink->row;
		const SparseDoubleElement *ptr = destTemp->colIndex[rowForColI]->colLink;
		while(ptr!=NULL)
		{
			const int row = ptr->row;
//			setSparseDoubleMatrix(destTemp2,ptr->data,row,i);
			setFastSparseDoubleMatrix(dest,ptr->data,row,i,&rowCache2[row],&colCache2);
			ptr = ptr->colLink;
		}
	}
	retMempoolSet(rowCache2,sizeof(SparseDoubleElement *)*dest->totalRow);
	
	freeSparseDoubleMatrix(destTemp);

//	copySparseDoubleMatrix(dest,destTemp2);
//	freeSparseDoubleMatrix(destTemp2);
}




static void updateInvL(SparseDoubleMatrix *invL,const double scale,const int baseRow,const int targetRow)
{
	SparseDoubleElement *baseRowPtr = invL->rowIndex[baseRow]->rowLink;
	while(baseRowPtr!=NULL)
	{
		const double elementSrc = baseRowPtr->data;
		const double elementTarget = getSparseDoubleMatrix(invL,targetRow,baseRowPtr->col);
		const double element = scale*elementSrc + elementTarget;
		setSparseDoubleMatrix(invL,element,targetRow,baseRowPtr->col);
		baseRowPtr = baseRowPtr->rowLink;
	}
}




// inverse lowerTriangular
void invLTSparseDoubleMatrix(SparseDoubleMatrix *invL, const SparseDoubleMatrix *lSrc)
{
	SparseDoubleMatrix *l = createSparseDoubleMatrix(lSrc->totalRow,lSrc->totalCol);
	copySparseDoubleMatrix(l,lSrc);
	identitySparseDoubleMatrix(invL);
	int i;
	for(i = 0; i<l->totalRow;i++)
	{
		SparseDoubleElement *eachRowL = l->rowIndex[i]->rowLink;
		const double base = eachRowL->data;
		eachRowL = eachRowL->colLink;
		while(eachRowL != NULL)
		{
			double element = eachRowL->data;
			double scale = -1.0 * (element/base);
//			fprintf(stderr,"scale:%g,baseRow:%d,targetRow:%d\n",scale,i,eachRowL->row);
			// update invL
			updateInvL(invL,scale,i,eachRowL->row);
//			dumpSparseDoubleMatrix(stderr,invL);
			// updata L
			if(scale!=0)
			{
				SparseDoubleElement *inRow = l->rowIndex[i]->rowLink->rowLink;
				while(inRow != NULL)
				{
					double elementIn = scale * getSparseDoubleMatrix(l,i,inRow->col);
					elementIn = getSparseDoubleMatrix(l,eachRowL->row,inRow->col) + elementIn;
					setSparseDoubleMatrix(l,elementIn,eachRowL->row,inRow->col);
					inRow = inRow->rowLink;
				}
			}
			SparseDoubleElement *next = eachRowL->colLink;
			delSparseDoubleMatrix(l,eachRowL->row,i);
			eachRowL = next;
		}
	}
	freeSparseDoubleMatrix(l);
}





SparseDoubleMatrix *read_csr_pid_SparseDoubleMatrix(const char *filename,const int pid)
{
	// init
	int i;
	FILE *fp = fopen(filename,"rb");
	int *row_ptr = NULL;
	int totalRow;
	int totalCol;
	long long nnz;
	SparseDoubleElement **rowCache = NULL;
	SparseDoubleElement **colCache = NULL;
	// head
	fread(&totalRow,sizeof(int),1,fp);	
	fread(&totalCol,sizeof(int),1,fp);	
	fread(&nnz,sizeof(long long),1,fp);	
	SparseDoubleMatrix *ptr = createPidSparseDoubleMatrix(totalRow,totalCol,pid);
	rowCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*ptr->totalRow,pid);
	for(i=0;i<ptr->totalRow;i++) rowCache[i] = NULL;
	colCache = getPidMempoolSet(sizeof(SparseDoubleElement *)*ptr->totalCol,pid);
	for(i=0;i<ptr->totalCol;i++) colCache[i] = NULL;
	// row_ptr
	row_ptr = getPidMempoolSet(sizeof(int)*(ptr->totalRow+1),pid);
	fread(row_ptr,sizeof(int),ptr->totalRow+1,fp);
	// col
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		int *col_list = getPidMempoolSet(sizeof(int)*length,pid);
		fread(col_list,sizeof(int),length,fp);
		for(j=0;j<length;j++)
		{
			const int row = i;
			const int col = col_list[j];
			setFastPidSparseDoubleMatrix(ptr,0.0,row,col,&rowCache[row],&colCache[col],pid);
		}
		retPidMempoolSet(col_list,sizeof(int)*length,pid);
	}
	// values
	for(i=0;i<ptr->totalRow;i++)
	{
		int j = 0;
		const int length = row_ptr[i+1] - row_ptr[i];
		double *val_list = getPidMempoolSet(sizeof(double)*length,pid);
		fread(val_list,sizeof(double),length,fp);
		SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			eachRow->data = val_list[j];
			j++;
			eachRow = eachRow->rowLink;
		}
		retPidMempoolSet(val_list,sizeof(double)*length,pid);
	}
	// final
	fclose(fp);
	retPidMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1),pid);
	retPidMempoolSet(rowCache,sizeof(SparseDoubleElement *)*ptr->totalRow,pid);
	retPidMempoolSet(colCache,sizeof(SparseDoubleElement *)*ptr->totalCol,pid);
	return ptr;
}





void write_csr_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr)
{
	// init
	int i;
	FILE *fp = fopen(filename,"wb");
	int *row_ptr = getMempoolSet(sizeof(int)*(ptr->totalRow+1));
	// head
	fwrite(&ptr->totalRow,sizeof(int),1,fp);
	fwrite(&ptr->totalCol,sizeof(int),1,fp);
	fwrite(&ptr->nnz,sizeof(long long),1,fp);
	// row_ptr
	row_ptr[0] = 0;
	for(i=0;i<ptr->totalRow;i++)
	{
		row_ptr[i+1] = row_ptr[i];
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			row_ptr[i+1]++;
			eachRow = eachRow->rowLink;
		}
	}
	fwrite(row_ptr,sizeof(int),ptr->totalRow+1,fp);
	retMempoolSet(row_ptr,sizeof(int)*(ptr->totalRow+1));
	// col
	// values	
	long long col_index = 0;
	int *col_list = malloc(sizeof(int)*ptr->nnz);
	double *val_list = malloc(sizeof(double)*ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		const SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			col_list[col_index] = eachRow->col;
			val_list[col_index] = eachRow->data;
			col_index++;
			eachRow = eachRow->rowLink;
		}
	}
	fwrite(col_list,sizeof(int),ptr->nnz,fp);
	fwrite(val_list,sizeof(double),ptr->nnz,fp);
	free(col_list);
	free(val_list);
	// final
	fclose(fp);
}




SparseDoubleMatrix *read_ind_SparseDoubleMatrix(const char *filename)
{
	// init
	long long i;
	FILE *fp = fopen(filename,"r");
	int totalRow;
	int totalCol;
	long long nnz;
	SparseDoubleElement **rowCache = NULL;
	SparseDoubleElement **colCache = NULL;
	// head
	fscanf(fp,"%d %d %lld\n",&totalRow,&totalCol,&nnz);
	SparseDoubleMatrix *ptr = createSparseDoubleMatrix(totalRow,totalCol);

	rowCache = getMempoolSet(sizeof(SparseDoubleElement *)*ptr->totalRow);
	for(i=0;i<ptr->totalRow;i++) rowCache[i] = NULL;
	colCache = getMempoolSet(sizeof(SparseDoubleElement *)*ptr->totalCol);
	for(i=0;i<ptr->totalCol;i++) colCache[i] = NULL;

	int insRow;
	int insCol;
	double val;
	for(i=0;i<nnz;i++)
	{
		fscanf(fp,"%d %d %lf\n",&insRow,&insCol,&val);
		setFastSparseDoubleMatrix(ptr,val,insRow,insCol,&rowCache[insRow],&colCache[insCol]);
	}

	retMempoolSet(rowCache,sizeof(SparseDoubleElement *)*ptr->totalRow);
	retMempoolSet(colCache,sizeof(SparseDoubleElement *)*ptr->totalCol);
	fclose(fp);
	return ptr;
}



void write_ind_SparseDoubleMatrix(const char *filename,const SparseDoubleMatrix *ptr)
{
	int i;
	FILE *fp = fopen(filename,"w");
	fprintf(fp,"%d %d %lld\n",ptr->totalRow,ptr->totalCol,ptr->nnz);
	for(i=0;i<ptr->totalRow;i++)
	{
		SparseDoubleElement *eachRow = ptr->rowIndex[i]->rowLink;
		while(eachRow!=NULL)
		{
			fprintf(fp,"%d %d %lf\n",i,eachRow->col,eachRow->data);
			eachRow = eachRow->rowLink;
		}
	}
	fclose(fp);
}





double colNormSparseDoubleMatrix(const SparseDoubleMatrix *ptr,const int col)
{
	double sum = 0.0;
	const SparseDoubleElement *colPtr = ptr->colIndex[col]->colLink;
	while(colPtr!=NULL)
	{
//		sum += pow(colPtr->data,2);
		sum += colPtr->data * colPtr->data;
		colPtr = colPtr->colLink;
	}
	return sqrt(sum);
}



// ===================================================

CSR_SparseDoubleMatrix *create_CSR_SparseDoubleMatrix(int row, int col, long long nnz)
{
	CSR_SparseDoubleMatrix *mtx = malloc(sizeof(CSR_SparseDoubleMatrix));

	mtx->totalRow = row;
	mtx->totalCol = col;
	mtx->nnz = nnz;

	mtx->rowPtr = malloc(sizeof(int)*(row+1));
	memset(mtx->rowPtr,0,sizeof(int)*(row+1));
	mtx->col = malloc(sizeof(int)*nnz);
	memset(mtx->col,0,sizeof(int)*nnz);
	mtx->val = malloc(sizeof(double)*nnz);
	memset(mtx->val,0,sizeof(double)*nnz);

	return mtx;
}




void free_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	free(mtx->rowPtr);
	free(mtx->col);
	free(mtx->val);
	free(mtx);
}




CSR_SparseDoubleMatrix *read_to_CSR_SparseDoubleMatrix(const char *filename)
{
	int i;
	int totalRow;
	int totalCol;
	long long nnz;
	FILE *fp = fopen(filename,"rb");
	// head
	fread(&totalRow,sizeof(int),1,fp);	
	fread(&totalCol,sizeof(int),1,fp);	
	fread(&nnz,sizeof(long long),1,fp);	
	
	CSR_SparseDoubleMatrix *mtx = create_CSR_SparseDoubleMatrix(totalRow,totalCol,nnz);
	fread(mtx->rowPtr,sizeof(int),totalRow+1,fp);
	fread(mtx->col,sizeof(int),nnz,fp);
	fread(mtx->val,sizeof(double),nnz,fp);

	return mtx;
}




/*
void pseudoWrite_CSR_SparseDoubleMatrix(CSR_SparseDoubleMatrix *mtx)
{
	if(mtx->count == 1) // no others occupy this matrix
	{
		free_CSR_SparseDoubleMatrix(mtx);
	}
	else if(mtx->count > 1) 
	{
		fprintf(stderr,"log: pseudo write counter = %d, skip write to file.\n",mtx->count);
		mtx->count--;
	}
	else
	{
		fprintf(stderr,"error: pseudo write counter = %d, illegal.\n",mtx->count);
	}
}
*/


// ===================================================



/*
// use biconjugate gradient stable method to solve Ax = b
void bcgsSparseQuadMatrix(QuadMatrix *x,const SparseQuadMatrix *A,const QuadMatrix *b, const QuadMatrix *xInit)
{
	const int nodeNum = A->totalRow;
	const int gvNum = A->gvNum;
	const double threshold = 0.001; // if deltax < 0.1%, then leave

	QuadMatrix *x_prev = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *x_next = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *r_init = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *r_prev = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *r_next = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *v_prev = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *v_next = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *p_prev = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *p_next = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *s = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *t = createQuadMatrix(nodeNum,1,gvNum);
	
	QuadElement *rho_prev = createQuadElement(gvNum);
	QuadElement *rho_next = createQuadElement(gvNum);
	QuadElement *alpha = createQuadElement(gvNum);
	QuadElement *beta = createQuadElement(gvNum);
	QuadElement *w_prev = createQuadElement(gvNum);
	QuadElement *w_next = createQuadElement(gvNum);

	setZeroQuadMatrix(x_prev);
	setZeroQuadMatrix(x_next);
	setZeroQuadMatrix(v_next);
	setZeroQuadMatrix(p_next);


	rho_next->m = 1.0;
	alpha->m = 1.0;
	beta->m = 1.0;
	w_next->m = 1.0;
	
	QuadMatrix *rTemp1 = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *rTemp2 = createQuadMatrix(nodeNum,1,gvNum);
	QuadElement *betaTemp1 = createQuadElement(gvNum);
	QuadElement *betaTemp2 = createQuadElement(gvNum);
	QuadMatrix *pTemp1 = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *pTemp2 = createQuadMatrix(nodeNum,1,gvNum);
	QuadElement *alphaTemp = createQuadElement(gvNum);
	QuadMatrix *sTemp = createQuadMatrix(nodeNum,1,gvNum);
	QuadElement *wTemp1 = createQuadElement(gvNum);
	QuadElement *wTemp2 = createQuadElement(gvNum);
	QuadMatrix *xTemp1 = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *xTemp2 = createQuadMatrix(nodeNum,1,gvNum);
	QuadMatrix *rTemp = createQuadMatrix(nodeNum,1,gvNum);


	mulVecSparseQuadMatrix(rTemp1,A,xInit);
	subQuadMatrix(r_init,b,rTemp1);
	copyQuadMatrix(r_next,r_init);
	copyQuadMatrix(x_next,xInit);

	int i;

	for(i=0;i<nodeNum;i++)
	{
		copyQuadMatrix(x_prev,x_next);
		copyQuadMatrix(r_prev,r_next);
		copyQuadMatrix(v_prev,v_next);
		copyQuadMatrix(p_prev,p_next);
		copyQuadElement(w_prev,w_next);
		copyQuadElement(rho_prev,rho_next);

		innerQuadMatrix(rho_next,r_init,r_prev);

		divQuadElement(betaTemp1,rho_next,rho_prev);
		divQuadElement(betaTemp2,alpha,w_prev);
		mulQuadElement(beta,betaTemp1,betaTemp2);

		scaleqQuadMatrix(pTemp1,w_prev,v_prev);
		subQuadMatrix(pTemp2,p_prev,pTemp1);
		scaleqQuadMatrix(pTemp1,beta,pTemp2);
		addQuadMatrix(p_next,r_prev,pTemp1,1);

		mulVecSparseQuadMatrix(v_next,A,p_next);

		innerQuadMatrix(alphaTemp,r_init,v_next);
		divQuadElement(alpha,rho_next,alphaTemp);

		scaleqQuadMatrix(sTemp,alpha,v_next);
		subQuadMatrix(s,r_prev,sTemp);

		mulVecSparseQuadMatrix(t,A,s);

		innerQuadMatrix(wTemp1,t,s);
		innerQuadMatrix(wTemp2,t,t);
		divQuadElement(w_next,wTemp1,wTemp2);

		scaleqQuadMatrix(xTemp1,alpha,p_next);
		scaleqQuadMatrix(xTemp2,w_next,s);
		addQuadMatrix(x_next,x_prev,xTemp1,1);
		addQuadMatrix(x_next,x_next,xTemp2,1);

		if(mRatioQuadMatrix(x_prev,x_next) < threshold) 
		{
			break;
		}
		else
		{
			scaleqQuadMatrix(rTemp,w_next,t);
			subQuadMatrix(r_next,s,rTemp);
		}
	}
	fprintf(stderr,"converge at step: %d\n",i);
	copyQuadMatrix(x,x_next);

	freeQuadMatrix(x_prev);
	freeQuadMatrix(x_next);
	freeQuadMatrix(r_init);
	freeQuadMatrix(r_prev);
	freeQuadMatrix(r_next);
	freeQuadMatrix(v_prev);
	freeQuadMatrix(v_next);
	freeQuadMatrix(p_prev);
	freeQuadMatrix(p_next);
	freeQuadMatrix(s);
	freeQuadMatrix(t);

	freeQuadElement(rho_prev);
	freeQuadElement(rho_next);
	freeQuadElement(alpha);
	freeQuadElement(beta);
	freeQuadElement(w_prev);
	freeQuadElement(w_next);

	freeQuadMatrix(rTemp1);
	freeQuadMatrix(rTemp2);
	freeQuadElement(betaTemp1);
	freeQuadElement(betaTemp2);
	freeQuadMatrix(pTemp1);
	freeQuadMatrix(pTemp2);
	freeQuadElement(alphaTemp);
	freeQuadMatrix(sTemp);
	freeQuadElement(wTemp1);
	freeQuadElement(wTemp2);
	freeQuadMatrix(xTemp1);
	freeQuadMatrix(xTemp2);
	freeQuadMatrix(rTemp);

}


*/



/*
// abs (a-b / a)
static double mRatioQuadMatrix(const QuadMatrix *a, const QuadMatrix *b)
{
	double ret = 0.0;
	int i,j;
	const int row = a->row;
	const int col = a->col;
	double aVal,bVal;

//	dumpQuadMatrix(a);
//	dumpQuadMatrix(b);

//	fprintf(stderr,"before for\n");
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			const QuadElement *ptrA = getPtrEntryQuadMatrix(a,i,j);	
			const QuadElement *ptrB = getPtrEntryQuadMatrix(b,i,j);	
			if(ptrA == NULL) aVal = 0;
			else aVal = ptrA->m;
			if(ptrB == NULL) bVal = 0;
			else bVal = ptrB->m;
			ret = ret +  fabs((aVal-bVal) / aVal);
//			fprintf(stderr,"a = %p b= %p \n",ptrA,ptrB);
//			exit(0);
		}
	}

//	fprintf(stderr,"ret = %lf\n",ret);
	return ret;
}
*/


