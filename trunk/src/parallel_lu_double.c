#include "parallel_lu_double.h"


//static int flagEachNode[16]; // one based
static int *flagEachNode; // one based
// the node[i] is ready to be processed




static pthread_mutex_t donelist_mutex = PTHREAD_MUTEX_INITIALIZER;

// assume the # of leaf of tree node is 8
static ParallelDoneListDouble *createDoneList(ALUDouble **alu, const ParallelETree *tree, SparseDoubleMatrix *a, const int treeInternalSize,const int *postorder)
{
	ParallelDoneListDouble *ptr = getMempoolSet(sizeof(ParallelDoneListDouble));
	ptr->treeInternalSize = treeInternalSize;

	// used to check if all the necessary nodes are done in setNode function
	ptr->orderList = postorder;

	// done
	ptr->done = getMempoolSet(sizeof(int)*(ptr->treeInternalSize+1));
	memset(ptr->done,0,sizeof(int)*(ptr->treeInternalSize+1));

	// eachnodeCurrentDone
	ptr->eachNodeCurrentDone = getMempoolSet(sizeof(int)*(ptr->treeInternalSize+1));
	memset(ptr->eachNodeCurrentDone,0,sizeof(int)*(ptr->treeInternalSize+1));

	ptr->alu = alu;
	ptr->tree = tree;
	ptr->a = a;
	return ptr;
}



static void freeDoneList(ParallelDoneListDouble *ptr)
{
	retMempoolSet(ptr->eachNodeCurrentDone,sizeof(int)*(ptr->treeInternalSize+1));
	retMempoolSet(ptr->done,sizeof(int)*(ptr->treeInternalSize+1));
//	retMempoolSet(ptr->orderList,sizeof(int)*ptr->treeInternalSize);
	retMempoolSet(ptr,sizeof(ParallelDoneListDouble));
}




static int checkNextParallelDoneList(ParallelDoneListDouble *ptr,const int current)
{
	int nextIndex = ptr->orderList[current+1];
	return ptr->done[nextIndex];
}




static void set_1_sort_ParallelDoneList(ParallelDoneListDouble *ptr,ToDoList *todoList,const int index)
{
	pthread_mutex_lock(&donelist_mutex);
	ptr->done[index] = 1;

	if(index%2 == 1) // index is in the right
	{
		if(ptr->done[index-1] == 1) // left is complete
		{
			sortParentsToHeadToDoList(todoList,index);
		}
	}
	else
	{
		sortParentsToHeadToDoList(todoList,index);
	}
	pthread_mutex_unlock(&donelist_mutex);
}





static int isCompleteParallelDoneList(ParallelDoneListDouble *doneList)
{
	pthread_mutex_lock(&donelist_mutex);
	int ret = 1;
	int i;
	for(i=1;i<doneList->treeInternalSize;i++)
	{
		if(doneList->done[i] == 0)
		{
			ret = 0;
			break;
		}
	}
	pthread_mutex_unlock(&donelist_mutex);
	return ret;
}

// =================================================

// N is tree node
// for each tree node N, they have their own baseRow and baseCol
// rowFix = *baseRow + row 
// colFox = *baseCol + col
static void getBaseOfPartialU(const SparseDoubleMatrix *u, const int srcOrder, const ParallelDoneListDouble *ptr, int *baseRow, int *baseCol, const int N,enum OOCFlag oocFlag, struct OOCInfo *oocInfoList)
{
	const int uRow = u->totalRow;
	const int uCol = u->totalCol;

	ALUDouble **aluList = ptr->alu;
	*baseRow = 0;
	*baseCol = 0;
	int i = 0;
	while(ptr->orderList[i]!=N)
	{
		const int order = ptr->orderList[i];

		*baseRow = *baseRow + aluList[order]->uRow;
		i++;
	}
	*baseCol = uCol - aluList[srcOrder]->uCol;
}






// N is tree node
// for each tree node N, they have their own baseRow and baseCol
// rowFix = *baseRow + row 
// colFox = *baseCol + col
static void getBaseOfPartialL(const SparseDoubleMatrix *l, const ParallelDoneListDouble *ptr, int *baseRow, int *baseCol, const int N,enum OOCFlag oocFlag, struct OOCInfo *oocInfoList)
{
	const int lRow = l->totalRow;
	const int lCol = l->totalCol;
	const int treeInternalSize = ptr->treeInternalSize;
	int i;
	ALUDouble **aluList = ptr->alu;
	
	int *eachCol = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=1;i<=treeInternalSize;i++) 
	{
		eachCol[i] = aluList[i]->lCol;
	}

	*baseRow = 0;
	*baseCol = 0;
	i = 0;
	while(ptr->orderList[i]!=N)
	{
		const int order = ptr->orderList[i];
		if(aluList[order]->lRow != aluList[order]->lCol) // the cross node
		{
			*baseCol = *baseCol - eachCol[order*2] - eachCol[order*2+1];
		}
		*baseRow = *baseRow + aluList[order]->lRow;
		*baseCol = *baseCol + aluList[order]->lCol;
		i++;
	}
		
	const int order = ptr->orderList[i];
	
	const SparseDoubleMatrix *partialL = aluList[order]->l;
	if(aluList[order]->lRow != aluList[order]->lCol) // the cross node
	{
		*baseCol = *baseCol - eachCol[order*2] - eachCol[order*2+1];
	}

	retMempoolSet(eachCol,sizeof(int)*(treeInternalSize+1));
}






// =================================================

static ParallelLUDoubleShareData *createParallelLUDoubleShareData(ParallelDoneListDouble *doneList, ToDoList *todolist, const int pid, const int N,const int rootCurrentBegin,const int currentEnd,pthread_mutex_t *mutex, pthread_cond_t *cond,const double *colNormA, const double tol,const int *baseRowL,const int *baseColL, const int *baseRowU, const int *baseColU,struct OOCInfo *oocInfoList,enum OOCFlag oocFlag)
{
	ParallelLUDoubleShareData *ptr = getMempoolSet(sizeof(ParallelLUDoubleShareData));

	ptr->doneList = doneList;
	ptr->todolist = todolist;
	ptr->pid = pid;
	ptr->N = N;
	ptr->rootCurrentBegin = rootCurrentBegin;
	ptr->currentEnd = currentEnd;
	ptr->mutex = mutex;
	ptr->cond = cond;
	ptr->colNormA = colNormA;
	ptr->tol = tol;

	ptr->baseRowL = baseRowL;
	ptr->baseColL = baseColL;
	ptr->baseRowU = baseRowU;
	ptr->baseColU = baseColU;
	ptr->oocInfoList = oocInfoList;
	ptr->oocFlag = oocFlag;

	return ptr;
}




static void freeParallelLUDoubleShareData(ParallelLUDoubleShareData *ptr)
{
	retMempoolSet(ptr,sizeof(ParallelLUDoubleShareData));
}




// =================================================



static ALUDouble **createALUList(const ParallelETree *tree,const int totalRow,const int threadNum,enum OOCFlag oocFlag)
{
	int i;
	ALUDouble **aluList = getMempoolSet(sizeof(ALUDouble *)*tree->size);
	memset(aluList,0,sizeof(ALUDouble *)*tree->size);
	for(i=0;i<tree->size;i++)
	{
		if(tree->node[i] == NULL)
		{
			continue;
		}
		else
		{
			aluList[i] = getPidMempoolSet(sizeof(ALUDouble),i);
			int aRow = tree->node[i]->rowEnd - tree->node[i]->rowBegin + 1; // the row of cross node
			int aCol;
			if(tree->node[i]->type == cross)
			{
				const int leftIndex = 2*i;
				const int rightIndex = 2*i + 1;
				const int aRowNew = aRow + (tree->node[i]->doneRowEnd - tree->node[i]->doneRowBegin +1);
				if(i%2==0) // left node
				{
					int parentIndex = GSL_MAX(1,i/2);
					if( parentIndex == 1) // it is already root
					{
						aCol = totalRow;
					}
					else
					{
						aCol = aluList[parentIndex]->uCol;
					}
				}
				else // right node
				{
					aCol = totalRow - tree->node[i]->doneRowBegin;
				}

				aluList[i]->aRow = aRowNew;
				aluList[i]->aCol = aCol;
				aluList[i]->lRow = aRow;
				aluList[i]->lCol = aRowNew;
				aluList[i]->uRow = aRow;
				aluList[i]->uCol = aCol;
/*
				char name[16];
				aluList[i]->l = createPidSparseDoubleMatrix(aRow,aRowNew,i);
				if(oocFlag == ooc)
				{
					sprintf(name,"lnode%d",i);
					write_csr_SparseDoubleMatrix(name,aluList[i]->l);
					freePidSparseDoubleMatrix(aluList[i]->l,i);
					aluList[i]->l = NULL;
				}
				aluList[i]->u = createPidSparseDoubleMatrix(aRow,aCol,i);
				if(oocFlag == ooc)
				{
					sprintf(name,"unode%d",i);
					write_csr_SparseDoubleMatrix(name,aluList[i]->u);
					freePidSparseDoubleMatrix(aluList[i]->u,i);
					aluList[i]->u = NULL;
				}
*/
			}
			else if(tree->node[i]->type == lu)
			{
				if(i%2 == 0) // left node
				{
					const int rootIndex = i/2;
					aCol = totalRow - tree->node[i]->rowBegin;
				}
				else
				{
					aRow = tree->node[i]->rowEnd - tree->node[i]->rowBegin + 1;
					const int rootIndex = i/2;
					const int leftIndex = i-1;
					const int leftRow = tree->node[leftIndex]->rowEnd - tree->node[leftIndex]->rowBegin + 1;
					const int rootCol = tree->node[rootIndex]->rowEnd - tree->node[rootIndex]->doneRowBegin + 1;
					aCol = totalRow - tree->node[i]->rowBegin;
				}
				
				aluList[i]->aRow = aRow;
				aluList[i]->aCol = aCol;
				aluList[i]->lRow = aRow;
				aluList[i]->lCol = aRow;
				aluList[i]->uRow = aRow;
				aluList[i]->uCol = aCol;
/*
				char name[16];
				aluList[i]->l = createPidSparseDoubleMatrix(aRow,aRow,i);
				if(oocFlag == ooc)
				{
					sprintf(name,"lnode%d",i);
					write_csr_SparseDoubleMatrix(name,aluList[i]->l);
					freePidSparseDoubleMatrix(aluList[i]->l,i);
					aluList[i]->l = NULL;
				}
				aluList[i]->u = createPidSparseDoubleMatrix(aRow,aCol,i);
				if(oocFlag == ooc)
				{
					sprintf(name,"unode%d",i);
					write_csr_SparseDoubleMatrix(name,aluList[i]->u);
					freePidSparseDoubleMatrix(aluList[i]->u,i);
					aluList[i]->u = NULL;
				}
*/
			}
			else
			{
				fprintf(stderr,"in ALU list: undefine\n");
				exit(0);
			}
//			fprintf(stderr,"i:%d, row:%d, col:%d\n",i,aluList[i]->l->totalRow,aluList[i]->l->totalCol);
//			fprintf(stderr,"i:%d, rowBegin:%d, rowEnd:%d, aRow:%d, aCol:%d\n",i,tree->node[i]->rowBegin,tree->node[i]->rowEnd,aRow,aCol);
		}
		clearPidMempoolSet(i);
	}
	return aluList;
}





static void freePidALU(ALUDouble *alu,const int pid)
{
	if(alu!=NULL)
	{
		if(alu->l!=NULL) freePidSparseDoubleMatrix(alu->l,pid);
		if(alu->u!=NULL) freePidSparseDoubleMatrix(alu->u,pid);
		retPidMempoolSet(alu,sizeof(ALUDouble),pid);
		clearPidMempoolSet(pid);
	}
}


// =================================================


static void *threadLU(void *par)
{
	time_t t1,t2;
	ParallelLUDoubleShareData *ptr = (ParallelLUDoubleShareData *)par;
	const int treeNodeID = ptr->N;	
	ALUDouble **alu = ptr->doneList->alu;
	const int pid = ptr->pid;

	flagEachNode[treeNodeID] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);

	alu[treeNodeID]->l = createPidSparseDoubleMatrix(alu[treeNodeID]->lRow,alu[treeNodeID]->lCol,pid);
	alu[treeNodeID]->u = createPidSparseDoubleMatrix(alu[treeNodeID]->uRow,alu[treeNodeID]->uCol,pid);

	if(ptr->oocFlag == ooc)
	{
//		loadByOOCInfo(ptr->oocInfoList,treeNodeID,'L',alu);
//		loadByOOCInfo(ptr->oocInfoList,treeNodeID,'U',alu);
		readByOOCInfo(ptr->oocInfoList,0,0,alu); // a
		ptr->doneList->a = ptr->oocInfoList->a;
	}
	
	fprintf(stderr,"lu begin: %d\n",treeNodeID);
	time(&t1);
	const int ltRowSrc = ptr->doneList->tree->node[treeNodeID]->rowBegin;
	const int ltColSrc = ltRowSrc;
	const int rbRowSrc = ptr->doneList->tree->node[treeNodeID]->rowEnd;
	const int rbColSrc = ltColSrc + alu[treeNodeID]->u->totalCol -1;
	
	SparseDoubleMatrix *a = createPidSparseDoubleMatrix(alu[treeNodeID]->u->totalRow,alu[treeNodeID]->u->totalCol,pid);
	getPidSubSparseDoubleMatrix(a,ptr->doneList->a,ltRowSrc,ltColSrc,rbRowSrc,rbColSrc,pid);
	
	if(ptr->oocFlag == ooc)  pseudoWriteByOOCInfo(ptr->oocInfoList,0,0,alu); // a;

	luPidSparseDoubleMatrix(alu[treeNodeID]->l,alu[treeNodeID]->u,a,pid);
	freePidSparseDoubleMatrix(a,pid);
	time(&t2);
	fprintf(stderr,"lu %d time:%g\n",treeNodeID,difftime(t2,t1));

	clearPidMempoolSet(pid);

	decActiveThread();
	if(ptr->oocFlag == ooc)
	{
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'L',alu);
		writeByOOCInfo(ptr->oocInfoList,treeNodeID,'U',alu);
	}
	set_1_sort_ParallelDoneList(ptr->doneList,ptr->todolist,treeNodeID);
}



//====================================================================


// i and j is the original index of U in setNode (not offset fixed)
static int dropPartialUij(const double Uij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowU, const int *baseColU,const int current)
{
	if(tol == 0) return 0;
	const int fixI = baseRowU[current] + i;
	const int fixJ = baseColU[current] + j;
	if(fixI == fixJ) return 0;

//	printf("current = %d, %d %d\n",current,fixI,fixJ);
	const double norm = colNormA[fixJ];
	if( fabs(Uij) >= tol*norm ) return 0;
	else return 1;
}



// i and j is the original index of L in setNode (not offset fixed)
// fixJ = baseColL[N] + j
//
// fixURowInd = baseRowU[current] + x
// fixUColInd = baseColU[current] + y
//
// to get U(fixJ,fixJ) <=> to solve x and y then seach partialU(x,y)
// fixJ = fixURowInd = fixUColInd
static int dropPartialLij(const double Lij,const double tol,const double *colNormA,const SparseDoubleMatrix *a,const SparseDoubleMatrix *partialU,const int i,const int j,const int *baseRowL, const int *baseColL, const int *baseRowU, const int *baseColU, const int current,const int N)
{
	if(tol == 0) return 0;

	const int fixI = baseRowL[N] + i;
	const int fixJ = baseColL[N] + j;
	if(fixI == fixJ) return 0;

	const double norm = colNormA[fixJ];
	const int x = fixJ - baseRowU[current];
	const int y = fixJ - baseColU[current];
	const double u_fixJ_fixJ = getSparseDoubleMatrix(partialU,x,y);
	
	if( fabs(Lij * u_fixJ_fixJ) >= tol * norm) return 0;
	else return 1;

}





static void *setNode(void *par)
{
	int i,j;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUDoubleShareData *ptr = (ParallelLUDoubleShareData *)par;
	ParallelDoneListDouble *doneList = ptr->doneList;
	ToDoList *todolist = ptr->todolist;

	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
	fprintf(stderr,"setNode begin: %d\n",ptr->N);

	const int N = ptr->N;
	const int rootCurrentBegin = ptr->rootCurrentBegin;
	const int currentEnd = ptr->currentEnd;;
	const int pid = ptr->pid;

	const double tol = ptr->tol;
	const double *colNormA = ptr->colNormA;

	if(ptr->oocFlag == ooc)
	{
//		loadByOOCInfo(ptr->oocInfoList,N,'L',doneList->alu);
//		loadByOOCInfo(ptr->oocInfoList,N,'U',doneList->alu);
		readByOOCInfo(ptr->oocInfoList,0,0,doneList->alu); // a
		doneList->a = ptr->oocInfoList->a;
	}
	
	doneList->alu[N]->l = createPidSparseDoubleMatrix(doneList->alu[N]->lRow,doneList->alu[N]->lCol,pid);
	doneList->alu[N]->u = createPidSparseDoubleMatrix(doneList->alu[N]->uRow,doneList->alu[N]->uCol,pid);

	const int rowLink = doneList->tree->node[N]->rowEnd - doneList->tree->node[N]->rowBegin +1;
	const int colLink = doneList->alu[N]->u->totalCol;
	const int rowCrossBegin = doneList->tree->node[N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[N]->doneRowBegin;
	const int *baseRowL = ptr->baseRowL;
	const int *baseColL = ptr->baseColL;
	const int *baseRowU = ptr->baseRowU;
	const int *baseColU = ptr->baseColU;
	
	SparseDoubleMatrix *l_link = doneList->alu[N]->l;
	SparseDoubleMatrix *u_link = doneList->alu[N]->u;
	getPidSubSparseDoubleMatrix(u_link,doneList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,doneList->a->totalCol-1,pid);

	if(ptr->oocFlag == ooc)  pseudoWriteByOOCInfo(ptr->oocInfoList,0,0,doneList->alu); // a;

	SparseDoubleElement **l_link_row_cache = NULL;
	SparseDoubleElement *l_link_col_cache = NULL;
	SparseDoubleElement *u_link_row_cache = NULL;
	SparseDoubleElement **u_link_col_cache = NULL;
	SparseDoubleElement **del_row_cache = NULL;
	SparseDoubleElement *del_col_cache = NULL;
	
	int idle_counter = 0;
	int rootCurrent = rootCurrentBegin;
	int current = -1;
	int base = 0;
	while(current != currentEnd)
	{
		if(checkNextParallelDoneList(doneList,rootCurrent) == 1)
		{
			l_link_row_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*l_link->totalRow,pid);
			u_link_col_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u_link->totalCol,pid);
			del_row_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u_link->totalRow,pid);
	
			for(i=0;i<l_link->totalRow;i++) l_link_row_cache[i] = NULL;
			l_link_col_cache = NULL;
			u_link_row_cache = NULL;
			for(i=0;i<u_link->totalCol;i++) u_link_col_cache[i] = NULL;
			for(i=0;i<u_link->totalRow;i++) del_row_cache[i] = NULL;
			del_col_cache = NULL;

			idle_counter = 0;
			rootCurrent++;
			current = doneList->orderList[rootCurrent];
			fprintf(stderr,"setNode %d process node %d\n",N,current);
			const int left = current*2;
			const int right = left+1;
			int offset; // how many rows of the cross term or the lu term
	
			int col_offset = 0;
			if(current > doneList->treeInternalSize/2 )  // lu node
			{
				if(ptr->oocFlag == ooc) readByOOCInfo(ptr->oocInfoList,current,'U',doneList->alu);
				
				offset = doneList->alu[current]->uRow;
			}
			else  // internal node
			{
				const int leftRow = doneList->alu[left]->aRow;
				const int rightRow = doneList->alu[right]->aRow;
				offset = doneList->alu[current]->uRow;
				if(ptr->oocFlag == ooc) readByOOCInfo(ptr->oocInfoList,current,'U',doneList->alu);
				
				col_offset = -1 * (leftRow + rightRow);
			}
		
			CSR_SparseDoubleMatrix *partial_u = doneList->alu[current]->csr_u;

			int i,j;
			double scale = 0.0;
			double newLinkU = 0.0;
			double scaledPartialUData = 0.0;
			for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
			{
				SparseDoubleElement *linkUCache1 = NULL;
				for(j=0;j<u_link->totalRow;j++) // sweep u_link row by row
				{
					int partial_u_ptr = partial_u->rowPtr[i];
					scale = getFastColSparseDoubleMatrix(u_link,j,base+i,&linkUCache1)/partial_u->val[partial_u_ptr];
					if(scale == 0.0)
					{
						delFastPidSparseDoubleMatrix(u_link,j,base+i,&del_row_cache[j],&del_col_cache,pid);
						continue;
					}
					setFastPidSparseDoubleMatrix(l_link,scale,j,base+partial_u->col[partial_u_ptr]+col_offset,&l_link_row_cache[j],&l_link_col_cache,pid);
					
					SparseDoubleElement *linkUCache2 = NULL;
					int k;
					for(k=partial_u->rowPtr[i];k<partial_u->rowPtr[i+1];k++)
					{
						const int row = j;
						const int col = base + partial_u->col[k] + col_offset;
						const double data = partial_u->val[k];
						double linkU = getFastRowSparseDoubleMatrix(u_link,row,col,&linkUCache2);
						scaledPartialUData = scale * data;
						newLinkU = linkU - scaledPartialUData;
						setFastPidSparseDoubleMatrix(u_link,newLinkU,row,col,&u_link_row_cache,&u_link_col_cache[col],pid);
					}
					delFastPidSparseDoubleMatrix(u_link,j,base + i,&del_row_cache[j],&del_col_cache,pid);
				}
			}
			base = base + offset;
			doneList->eachNodeCurrentDone[N] = current;

			retPidMempoolSet(l_link_row_cache,sizeof(SparseDoubleElement *)*l_link->totalRow,pid);
			retPidMempoolSet(u_link_col_cache,sizeof(SparseDoubleElement *)*u_link->totalCol,pid);
			retPidMempoolSet(del_row_cache,sizeof(SparseDoubleElement *)*u_link->totalRow,pid);

			if(ptr->oocFlag == ooc) pseudoWriteByOOCInfo(ptr->oocInfoList,current,'U',doneList->alu);
			
		}
		else
		{
			idle_counter++;
			if(idle_counter > 5) sleep(idle_counter);
			if(ptr->oocFlag == ooc)
			{
				fprintf(stderr,"node %d is going to sleep\n",N);
				decActiveThread();
				
				writeByOOCInfo(ptr->oocInfoList,N,'L',doneList->alu);
				writeByOOCInfo(ptr->oocInfoList,N,'U',doneList->alu);

				clearPidMempoolSet(pid);
				flagEachNode[N] = 1;
				pushBackToDoList(todolist,N);
				pthread_cond_wait(ptr->cond,ptr->mutex);
				
				loadByOOCInfo(ptr->oocInfoList,N,'L',doneList->alu);
				loadByOOCInfo(ptr->oocInfoList,N,'U',doneList->alu);
				l_link = doneList->alu[N]->l;
				u_link = doneList->alu[N]->u;

				fprintf(stderr,"node %d is continue to work\n",N);
				time(&t1);
			}
			else // ic
			{
				fprintf(stderr,"node %d is going to sleep\n",N);
				decActiveThread();
				clearPidMempoolSet(pid);
				flagEachNode[N] = 1;
				pushBackToDoList(todolist,N);
				pthread_cond_wait(ptr->cond,ptr->mutex);
				fprintf(stderr,"node %d is continue to work\n",N);
				time(&t1);
			}
		}
	}
	time(&t2);
	fprintf(stderr,"set node %d skew time:%g\n",N,difftime(t2,t1));
	fprintf(stderr,"node %d enter the final stage\n",N);

	l_link_row_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*l_link->totalRow,pid);
	u_link_col_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u_link->totalCol,pid);
	del_row_cache = getPidMempoolSet(sizeof(SparseDoubleElement *)*u_link->totalRow,pid);

	for(i=0;i<l_link->totalRow;i++) l_link_row_cache[i] = NULL;
	l_link_col_cache = NULL;
	u_link_row_cache = NULL;
	for(i=0;i<u_link->totalCol;i++) u_link_col_cache[i] = NULL;
	for(i=0;i<u_link->totalRow;i++) del_row_cache[i] = NULL;
	del_col_cache = NULL;

	double scale = 0.0;
	double newLinkU = 0.0;
	double scaledPartialUData = 0.0;
	for(i=0;i<u_link->totalRow;i++) // for each row in link (to be pivot row)
	{
		SparseDoubleElement *linkUCache1 = NULL;
		for(j=i+1;j<u_link->totalRow;j++) // for each row under the pivot row
		{
			const SparseDoubleElement *pivotRowPtr = u_link->rowIndex[i]->rowLink;
			scale = getFastColSparseDoubleMatrix(u_link,j,base+i,&linkUCache1) / pivotRowPtr->data;
			if(scale == 0.0)
			{
				delFastPidSparseDoubleMatrix(u_link,j,base + i,&del_row_cache[j],&del_col_cache,pid);
				continue;
			}
//			if(!dropPartialLij(scale,tol,colNormA,doneList->a,u_link,j,pivotRowPtr->col,baseRowL,baseColL,baseRowU,baseColU,N,N))
			{
				setFastPidSparseDoubleMatrix(l_link,scale,j,pivotRowPtr->col,&l_link_row_cache[j],&l_link_col_cache,pid);
			}

			SparseDoubleElement *linkUCache2 = NULL;
			while(pivotRowPtr!=NULL) // for each col in pivot row
			{
				const int row = j;
				const int col = pivotRowPtr->col;
				const double data = pivotRowPtr->data;
				double linkU = getFastRowSparseDoubleMatrix(u_link,row,col,&linkUCache2);
				scaledPartialUData = scale * data;
				newLinkU = linkU - scaledPartialUData;
//				if(!dropPartialUij(newLinkU,tol,colNormA,doneList->a,u_link,row,col,baseRowU,baseColU,N))
				{
					setFastPidSparseDoubleMatrix(u_link,newLinkU,row,col,&u_link_row_cache,&u_link_col_cache[col],pid);
				}

				pivotRowPtr = pivotRowPtr->rowLink;
			}
			delFastPidSparseDoubleMatrix(u_link,j,base + i,&del_row_cache[j],&del_col_cache,pid);
		}
	}

	retPidMempoolSet(l_link_row_cache,sizeof(SparseDoubleElement *)*l_link->totalRow,pid);
	retPidMempoolSet(u_link_col_cache,sizeof(SparseDoubleElement *)*u_link->totalCol,pid);
	retPidMempoolSet(del_row_cache,sizeof(SparseDoubleElement *)*u_link->totalRow,pid);
	// ====================================================================
	time(&t3);
	fprintf(stderr,"node %d final time: %g\n",N,difftime(t3,t2));
			
	decActiveThread();
	if(ptr->oocFlag == ooc)
	{
		writeByOOCInfo(ptr->oocInfoList,N,'L',doneList->alu);
		writeByOOCInfo(ptr->oocInfoList,N,'U',doneList->alu);
	}
	clearPidMempoolSet(pid);
	doneList->eachNodeCurrentDone[N] = N;
	set_1_sort_ParallelDoneList(doneList,ptr->todolist,N);
	pthread_exit(0);
}



// =================================================


static void threadHandlerNew(void *par)
{
	struct ThreadHandlerParDouble *thread_handle_ptr = par;

	ParallelLUDoubleShareData **list = thread_handle_ptr->list;
	ParallelDoneListDouble *doneList = list[0]->doneList;
	ToDoList *todolist = list[0]->todolist; 
	const int threadNum = thread_handle_ptr->threadNum;

//	dumpToDoList(stdout,todolist);
	initActiveThread();
	while(!isCompleteParallelDoneList(doneList))
	{
		while(getActiveThread() < threadNum)
		{
			const int indexInToDoList = getFirstToDoList(todolist);
			if(indexInToDoList == -1) // already empty
			{
				if(isCompleteParallelDoneList(doneList))
				{
					break;
				}
				else
				{
					sleep(1);
				}
			}
			else
			{
				safeWaitFlagMatrix(flagEachNode,indexInToDoList);
//				fprintf(stderr,"wake up node %d\n",indexInToDoList);	
				incActiveThread();
//				dumpToDoList(stderr,todolist);
				pthread_cond_signal(list[indexInToDoList-1]->cond);
			}
		}
		usleep(500000);
	}
}





//=====================================================================



static void assembleU(SparseDoubleMatrix *u, ParallelDoneListDouble *ptr, const int orderListSize, const int *baseRow, const int *baseCol,struct OOCInfo* oocInfoList,enum OOCFlag oocFlag)
{
	clearSparseDoubleMatrix(u);
	const int uRow = u->totalRow;
	const int uCol = u->totalCol;
	ALUDouble **aluList = ptr->alu;
	
	int i;	
	for(i=0;i<orderListSize;i++)
	{
		if(oocFlag == ic)
		{
			const int order = ptr->orderList[i];
			const SparseDoubleMatrix *partialU = aluList[order]->u;
			mergeSparseDoubleMatrix(u,partialU,uRow,uCol,baseRow[order],baseCol[order]);
			freePidSparseDoubleMatrix(aluList[order]->u,order);
			aluList[order]->u = NULL;
			clearPidMempoolSet(order);
		}
		else if(oocFlag == ooc)
		{
			const int order = ptr->orderList[i];
			SparseDoubleMatrix *dummyMatrix = createPidSparseDoubleMatrix(u->totalRow,u->totalCol,order);
			loadByOOCInfo(oocInfoList,order,'U',aluList);
			mergePidSparseDoubleMatrix(dummyMatrix,aluList[order]->u,uRow,uCol,baseRow[order],baseCol[order],order);
			freePidSparseDoubleMatrix(aluList[order]->u,order);
			aluList[order]->u = dummyMatrix;
			writeByOOCInfo(oocInfoList,order,'U',aluList);
		}
		else
		{
			fprintf(stderr,"error at assemble U\n");
			exit(0);
		}
	}
}







static void assembleL(SparseDoubleMatrix *l, ParallelDoneListDouble *ptr, const int orderListSize, const int *baseRow, const int *baseCol,struct OOCInfo * oocInfoList,enum OOCFlag oocFlag)
{
	identitySparseDoubleMatrix(l);
	const int lRow = l->totalRow;
	const int lCol = l->totalCol;
	ALUDouble **aluList = ptr->alu;

	int i;
	for(i=0;i<orderListSize;i++)
	{
		if(oocFlag == ic)
		{
			const int order = ptr->orderList[i];
			const SparseDoubleMatrix *partialL = aluList[order]->l;
			mergeSparseDoubleMatrix(l,partialL,lRow,lCol,baseRow[order],baseCol[order]);	
			freePidSparseDoubleMatrix(aluList[order]->l,order);
			aluList[order]->l = NULL;
			clearPidMempoolSet(order);
		}
		else if(oocFlag == ooc)
		{
			const int order = ptr->orderList[i];
			SparseDoubleMatrix *dummyMatrix = createPidSparseDoubleMatrix(l->totalRow,l->totalCol,order);
			identityPidSparseDoubleMatrix(dummyMatrix,order);
			loadByOOCInfo(oocInfoList,order,'L',aluList);
			mergePidSparseDoubleMatrix(dummyMatrix,aluList[order]->l,lRow,lCol,baseRow[order],baseCol[order],order);			
			freePidSparseDoubleMatrix(aluList[order]->l,order);
			aluList[order]->l = dummyMatrix;
			writeByOOCInfo(oocInfoList,order,'L',aluList);
		}
		else
		{
			fprintf(stderr,"error at assemble U\n");
			exit(0);
		}
	}
}



// =================================================


static int findBeforeBegin(const int *postorder_inv,int key,const int treeInternalSize)
{
	if(key > treeInternalSize/2)
	{
		return 0;
	}
	else
	{
		int i = key;
		while(key < treeInternalSize)
		{
			// the left child
			key = key * 2;
		}
		key = key / 2;
		return postorder_inv[key]-1;
	}
}


// =================================================


struct OOCInfo * parallelLUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,enum OOCFlag oocFlag)
{
	return parallelILUDouble(l,u,tree,a,p,threadNum,0.0,oocFlag);
}




struct OOCInfo * parallelILUDouble(SparseDoubleMatrix *l,SparseDoubleMatrix *u, ParallelETree *tree, SparseDoubleMatrix *a, const SparseDoubleMatrix *p, const int threadNum,const double tol,enum OOCFlag oocFlag)
{
	assert(tol==0);
	const int aRow = a->totalRow;
	const int aCol = a->totalCol;

	const int treeInternalSize = (tree->size/2)-1;
	// pre-processing
	int i;
	int status;
	
	// the super nodal tree & its postorder traversal
	int *postorder = getMempoolSet(sizeof(int)*treeInternalSize);
	int *arrayTree = getMempoolSet(sizeof(int)*tree->size);
	memset(arrayTree,0,sizeof(int)*tree->size);
	for(i=1;i<=treeInternalSize;i++) arrayTree[i] = i;
	arrayToPostOrder(postorder,arrayTree,treeInternalSize);
	retMempoolSet(arrayTree,sizeof(int)*tree->size);

	// postorder inverse table
	int *postorder_inv = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=0;i<treeInternalSize;i++)
	{	
		const int index = postorder[i];
		postorder_inv[index] = i;
	}


	// used to communcate the state of each threads
	flagEachNode = getMempoolSet(sizeof(int)*treeInternalSize+1);
	memset(flagEachNode,0,sizeof(int)*treeInternalSize+1);

	getDoneRowInfoNew(tree,postorder,postorder_inv,treeInternalSize);
	
	ALUDouble **aluList = createALUList(tree,aRow,threadNum,oocFlag);
	ParallelDoneListDouble *doneList = createDoneList(aluList,tree,a,treeInternalSize,postorder);

	// used in the ooc mode
	struct OOCInfo *oocInfoList = createOOCInfoList(aluList,treeInternalSize,postorder,oocFlag);
	write_csr_SparseDoubleMatrix(oocInfoList[0].name,a); // a
	freeSparseDoubleMatrix(a);
	clearPidMempoolSet(0);
	
	ToDoList *todolist = createToDoList(treeInternalSize,postorder);

	// set the par entries ... used for the parallel lu
	// mutex and cv list for node
	pthread_mutex_t mutex[treeInternalSize];
	pthread_cond_t cond[treeInternalSize];
	for(i=0;i<treeInternalSize;i++)
	{
		pthread_mutex_init(&mutex[i],NULL);
		pthread_cond_init(&cond[i],NULL);
	}

	// used for ilu
	double *colNormA = NULL;
/*
	if(tol!=0.0)
	{
		colNormA = getMempoolSet(sizeof(double)*a->totalCol);
		for(i=0;i<a->totalCol;i++) 	colNormA[i] = colNormSparseDoubleMatrix(a,i);
	}
*/

	int *baseRowL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseColL = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseRowU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	int *baseColU = getMempoolSet(sizeof(int)*(treeInternalSize+1));
	for(i=0;i<treeInternalSize;i++)
	{
		const int order = postorder[i];
		getBaseOfPartialL(l,doneList,&baseRowL[order],&baseColL[order],order,oocFlag,oocInfoList);
		getBaseOfPartialU(u,order,doneList,&baseRowU[order],&baseColU[order],order,oocFlag,oocInfoList);
	}


	ParallelLUDoubleShareData **parList = getMempoolSet(sizeof(ParallelLUDoubleShareData *)*treeInternalSize);
	
	if(threadNum == 1)
	{
		for(i=0;i<treeInternalSize;i++)
		{
			const int ind = findBeforeBegin(postorder_inv,i+1,treeInternalSize);
			if(i<treeInternalSize/2)
				parList[i] = createParallelLUDoubleShareData(doneList,todolist,1,i+1,ind,2*(i+1)+1,&mutex[i],&cond[i],colNormA,tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag);
			else parList[i] = createParallelLUDoubleShareData(doneList,todolist,1,i+1,0,0,&mutex[i],&cond[i],colNormA,tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag);
		}
	}
	else
	{
		for(i=0;i<treeInternalSize;i++)
		{
			const int ind = findBeforeBegin(postorder_inv,i+1,treeInternalSize);
			if(i<treeInternalSize/2)
				parList[i] = createParallelLUDoubleShareData(doneList,todolist,i+1,i+1,ind,2*(i+1)+1,&mutex[i],&cond[i],colNormA,tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag);
			else parList[i] = createParallelLUDoubleShareData(doneList,todolist,i+1,i+1,0,0,&mutex[i],&cond[i],colNormA,tol,baseRowL,baseColL,baseRowU,baseColU,oocInfoList,oocFlag);
		}
	}

	// send out to each pthread
	pthread_t pid[treeInternalSize];
	for(i=0;i<treeInternalSize;i++)
	{
		if(i<treeInternalSize/2) pthread_create(&pid[i],NULL,setNode,parList[i]);
		else pthread_create(&pid[i],NULL,threadLU,parList[i]);
	}
	sleep(1);

	struct ThreadHandlerParDouble thread_handle;
	thread_handle.list = parList;
	thread_handle.threadNum = threadNum;
	threadHandlerNew(&thread_handle);

	for(i=0;i<treeInternalSize;i++)
	{
		pthread_join(pid[i],NULL);
	}

	// copy to the result
	fprintf(stderr,"copy the result\n");
	fprintf(stderr,"berfore the assemble\n");
	assembleL(l,doneList,treeInternalSize,baseRowL,baseColL,oocInfoList,oocFlag);
	assembleU(u,doneList,treeInternalSize,baseRowU,baseColU,oocInfoList,oocFlag);
	fprintf(stderr,"after the assemble\n");
	
	retMempoolSet(baseRowL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseColL,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseRowU,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(baseColU,sizeof(int)*(treeInternalSize+1));

	if(tol!=0.0)  retMempoolSet(colNormA,sizeof(double)*aCol);
	retMempoolSet(postorder,sizeof(int)*(treeInternalSize));
	retMempoolSet(postorder_inv,sizeof(int)*(treeInternalSize+1));
	retMempoolSet(flagEachNode,sizeof(int)*(treeInternalSize+1));
	for(i=0;i<tree->size;i++) freePidALU(aluList[i],i);
	for(i=0;i<treeInternalSize;i++) freeParallelLUDoubleShareData(parList[i]);
	retMempoolSet(parList,treeInternalSize*sizeof(ParallelLUDoubleShareData*));
	retMempoolSet(aluList,sizeof(ALUDouble *)*tree->size);
	freeDoneList(doneList);
	freeToDoList(todolist);

	return oocInfoList;
}



//=================================

struct OOCInfo *createOOCInfoList(ALUDouble **aluList, const int treeInternalSize, const int *postorder,enum OOCFlag oocFlag)
{
	const int totalLength = 1 + 2*treeInternalSize;
	struct OOCInfo *ptr = getMempoolSet(sizeof(struct OOCInfo)*(totalLength));
	memset(ptr,0,sizeof(struct OOCInfo)*(1+2*treeInternalSize));

	ptr[0].totalLength = totalLength;
	ptr[0].count = 0;
	ptr[0].postorder = getMempoolSet(sizeof(int)*treeInternalSize);
	memcpy(ptr[0].postorder,postorder,sizeof(int)*treeInternalSize);
	sprintf(ptr[0].name,"a.mtx");
	pthread_mutex_init(&ptr[0].mutex, NULL);
	int i;
	
	int rowBegin = 0;
	for(i=0;i<treeInternalSize;i++)
	{
		const int index = postorder[i];
		ptr[index].node = postorder[i];
		ptr[index].type = 'L';
		sprintf(ptr[index].name,"lnode%d",postorder[i]);
		ptr[index].rowBegin = rowBegin;
		rowBegin += aluList[index]->lRow;
		ptr[index].rowEnd = rowBegin;
		ptr[index].count = 0;
		ptr[index].totalLength = totalLength;
		pthread_mutex_init(&ptr[index].mutex, NULL);
	}

	rowBegin = 0;
	for(i=0;i<treeInternalSize;i++)
	{
		const int index = treeInternalSize + postorder[i];
		ptr[index].node = postorder[i];
		ptr[index].type = 'U';
		sprintf(ptr[index].name,"unode%d",postorder[i]);
		ptr[index].rowBegin = rowBegin;
		rowBegin += aluList[postorder[i]]->uRow;
		ptr[index].rowEnd = rowBegin;
		ptr[index].count = 0;
		ptr[index].totalLength = totalLength;
		pthread_mutex_init(&ptr[index].mutex, NULL);
	}

	return ptr;
}


void freeOOCInfoList(struct OOCInfo *ptr)
{
	const int totalLength = ptr[0].totalLength;
	int i;
	for(i=0;i<totalLength;i++) pthread_mutex_destroy(&ptr[i].mutex); // 0 is dummy
	retMempoolSet(ptr[0].postorder,sizeof(int)*(totalLength-1)/2);
	retMempoolSet(ptr,sizeof(struct OOCInfo)*totalLength);
}



void readByOOCInfo(struct OOCInfo *ptr,const int node,const char type, ALUDouble **aluList)
{	
	int index;
	if(node == 0) index = 0;
	else if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0)
	{
		if(index == 0)
		{
			ptr->a = read_csr_pid_SparseDoubleMatrix(ptr[index].name,node);
		}
		else
		{
			if(type == 'L') aluList[node]->l = read_csr_pid_SparseDoubleMatrix(ptr[index].name,node);
			else // U
			{
				aluList[node]->csr_u = read_to_CSR_SparseDoubleMatrix(ptr[index].name); // read as csr format
			}
		}
		fprintf(stderr,"log: read %d from file successfully.\n",node);
	}
	else if(ptr[index].count > 0)
	{
		fprintf(stderr,"log: read %d ... already in mem.\n",node);
	}
	else
	{
		fprintf(stderr,"error: read %d from file failed\n",node);
		exit(0);
	}
	ptr[index].count++;
	pthread_mutex_unlock(&ptr[index].mutex);
}




// logically, this function is used when the matrix is not modified since it is read to memory
void pseudoWriteByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{
	int index;
	if(node == 0) index = 0;
	else if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;
	
	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 1) // no others occupy this file
	{
		if(index == 0)
		{
			freePidSparseDoubleMatrix(ptr->a,node);
			ptr->a = NULL;
		}
		else
		{
			if(type == 'L')
			{
				freePidSparseDoubleMatrix(aluList[node]->l,node);
				aluList[node]->l = NULL;
			}
			else // U
			{
				free_CSR_SparseDoubleMatrix(aluList[node]->csr_u);
			}
		}
		clearPidMempoolSet(node);
		fprintf(stderr,"log: pseudo write %d to file successfully.\n",node);
	}
	else if(ptr[index].count > 1)
	{
		// nothing ... just let it go ...
		fprintf(stderr,"log: pseudo write counter = %d, skip write %d to file.\n",ptr[index].count,node);
	}
	else
	{
		fprintf(stderr,"error: pseudo write %d ... counter illegal....\n",node);
		exit(0);
	}
	ptr[index].count--;
	pthread_mutex_unlock(&ptr[index].mutex);
	
}




void loadByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{	
	int index;
	if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0)
	{
		if(type == 'L') aluList[node]->l = read_csr_pid_SparseDoubleMatrix(ptr[index].name,node);
		else  aluList[node]->u = read_csr_pid_SparseDoubleMatrix(ptr[index].name,node);
		fprintf(stderr,"log: load %d from file successfully.\n",node);
	}
	else
	{
		fprintf(stderr,"error: load %d ... can not load the same file twice....\n",node);
		exit(0);
	}
	pthread_mutex_unlock(&ptr[index].mutex);
}



// the file will be overwritten by this function.
// typicallly, load <=> write , read <=> pseudoWrite
void writeByOOCInfo(struct OOCInfo *ptr,const int node,const char type,ALUDouble **aluList)
{
	int index;
	if(type == 'L') index = node;
	else index = ptr[0].totalLength/2 + node;

	pthread_mutex_lock(&ptr[index].mutex);
	if(ptr[index].count == 0) // no others occupy this file
	{
		if(type == 'L')
		{
			write_csr_SparseDoubleMatrix(ptr[index].name,aluList[node]->l);
			freePidSparseDoubleMatrix(aluList[node]->l,node);
			aluList[node]->l = NULL;
		}
		else
		{
			write_csr_SparseDoubleMatrix(ptr[index].name,aluList[node]->u);
			freePidSparseDoubleMatrix(aluList[node]->u,node);
			aluList[node]->u = NULL;
		}
		clearPidMempoolSet(node);
		fprintf(stderr,"log: write %d to file successfully.\n",node);
	}
/*
	else if(ptr[index].count > 0)
	{
		fprintf(stderr,"log: when writing the matrix %d to file. alread read by other nodes.\n",node);
	}
*/	else // count != 0
	{
		fprintf(stderr,"error: write node %d failed\n",node);
	}
	pthread_mutex_unlock(&ptr[index].mutex);
	
}




//      ax = b
// -> plux = b
// ->  lux = pb  , y = ux 
// ->   ly = pb
// solve y , then solve ux = y
void oocTriSolveSparseDoubleMatrix(double *x, struct OOCInfo *oocInfoList, const SparseDoubleMatrix *p,const SparseDoubleMatrix *pTrans,const double *b)
{
	memset(x,0,sizeof(double)*p->totalRow);
	double *y = getMempoolSet(sizeof(double)*p->totalCol);
	double *pb = getMempoolSet(sizeof(double)*p->totalRow);
	double sum = 0.0;
	double element = 0.0;
	const int treeInternalSize = (oocInfoList[0].totalLength-1)/2;
	const int *postorder = oocInfoList[0].postorder;
	const int base = oocInfoList[0].totalLength/2; // used for u matrix

	int i;
	for(i=0;i<p->totalRow;i++)
	{
		const int col = p->rowIndex[i]->rowLink->col;
		pb[i] = b[col];
	}


	int postorderCounter = 0;
	SparseDoubleElement *eachRow;
	SparseDoubleElement *inRow;
	SparseDoubleMatrix *l = read_csr_pid_SparseDoubleMatrix(oocInfoList[postorder[postorderCounter]].name,0);
//	fprintf(stderr,"read %s done\n",oocInfoList[postorder[postorderCounter]].name);
//	fprintf(stderr,"begin:%d end:%d\n",oocInfoList[postorder[postorderCounter]].rowBegin,oocInfoList[postorder[postorderCounter]].rowEnd);
	// solve ly = pb
	for(i=0;i<l->totalRow;i++)
	{
		if( oocInfoList[postorder[postorderCounter]].rowEnd == i)
		{
			freeSparseDoubleMatrix(l);
			postorderCounter++;	
			l = read_csr_pid_SparseDoubleMatrix(oocInfoList[postorder[postorderCounter]].name,0);
//			fprintf(stderr,"read %s done\n",oocInfoList[postorder[postorderCounter]].name);
//			fprintf(stderr,"begin:%d end:%d\n",oocInfoList[postorder[postorderCounter]].rowBegin,oocInfoList[postorder[postorderCounter]].rowEnd);
		}
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
	freeSparseDoubleMatrix(l);


	postorderCounter = treeInternalSize - 1;
	SparseDoubleMatrix *u = read_csr_pid_SparseDoubleMatrix(oocInfoList[base+postorder[postorderCounter]].name,0);
	double *xTemp = getMempoolSet(sizeof(double)*u->totalRow);
//	fprintf(stderr,"read %d done\n",base + postorder[postorderCounter]);
//	fprintf(stderr,"read %d done\n",postorder[postorderCounter]);
//	fprintf(stderr,"read %s done\n",oocInfoList[base+postorder[postorderCounter]].name);
//	fprintf(stderr,"begin:%d end:%d\n",oocInfoList[base+postorder[postorderCounter]].rowBegin,oocInfoList[base+postorder[postorderCounter]].rowEnd);
	// solve ux = y
	for(i=u->totalRow-1;i>=0;i--)
	{
		if( oocInfoList[base + postorder[postorderCounter]].rowBegin-1 == i)
		{
			freeSparseDoubleMatrix(u);
			postorderCounter--;
			u = read_csr_pid_SparseDoubleMatrix(oocInfoList[base+postorder[postorderCounter]].name,0);
//			fprintf(stderr,"read %s done\n",oocInfoList[base+postorder[postorderCounter]].name);
//			fprintf(stderr,"begin:%d end:%d\n",oocInfoList[base+postorder[postorderCounter]].rowBegin,oocInfoList[base+postorder[postorderCounter]].rowEnd);
		}

		sum = 0;
		eachRow = u->rowIndex[i]->rowLink;
		inRow = eachRow->rowLink;
		while(inRow != NULL)
		{
			const double uij = inRow->data;
			const double xj = xTemp[inRow->col];
			element = uij*xj;
			sum = sum + element;
			inRow = inRow->rowLink;
		}
		element = y[i];
		element = (element-sum) / u->rowIndex[i]->rowLink->data;
		xTemp[i] = element;
	}
	retMempoolSet(y,sizeof(double)*u->totalCol);

//	double *xTemp = getMempoolSet(sizeof(double)*u->totalRow);
	for(i=0;i<u->totalRow;i++)
	{
		int col = pTrans->rowIndex[i]->rowLink->col;
		x[i] = xTemp[col];
	}
//	memcpy(x,xTemp,sizeof(double)*u->totalRow);
	retMempoolSet(xTemp,sizeof(double)*u->totalRow);
	freeSparseDoubleMatrix(u);
}














/*
static void *setNodeDense(void *par)
{
	int i,j,k;
	time_t t1,t2,t3;
	time(&t1);
	ParallelLUDoubleShareData *ptr = (ParallelLUDoubleShareData *)par;
	ParallelDoneListDouble *doneList = ptr->doneList;
	ToDoList *todolist = ptr->todolist;

	flagEachNode[ptr->N] = 1;
	pthread_cond_wait(ptr->cond,ptr->mutex);
	fprintf(stderr,"setNode begin: %d\n",ptr->N);

	const int N = ptr->N;
	const int L = 2*N;
	const int R = L+1;
	const int rootCurrentBegin = ptr->rootCurrentBegin;
	const int currentEnd = ptr->currentEnd;;
	const int pid = ptr->pid;

	const int rowLink = doneList->tree->node[N]->rowEnd - doneList->tree->node[N]->rowBegin +1;
	const int colLink = doneList->alu[N]->u->totalCol;
	const int rowCrossBegin = doneList->tree->node[N]->rowBegin;
	const int rowCrossEnd = doneList->tree->node[N]->rowEnd;
	const int currentDoneRowBegin = doneList->tree->node[N]->doneRowBegin;
	
	SparseDoubleMatrix *l_link = createPidSparseDoubleMatrix(rowLink,doneList->alu[N]->l->totalCol,pid);
	SparseDoubleMatrix *u_link = createPidSparseDoubleMatrix(rowLink,colLink,pid);
	getPidSubSparseDoubleMatrix(u_link,doneList->a,rowCrossBegin,currentDoneRowBegin,rowCrossEnd,doneList->a->totalCol-1,pid);
	
	const int l_link_dense_row = rowLink;
	const int l_link_dense_col = doneList->alu[N]->l->totalCol;
	fprintf(stderr,"allocate dense in %d before, size = %d\n",N,sizeof(double)*l_link_dense_row*l_link_dense_col);
	double *l_link_dense = getPidMempoolSet(sizeof(double)*l_link_dense_row*l_link_dense_col,pid);
	memset(l_link_dense,0,sizeof(double)*l_link_dense_row*l_link_dense_col);
	fprintf(stderr,"allocate dense in %d sucessful\n",N);
	
	const int u_link_dense_row = rowLink;
	const int u_link_dense_col = colLink;
	double *u_link_dense = getPidMempoolSet(sizeof(double)*u_link_dense_row*u_link_dense_col,pid);
	memset(u_link_dense,0,sizeof(double)*u_link_dense_row*u_link_dense_col);

	sparse2DenseDoubleMatrix(l_link_dense,l_link);
	clearPidSparseDoubleMatrix(l_link,pid);
	sparse2DenseDoubleMatrix(u_link_dense,u_link);
	clearPidSparseDoubleMatrix(u_link,pid);

	int rootCurrent = rootCurrentBegin;
	int current = -1;
	int base = 0;
	while(current != currentEnd)
	{
		if(checkNextParallelDoneList(doneList,rootCurrent) == 1)
		{
			rootCurrent++;
			current = doneList->orderList[rootCurrent];
//			printf("setNode %d process node %d\n",N,current);
			const int left = current*2;
			const int right = left+1;
			int offset; // how many rows of the cross term or the lu term
			SparseDoubleMatrix *partial_u;
			if(current > 7)
			{
				offset = doneList->alu[current]->u->totalRow;
				partial_u = doneList->alu[current]->u;
			}
			else
			{
				const int leftRow = doneList->alu[left]->u->totalRow;
				const int rightRow = doneList->alu[right]->u->totalRow;
				offset = doneList->alu[current]->u->totalRow - leftRow - rightRow;
				const int aTotal = doneList->a->totalCol;
				const SparseDoubleMatrix *const uTemp = doneList->alu[current]->u;
				partial_u = createPidSparseDoubleMatrix(uTemp->totalRow-leftRow-rightRow,doneList->a->totalCol-base,pid);
				getPidSubSparseDoubleMatrix(partial_u,uTemp,leftRow+rightRow, leftRow+rightRow,uTemp->totalRow-1,uTemp->totalCol-1,pid);
			}

			int i,j;
			double scale = 0.0;
			double newLinkU = 0.0;
			double scaledPartialUData = 0.0;
			for(i=0;i<partial_u->totalRow;i++) // sweep partial_u row by row
			{
				for(j=0;j<u_link_dense_row;j++) // sweep u_link_dense row by row
				{
					const SparseDoubleElement *partial_u_ptr = partial_u->rowIndex[i]->rowLink;
					assert(partial_u_ptr!=NULL);
					scale = getMyMatrix(u_link_dense,u_link_dense_row,u_link_dense_col,j,base+i) / partial_u_ptr->data;
					if(scale==0)
					{
						setMyMatrix(u_link_dense,0,u_link_dense_row,u_link_dense_col,j,base+i);
						continue;
					}
					setMyMatrix(l_link_dense,scale,l_link_dense_row,l_link_dense_col,j,base+partial_u_ptr->col);
					while(partial_u_ptr!=NULL) // sweep each element in partial row i
					{
						const int row = j;
						const int col = base + partial_u_ptr->col;
						const double data = partial_u_ptr->data;
						const double linkU = getMyMatrix(u_link_dense,u_link_dense_row,u_link_dense_col,row,col);
						scaledPartialUData = scale * data;
						newLinkU = linkU - scaledPartialUData;
						setMyMatrix(u_link_dense,newLinkU,u_link_dense_row,u_link_dense_col,row,col);

						partial_u_ptr = partial_u_ptr->rowLink;
					}
					setMyMatrix(u_link_dense,0,u_link_dense_row,u_link_dense_col,j,base+i);
				}
			}
			
			base = base + offset;
			if(current <=7)
			{
				freePidSparseDoubleMatrix(partial_u,pid);
			}
			doneList->eachNodeCurrentDone[N] = current;
			if(N==1)
			{
				freeUnNecessaryALU(doneList,current);
			}
		}
		else
		{
			fprintf(stderr,"node %d is going to sleep\n",N);
			decActiveThread();
			flagEachNode[N] = 1;
			pushBackToDoList(todolist,N);

			pthread_cond_wait(ptr->cond,ptr->mutex);
			fprintf(stderr,"node %d is continue to work\n",N);
			time(&t1);
		}
	}

	if(N==1) freeExcept1_2and3(doneList);

	time(&t2);
//	printf("set node %d skew time:%g\n",N,difftime(t2,t1));
	fprintf(stderr,"node %d enter the final stage\n",N);

	double scale = 0.0;
	double newLinkU = 0.0;
	double scaledPartialUData = 0.0;
	for(i=0;i<u_link_dense_row;i++) // for each row in link (to be pivot row)
	{
		double *u_link_i = &u_link_dense[i*u_link_dense_col];
		const double pivot = u_link_i[base+i];
		for(j=i+1;j<u_link_dense_row;j++) // for each row under the pivot row
		{
			double *l_link_j = &l_link_dense[j*l_link_dense_col];
			double *u_link_j = &u_link_dense[j*u_link_dense_col];
			const double eachRowHeadData = u_link_j[base+i];
			assert(pivot != 0);
			scale = eachRowHeadData / pivot;
			if(scale == 0)
			{
				u_link_j[base+i] = 0;
				continue;
			}
			l_link_j[base+i] = scale;
			for(k=base;k<u_link_dense_col;k++)
			{
				const double data1 = u_link_i[k];
				const double data2 = u_link_j[k];
				scaledPartialUData = scale *data1;
				newLinkU = data2 - scaledPartialUData;
				u_link_j[k] = newLinkU;
			}
			u_link_j[base+i] = 0;
		}
	}
	time_t begin,end;

	time(&begin);
	dense2PidSparseDoubleMatrix(l_link,l_link_dense,pid);
	dense2PidSparseDoubleMatrix(u_link,u_link_dense,pid);
	retPidMempoolSet(l_link_dense,sizeof(double)*l_link_dense_row*l_link_dense_col,pid);
	retPidMempoolSet(u_link_dense,sizeof(double)*u_link_dense_row*u_link_dense_col,pid);
	time(&end);
//	printf("sparse to dense time: %g\n",difftime(end,begin));

	// ====================================================================
	for(i=0;i<l_link->totalRow;i++)
	{
		const double element = 1.0;
		const int row = i;
		const int col = doneList->alu[L]->l->totalRow + doneList->alu[R]->l->totalRow + i;
		setPidSparseDoubleMatrix(l_link,element,row,col,pid);
	}

	time(&begin);
//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseDoubleMatrix(doneList->alu[N]->l,doneList->alu[L]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,0,0,N);
	mergePidSparseDoubleMatrix(doneList->alu[N]->l,doneList->alu[R]->l,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow,doneList->alu[L]->l->totalRow,N);
	mergePidSparseDoubleMatrix(doneList->alu[N]->l,l_link,doneList->alu[N]->l->totalRow,doneList->alu[N]->l->totalCol,doneList->alu[L]->l->totalRow+doneList->alu[R]->l->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseDoubleMatrix(l_link,pid);

//	pthread_mutex_lock(&setnode_merge_mutex);
	mergePidSparseDoubleMatrix(doneList->alu[N]->u,doneList->alu[L]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,0,0,N);
	mergePidSparseDoubleMatrix(doneList->alu[N]->u,doneList->alu[R]->u,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow,doneList->alu[L]->u->totalRow,N);
	mergePidSparseDoubleMatrix(doneList->alu[N]->u,u_link,doneList->alu[N]->u->totalRow,doneList->alu[N]->u->totalCol,doneList->alu[L]->u->totalRow+doneList->alu[R]->u->totalRow,0,N);
//	pthread_mutex_unlock(&setnode_merge_mutex);
	freePidSparseDoubleMatrix(u_link,pid);
	time(&end);
//	printf("merge time: %g\n",difftime(end,begin));
	
	if(N==1)
	{
		freePidALU(doneList->alu[2],2);
		freePidALU(doneList->alu[3],3);
		clearPidMempoolSet(2);
		clearPidMempoolSet(3);
	}
	// ====================================================================

	time(&t3);
//	printf("node %d final time: %g\n",N,difftime(t3,t2));
			
	set_1_sort_ParallelDoneList(doneList,ptr->todolist,N);
	decActiveThread();

	doneList->eachNodeCurrentDone[N] = N;
	pthread_exit(0);
}
*/



// =================================================




/*
static int isAllAncestorsDone(ParallelDoneListDouble *doneList, const int currentFinishNode)
{
	 int ret = 1;
	 int current = currentFinishNode;
	 while(current>1)
	 {
	 	if(doneList->done[current]!=1)
		{
			ret = 0;
			break;
		}
	 	current = current/2;
	 }
	 return ret;
}



static int indexInOrderList(ParallelDoneListDouble *doneList, const int key)
{
	int i;
	int index = -1;
	for(i=0;i<15;i++)
	{
		if(key == doneList->orderList[i])
		{
			index = i;
		}
	}
	return index;
}




static int isAllAncestorsPartialDone(ParallelDoneListDouble *doneList,const int currentFinishNode)
{
	int ret = 1;
	int current = currentFinishNode;
	int baseIndex = indexInOrderList(doneList,currentFinishNode);
	while(current>0)
	{
		int targetIndex = indexInOrderList(doneList,doneList->eachNodeCurrentDone[current]);
//		printf("current:%d , currentFinishNode:%d\n",current,currentFinishNode);
		if(targetIndex < baseIndex)
		{
			ret = 0;
			break;
		}
		current = current/2;
	}
//	printf("ret=%d\n",ret);
	return ret;
}




//static int freeALU_flag[16] = {0};

static void freeUnNecessaryALU(ParallelDoneListDouble *doneList, const int currentFinishNode)
{
//	if(currentFinishNode > 7) return;

	int i = 0;
	int current = doneList->orderList[i];
	while(1)
	{
		int parent = current/2;
		if(freeALU_flag[current]!=1 && doneList->done[parent]==1 && isAllAncestorsPartialDone(doneList,parent)) // not free yet
		{
//			freePidALU(doneList->alu[current],current);
			clearPidMempoolSet(current);
			freeALU_flag[current] = 1;
			fprintf(stderr,"free alu current:%d\n",current);
		}


		if(current == currentFinishNode)
		{
			break;
		}
		else
		{
			i++;
			current = doneList->orderList[i];	
		}
	}
	
}




static void freeExcept1_2and3(ParallelDoneListDouble *doneList)
{
	int i;
	for(i=4;i<doneList->tree->size;i++)
	{
		if(freeALU_flag[i]!=1)
		{
//			freePidALU(doneList->alu[i],i);
			clearPidMempoolSet(i);
			freeALU_flag[i] = 1;
		}
	}
}



//====================================================================


// get the "full U matrix under (include) node N in the ALU list tree
static SparseDoubleMatrix * getPartialU(const int N, ParallelDoneListDouble *doneList,const int pid)
{
	if(N>7)
	{
		const int row = doneList->alu[N]->u->totalRow;
		const int col = doneList->alu[N]->u->totalCol;
		SparseDoubleMatrix *ret = createPidSparseDoubleMatrix(row,col,pid);
		copyPidSparseDoubleMatrix(ret,doneList->alu[N]->u,pid);
		return ret;
	}
	else
	{
		// recursively get the result from L
		// recursively get the result from R
		// get the result from current
		// merge, free and return
		const int L = 2*N;
		const int R = 2*N + 1;
		
		// L part
		SparseDoubleMatrix *l_sub = getPartialU(L,doneList,pid);
		const int l_sub_row = l_sub->totalRow;
		const int l_sub_col = l_sub->totalCol;
		// R part
		SparseDoubleMatrix *r_sub = getPartialU(R,doneList,pid);
		const int r_sub_row = r_sub->totalRow;
		const int r_sub_col = r_sub->totalCol;
		// N part		
		const int n_sub_row = doneList->alu[N]->u->totalRow;
		const int n_sub_col = doneList->alu[N]->u->totalCol;
		SparseDoubleMatrix *n_sub = doneList->alu[N]->u;
		
		SparseDoubleMatrix *ret = createPidSparseDoubleMatrix(l_sub_row+r_sub_row+n_sub_row,n_sub_col,pid);
		mergePidSparseDoubleMatrix(ret,l_sub,ret->totalRow,ret->totalCol,0,0,pid);
		freePidSparseDoubleMatrix(l_sub,pid);
		mergePidSparseDoubleMatrix(ret,r_sub,ret->totalRow,ret->totalCol,l_sub_row,l_sub_row,pid);
		freePidSparseDoubleMatrix(r_sub,pid);
		mergePidSparseDoubleMatrix(ret,n_sub,ret->totalRow,ret->totalCol,l_sub_row+r_sub_row,0,pid);
		return ret;
	}
}
*/
