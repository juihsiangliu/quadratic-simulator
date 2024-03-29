#include "mempool.h"


static const int _2M = 2*1024*1024;
static const int _4M = 4*1024*1024;
static const int _8M = 8*1024*1024;


static int int_min(const int a,const int b)
{
	if(a<b) return a;
	else return b;
}


static int int_max(const int a,const int b)
{
	if(a<b) return b;
	else return a;
}

//==================================================




static void dumpMempool(FILE *fp,const Mempool *pool)
{
	fprintf(fp,"===============\n");
	fprintf(fp,"pool name: %s\n",pool->poolName);
	fprintf(fp,"batchNumOfBlock: %d\n",pool->batchNumOfBlock);
	fprintf(fp,"blockSize: %d\n",pool->blockSize);

	fprintf(fp,"*** freeList ***\n");
	dumpDqueue(fp,pool->dqueue);
}



static Mempool *createMempool(const char *poolName,const int blockSize,const int numOfBlock)
{
	int i;
	Mempool *pool = (Mempool *)malloc(sizeof(Mempool));
	strcpy(pool->poolName,poolName);
	pool->blockSize = blockSize;


	const int maxHeapSize = 2*int_max(numOfBlock,1);

	pool->dqueue = createDqueue(maxHeapSize);
	pool->batchNumOfBlock = int_max(numOfBlock,1);

/*
	for(i=0;i< numOfBlock ;i++)
	{
		void *ptr = malloc(blockSize);
		if(ptr == NULL) fprintf(stderr,"error: not enought memory\n");
		insertTailDqueue(pool->dqueue,ptr);
	}
*/

	pool->count = 0;
	pool->malloc_count = 0;
	pthread_mutex_init(&(pool->mutex), NULL);


	return pool;
}




static void freeMempool(Mempool *pool)
{	
	int i;
	Dqueue *dqueue = pool->dqueue;
	while(!isEmptyDqueue(dqueue))
	{
		void *mem = delTailDqueue(dqueue);
		free(mem);
		pool->malloc_count--;
	}
	freeDqueue(dqueue);
	pthread_mutex_destroy(&(pool->mutex));
	free(pool);
}




static void clearMempool(Mempool *pool)
{
	pthread_mutex_lock(&(pool->mutex));
	Dqueue *dqueue = pool->dqueue;
	while(!isEmptyDqueue(dqueue))
	{
		void *mem = delTailDqueue(dqueue);
		free(mem);
		pool->malloc_count--;
	}
	pthread_mutex_unlock(&(pool->mutex));
}




// return a mem block from pool to user 
static void *getMempool(Mempool *pool)
{
	pthread_mutex_lock(&(pool->mutex));

	if(isEmptyDqueue(pool->dqueue))
	{
		int i;
		for(i=0;i<pool->batchNumOfBlock;i++)
		{
			if(!isFullDqueue(pool->dqueue))
			{
				void *ptr = malloc(pool->blockSize);
				pool->malloc_count++;
				if(ptr == NULL) fprintf(stderr,"error: not enought memory\n");
				insertTailDqueue(pool->dqueue,ptr);
			}
			else break;
		}
	}
	void *ret = delTailDqueue(pool->dqueue);
	pool->count++;

	pthread_mutex_unlock(&(pool->mutex));

	return ret;
}



static void retMempool(Mempool *pool, void *data)
{
	pthread_mutex_lock(&pool->mutex);
	if(isFullDqueue(pool->dqueue))
	{
		void *head = delHeadDqueue(pool->dqueue);
		free(head);
		pool->malloc_count--;
	}
	insertTailDqueue(pool->dqueue,data);
	pool->count--;
	pthread_mutex_unlock(&(pool->mutex));
}


// ======================================================



static MempoolSet *createMempoolSetKernelNew(const int minBlockSize,const int numOfMempool,const int hugeBlockSize,const int *const numOfBlockList)
{
	MempoolSet *set = (MempoolSet *)malloc(sizeof(MempoolSet));
	set->poolSet = (Mempool **)malloc(sizeof(Mempool *)*numOfMempool);
	int i;
	set->minBlockSize = minBlockSize;
	set->maxBlockSize = minBlockSize;
	set->hugeBlockSize = hugeBlockSize;
	set->numOfPool = numOfMempool;
	set->lgMinBlockSize = log2(set->minBlockSize);

	for(i=0;i<numOfMempool;i++)
	{
		char name[32];
		int n;
		n = sprintf(name,"normal:%d",set->maxBlockSize);
		set->poolSet[i] = createMempool(name,set->maxBlockSize,numOfBlockList[i]);
		set->maxBlockSize = set->maxBlockSize*2;
	}
	set->maxBlockSize = set->maxBlockSize/2;

	set->huge = createMempool("huge",set->hugeBlockSize,0);


	return set;
}





static void freeMempoolSetKernel(MempoolSet *ptr)
{
	int i;
	for(i=0;i<ptr->numOfPool;i++) freeMempool(ptr->poolSet[i]);
	free(ptr->poolSet);
	freeMempool(ptr->huge);
	free(ptr);
}



static void clearMempoolSetKernel(MempoolSet *ptr)
{
	int i;
	for(i=0;i<ptr->numOfPool;i++) clearMempool(ptr->poolSet[i]);
	clearMempool(ptr->huge);
}





static void *getMempoolSetKernel(MempoolSet *ptr,const int blockSize)
{
	void *ret = NULL;
	if(blockSize > ptr->maxBlockSize && blockSize <= ptr->hugeBlockSize)
	{
//		ret = getMempool(ptr->huge);
		ret = malloc(blockSize);
		ptr->huge->count++;
	}
	else if(blockSize <= ptr->maxBlockSize)
	{
		const double lgBlockSize = log2(blockSize);
		int index = ceil(lgBlockSize - ptr->lgMinBlockSize);
		if(index < 0) index = 0;
		ret = getMempool(ptr->poolSet[index]);
	}
	else
	{
		fprintf(stderr,"Error in get MempoolSet (Huge) \n");
		fprintf(stderr,"request size:%d , hugeSize:%d\n",blockSize,ptr->hugeBlockSize);
		exit(0);
		// error
	}

//	memset(ret,0,sizeof(blockSize));
	return ret;
}





static void retMempoolSetKernel(MempoolSet *ptr,void *data, const int blockSize)
{
	if(blockSize > ptr->maxBlockSize && blockSize <= ptr->hugeBlockSize)
	{
//		retMempool(ptr->huge,data);
		free(data);
		ptr->huge->count--;
	}
	else
	{
		const double lgBlockSize = log2(blockSize);
		int index = ceil(lgBlockSize - ptr->lgMinBlockSize);
		if(index < 0) index = 0;
		retMempool(ptr->poolSet[index],data);
	}
}




static void dumpMempoolSetKernel(FILE *fp, const MempoolSet *set)
{
	int i;
	for(i=0;i<set->numOfPool;i++) dumpMempool(fp,set->poolSet[i]);
	dumpMempool(fp,set->huge);
}



//=============================================================





void createMempoolSet(const int minBlockSize,const int numOfMempool,const int hugeBlockSize,const int *const numOfBlockList,const int poolListSize)
{
	int i;
	_poolList = (MempoolSet **)malloc(sizeof(MempoolSet *)*poolListSize);
	for(i=0;i<poolListSize;i++)
	{
		_poolList[i] = createMempoolSetKernelNew(minBlockSize,numOfMempool,hugeBlockSize,numOfBlockList);
	}
	_pool = _poolList[0];
	_poolListSize = poolListSize;
}




void freeMempoolSet(void)
{
	int i;
	for(i=0;i<_poolListSize;i++)
	{
		freeMempoolSetKernel(_poolList[i]);
	}
	free(_poolList);
}



void clearPidMempoolSet(const int pid)
{
	clearMempoolSetKernel(_poolList[pid%_poolListSize]);
}




void *getMempoolSet(const int blockSize)
{
	void *ptr;
	ptr =  getMempoolSetKernel(_pool,blockSize);
	return ptr;
}




void *getPidMempoolSet(const int blockSize,const int pid)
{
	void *ptr;
	ptr =  getMempoolSetKernel(_poolList[pid%_poolListSize],blockSize);
	return ptr;
}



void retMempoolSet(void *data, const int blockSize)
{
	retMempoolSetKernel(_pool,data,blockSize);
}




void retPidMempoolSet(void *data, const int blockSize,const int pid)
{
	retMempoolSetKernel(_poolList[pid%_poolListSize],data,blockSize);
}



void dumpMempoolSet(FILE *fp)
{
	dumpMempoolSetKernel(fp,_pool);
}


void usageMempoolSet(FILE *fp)
{
	int j;

//	fprintf(stderr,"pool list size:%d\n",_poolListSize);
	for(j=0;j<_poolListSize;j++)
	{
		int i;
		fprintf(fp,"poolListIndex: %d\n",j);
		for(i=0;i<_poolList[j]->numOfPool;i++)
		{
//			fprintf(fp,"*\tmalloc count:%2d\n",_poolList[j]->poolSet[i]->malloc_count);
			fprintf(fp,"\tpool:%2d count:%2d\n",i,_poolList[j]->poolSet[i]->count);
		}
		fprintf(fp,"\tpool:hg count:%2d\n",_poolList[j]->huge->count);
		fprintf(fp,"\n");
/*
		fprintf(fp,"\n================================\n");
		fprintf(fp,"memory usage stats\n");
		const int MB = 1024*1024;
		int i;
		int total = 0;
		for(i=0;i<_poolList[j]->numOfPool;i++)
		{
			const int val = _poolList[j]->poolSet[i]->batchNumOfBlock * ((float)_pool->poolSet[i]->blockSize/MB);
			total += val;
			fprintf(fp,"memPool %d byte: %dMB\n",_poolList[j]->poolSet[i]->blockSize,val);
		}

		const int hugeVal = _poolList[j]->huge->batchNumOfBlock * (_poolList[j]->huge->blockSize/MB);
		total += hugeVal;
		fprintf(fp,"memPool huge %d byte: %dMB\n",_poolList[j]->hugeBlockSize,hugeVal);
		fprintf(fp,"total memPool usage: %dMB\n",total);
		fprintf(fp,"================================\n");
*/
	}
}





/*
const int num = 1024*1024;
void *ptr;

void *alloc(void *par)
{
	fprintf(stderr,"enter to alloc\n");
	int i;
	ptr = NULL;

	fprintf(stderr,"enter to alloc\n");
	while(1)
	{
		if(ptr == NULL)
		{
			ptr = getPidMempoolSet(1024*1024*1024,1);
			fprintf(stderr,"alloc\n");
		}
	}
	pthread_exit(0);
}

void *dealloc(void *par)
{
	int i;

	fprintf(stderr,"enter to dealloc\n");
	while(1)
	{
		if(ptr !=NULL) 
		{
			retPidMempoolSet(ptr,1024*1024*1024,1);
			ptr = NULL;
			fprintf(stderr,"dealloc\n");
		}
	}

	pthread_exit(0);
}


int main(void)
{
	int list[17] = {512,512,512,8,8,8,8,8,8,8,8,8,8,8,0,0,0};
	createMempoolSet(16,17,1.5*1024*1024*1024,list,2); // 0 is used for "main thread", threadNum+1 is used for "extra-root thread"

	pthread_t pid[2];
	pthread_create(&pid[0],NULL,alloc,NULL);
	pthread_create(&pid[1],NULL,dealloc,NULL);
		
	pthread_join(pid[0],NULL);
	pthread_join(pid[1],NULL);
	

	freeMempoolSet();
	return 0;
}

*/
