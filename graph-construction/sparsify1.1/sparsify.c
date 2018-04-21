#include <metis.h>

GraphType* getRandomLocalSparsifiedGraph(GraphType* g )
{
	GraphType* retGraph = (GraphType*)malloc(sizeof(GraphType));
	retGraph->nvtxs = g->nvtxs;
	retGraph->nedges = g->nedges;
	retGraph->xadj = (idxtype*)malloc(sizeof(idxtype)*(g->nvtxs+1));
	retGraph->adjncy = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));
	if ( g->adjwgt != NULL )
		retGraph->adjwgt = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));

	int i,j,k;
	j=0;
	for ( i=0; i<retGraph->nvtxs; i++)
	{
		retGraph->xadj[i] = j;
		int d = g->xadj[i+1] - g->xadj[i];
		if ( d <= 0 )
			continue;
		int numToRetain = (int) ceil(sqrt(d*1.0));
		if ( numToRetain > d )
			numToRetain	= d;
		for ( k=g->xadj[i]; k<g->xadj[i+1]; k++ )
		{
			// for each edge, flip a biased coin
			int r = RandomInRangeFast(d);
			if ( r < numToRetain ) // this test ensures proper
			// bias in the coin
			{
				retGraph->adjncy[j] = g->adjncy[k];
				if ( g->adjwgt != NULL )
					retGraph->adjwgt[j] = g->adjwgt[k];
				j++;
			}
			if ( j-retGraph->xadj[i] == numToRetain )
				break;
		}
	}
	retGraph->xadj[retGraph->nvtxs] = j;
	retGraph->nedges = j;

	printf("Num. edges in random locally sparsified graph:%d\n",
	retGraph->nedges);
	return retGraph;
}

GraphType* getExactGlobalSimSparsifiedGraph(GraphType* g, float
percentToRetain)
{
	int i,j, nvtxs = g->nvtxs;
	for (i=0; i<nvtxs; i++)
	{
		iidxsort(g->xadj[i+1]-g->xadj[i], g->adjncy + g->xadj[i]); 
	}

	Matrix* simMatrix = allocMatrix(g->nvtxs, g->nedges, 0, 0,
	0);

	simMatrix->xadj[0] = 0;
	int k=0;
	for ( i=0; i<nvtxs; i++ )
	{
		idxtype* base = g->adjncy + g->xadj[i];
		int degree = g->xadj[i+1] - g->xadj[i];
		k = simMatrix->xadj[i];
		for ( j=g->xadj[i]; j<g->xadj[i+1]; j++ )
		{
			int nbr = g->adjncy[j];
			if ( nbr == i )
				continue;
			int nbr_degree = g->xadj[nbr+1] - g->xadj[nbr];
			int intersect = sizeOfSetIntersect(base, degree,
			g->adjncy+g->xadj[nbr],nbr_degree);
			wgttype sim = (wgttype) (intersect*1.0)/(degree +
			nbr_degree - intersect);
			
			simMatrix->adjncy[k] = nbr;
			simMatrix->adjwgt[k] = sim;
			k++;
		}
		simMatrix->xadj[i+1] = k;
		if ( i > 0 && i % 20000 == 0 )
		{
			printf("Done with %d vertices\n", i);
			fflush(stdout);
		}
	}
	simMatrix->nnz = simMatrix->xadj[nvtxs];

	printf("Done calculating exact similarities.\n");
//	printf("Num edges in initial graph:%d\n", simMatrix->nnz);

	GraphType* retGraph = (GraphType*)malloc(sizeof(GraphType));
	retGraph->nvtxs = g->nvtxs;
	retGraph->nedges = simMatrix->xadj[nvtxs];
	retGraph->xadj = (idxtype*)malloc(sizeof(idxtype)*(g->nvtxs+1));
	retGraph->adjncy =
	(idxtype*)malloc(sizeof(idxtype)*(retGraph->nedges));
	retGraph->adjwgt =
	(idxtype*)malloc(sizeof(idxtype)*(retGraph->nedges));

	convertMatrixToGraph(simMatrix, &(retGraph->xadj),
	&(retGraph->adjncy), &(retGraph->adjwgt), 0);
//	printf("initial graph nedges: %d\n", retGraph->xadj[nvtxs]);
	int totalToRetain = (int)
	round(retGraph->xadj[nvtxs]*percentToRetain);
	exactPruneGraph(retGraph, totalToRetain, 0);

	freeMatrix(simMatrix);

	return retGraph;
}

GraphType* getExactSimSparsifiedGraph(GraphType* g, float
degreeExponent)
{
	int i,j, nvtxs = g->nvtxs;
	for (i=0; i<nvtxs; i++)
	{
		iidxsort(g->xadj[i+1]-g->xadj[i], g->adjncy + g->xadj[i]); 
	}

	Matrix* simMatrix = allocMatrix(g->nvtxs, g->nedges, 0, 0,
	0);

	simMatrix->xadj[0] = 0;
	int k=0;
	for ( i=0; i<nvtxs; i++ )
	{
		idxtype* base = g->adjncy + g->xadj[i];
		int degree = g->xadj[i+1] - g->xadj[i];
		k = simMatrix->xadj[i];
		for ( j=g->xadj[i]; j<g->xadj[i+1]; j++ )
		{
			int nbr = g->adjncy[j];
			if ( nbr == i )
				continue;
			int nbr_degree = g->xadj[nbr+1] - g->xadj[nbr];
			int intersect = sizeOfSetIntersect(base, degree,
			g->adjncy+g->xadj[nbr],nbr_degree);
			wgttype sim = (wgttype) (intersect*1.0)/(degree +
			nbr_degree - intersect);
			
			simMatrix->adjncy[k] = nbr;
			simMatrix->adjwgt[k] = sim;
			k++;
		}
		int d = k - simMatrix->xadj[i];

		int numToRetain = 0;
		if ( degreeExponent == 0.5 )
			numToRetain = (int) ceil(sqrt(d*1.0));
		else
			numToRetain = (int) ceil(pow(d*1.0, degreeExponent));
		if (numToRetain <= 0 )
			numToRetain = 1;
//		fprintf(stderr, "numToRetain for %d is %d\n", (i+1),
//		numToRetain);
		if ( numToRetain < d )
		{
			exactPrune(&d,
			simMatrix->adjncy+simMatrix->xadj[i],
			simMatrix->adjwgt+simMatrix->xadj[i], numToRetain,
			0);
		}
		// this pruning might retain more than numToRetain,
		// since it includes all values that are equal to the
		simMatrix->xadj[i+1] = simMatrix->xadj[i] + numToRetain;
		if ( i > 0 && i % 20000 == 0 )
		{
			printf("Done with %d vertices\n", i);
			fflush(stdout);
		}
	}
	simMatrix->nnz = simMatrix->xadj[nvtxs];

	printf("Done calculating exact similarities.\n");
	printf("Done retaining top edges.\n");
//	printf("Num edges in initial graph:%d\n", simMatrix->nnz);

	GraphType* retGraph = (GraphType*)malloc(sizeof(GraphType));
	retGraph->nvtxs = g->nvtxs;
	retGraph->nedges = simMatrix->xadj[nvtxs];
	retGraph->xadj = (idxtype*)malloc(sizeof(idxtype)*(g->nvtxs+1));
	retGraph->adjncy =
	(idxtype*)malloc(sizeof(idxtype)*(retGraph->nedges));
	retGraph->adjwgt =
	(idxtype*)malloc(sizeof(idxtype)*(retGraph->nedges));

	convertMatrixToGraph(simMatrix, &(retGraph->xadj),
	&(retGraph->adjncy), &(retGraph->adjwgt), 0);
//	printf("initial graph nedges: %d\n", retGraph->xadj[nvtxs]);

	return retGraph;
}


GraphType* getGlobalHashSparsifiedGraph(GraphType* g, 
 float pct, int nHashes)
{
	int numHashes = nHashes;
	int* randoms = (int*)malloc(sizeof(int)*2*numHashes);
	generateRandoms( 2*numHashes, randoms );

	Hashtable * ht = buildHashtable(g, numHashes, randoms);
	free(randoms);
	randoms = NULL;

	printf("Done minwise hashing.\n");
	fflush(stdout);

	GraphType* retGraph = (GraphType*)malloc(sizeof(GraphType));
	retGraph->nvtxs = g->nvtxs;
	retGraph->nedges = g->nedges;
	retGraph->xadj = (idxtype*)malloc(sizeof(idxtype)*(g->nvtxs+1));
	retGraph->adjncy = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));
	retGraph->adjwgt = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));

	int i, j, k;

	j = 0;
	for ( i=0; i<retGraph->nvtxs; i++ )
	{
		retGraph->xadj[i] = j;
		for ( k=g->xadj[i]; k<g->xadj[i+1]; k++ )
		{
			int numMatched = 0, l;
			int from = i, to = g->adjncy[k];
			if ( from == to )
			{
				// we won't consider self-loops
				continue;
			}
			idxtype *fromBase, *toBase;
			fromBase = ht->hashes + from*numHashes;
			toBase = ht->hashes + to*numHashes;
			for ( l=0; l<numHashes; l++ )
			{
				if ( fromBase[l] == toBase[l] )
					numMatched++;
			}
			retGraph->adjwgt[j] = numMatched;
			retGraph->adjncy[j] = to;
			j++;
		}	
	}
	retGraph->xadj[retGraph->nvtxs] = j;

	printf("Done estimating similarities for all edges.\n");
//	printf("Done building sim. weighted graph.\n");
//	printf("Total edges in sim. weighted graph:%d\n",
//	retGraph->xadj[retGraph->nvtxs]);
/*	if ( checkValidUndirectedGraph(retGraph) )
		printf("Sim. weighted graph is a valid undirected graph\n");
	else
		printf("Sim. weighted graph is a valid undirected graph\n");
*/
	fflush(stdout);

	freeHashtable(ht);
	int numToRetain = (int)round(retGraph->xadj[retGraph->nvtxs]*pct);
//	printf("numToRetain:%d\n", numToRetain);
	//exactPruneGraph(retGraph, numToRetain );
//	exactPruneGraph(retGraph, numToRetain, numHashes );
	exactPruneGraphCountingSort(retGraph, numToRetain, numHashes );

//	ddretainTopNeighborsPerNode(retGraph, ht, degreeExponent);
//	printf("Total edges in global sparsified graph:%d\n",
//	retGraph->xadj[retGraph->nvtxs]);

	return retGraph;
}


/* estimate similarity between adjacency lists of each node by
 * min-hashing each node numHashes number of times. (similarity =
 * number of min-hashes that match / numHashes) Retain those
 * edges whose incident nodes have similarity greater than
 * (threshold/numHashes) */
GraphType* getHashSparsifiedGraph(GraphType* g, int numHashes,
float degreeExponent)
{
	int* randoms = (int*)malloc(sizeof(int)*2*numHashes);
	generateRandoms( 2*numHashes, randoms );

	Hashtable * ht = buildHashtable(g, numHashes, randoms);
	free(randoms);
	randoms = NULL;

	printf("Done minwise hashing.\n");
	fflush(stdout);

	GraphType* retGraph = (GraphType*)malloc(sizeof(GraphType));
	retGraph->nvtxs = g->nvtxs;
	retGraph->nedges = g->nedges;
	retGraph->xadj = (idxtype*)malloc(sizeof(idxtype)*(g->nvtxs+1));
	retGraph->adjncy = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));
	retGraph->adjwgt = (idxtype*)malloc(sizeof(idxtype)*(g->nedges));

	int i, j, k;

	j = 0;
	for ( i=0; i<retGraph->nvtxs; i++ )
	{
		retGraph->xadj[i] = j;
		for ( k=g->xadj[i]; k<g->xadj[i+1]; k++ )
		{
			int numMatched = 0, l;
			int from = i, to = g->adjncy[k];
			if ( from == to )
			{
				// we won't consider self-loops
				continue;
			}
			idxtype *fromBase, *toBase;
			fromBase = ht->hashes + from*numHashes;
			toBase = ht->hashes + to*numHashes;
			for ( l=0; l<numHashes; l++ )
			{
				if ( fromBase[l] == toBase[l] )
					numMatched++;
			}
			retGraph->adjwgt[j] = numMatched;
			retGraph->adjncy[j] = to;
			j++;
		}	
	}
	retGraph->xadj[retGraph->nvtxs] = j;

	printf("Done estimating similarities for all edges.\n");
//	printf("Total edges in sim. weighted graph:%d\n",
//	retGraph->xadj[retGraph->nvtxs]);
/*	if ( checkValidUndirectedGraph(retGraph) )
		printf("Sim. weighted graph is a valid undirected graph\n");
	else
		printf("Sim. weighted graph is a valid undirected graph\n");
*/
	fflush(stdout);

	retainTopNeighborsPerNode(retGraph, ht, degreeExponent);
//	printf("Total edges in local sparsified graph:%d\n",
//	retGraph->xadj[retGraph->nvtxs]);

	freeHashtable(ht);
//	exactPruneGraph(retGraph, 0, threshold);
	return retGraph;
}

int countingSortThresholdEntireGraphByAdjwgt(int nedges, idxtype*
adjncy, idxtype* adjwgt, int position, int maxValue)
{
	
	int* counts = (int*)malloc(sizeof(int)*(maxValue+1));
	int* indexes = (int*)malloc(sizeof(int)*(maxValue+1));
	int i,j;

		// temp arrays to hold values as they are being sorted.
//	idxtype* temp_adjncy = idxmalloc(nedges, "");
//	idxtype* temp_adjwgt = idxmalloc(nedges, "");
	int ret = 0;

//	printf("position:%d\n", position);

			for( j=0; j<maxValue+1; j++)
			counts[j] = 0;
		// do actual counting
		for( j=0; j<nedges; j++ )
			counts[adjwgt[j]]++;
		// indexes for each of the possible values.
		// all values of j can only come after all values of j-1
		indexes[0] = 0;
		for( j=1; j<maxValue+1; j++)
		{
			indexes[j] = indexes[j-1] + counts[j-1];
			if ( position < indexes[j] )
			{
				ret = j-1;
				break;
			}
		}
//		printf("ret from countingSortThreshold:%d\n", ret); 
//		free(temp_adjncy);
//		free(temp_adjwgt);
		free(counts);
		free(indexes);
		return ret;
	

}
// we assume that minValue is 0.
void countingSortAdjListsByAdjwgt(GraphType* graph, int maxValue)
{
	int* counts = (int*)malloc(sizeof(int)*(maxValue+1));
	int* indexes = (int*)malloc(sizeof(int)*(maxValue+1));
	int i,j;

	// temp arrays to hold values as they are being sorted.
	idxtype* temp_adjncy = idxmalloc(graph->nvtxs, "");
	idxtype* temp_adjwgt = idxmalloc(graph->nvtxs, "");

	for (i=0; i<graph->nvtxs; i++)
	{
		// initialize counts to 0
		for( j=0; j<maxValue+1; j++)
			counts[j] = 0;
		// do actual counting
		for( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
			counts[graph->adjwgt[j]]++;
		// indexes for each of the possible values.
		// all values of j can only come after all values of j-1
		indexes[0] = 0;
		for( j=1; j<maxValue+1; j++)
			indexes[j] = indexes[j-1] + counts[j-1];

		// put elements into appropriate positions using indexes,
		// which is suitably maintained.
		for( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
		{
			temp_adjncy[indexes[graph->adjwgt[j]]] =
			graph->adjncy[j];
			temp_adjwgt[indexes[graph->adjwgt[j]]] =
			graph->adjwgt[j];
			indexes[graph->adjwgt[j]]++;
		}

		// copy sorted values back.
		for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
		{
			graph->adjncy[j] = temp_adjncy[j-graph->xadj[i]];
			graph->adjwgt[j] = temp_adjwgt[j-graph->xadj[i]];
		}

	}
	free(temp_adjncy);
	free(temp_adjwgt);

}


/* Sorts the neighbor list of each node by edge weight, and
 * retains the top (log d + 1) neighbors for each node, where d
 * is the degree. */
void retainTopNeighborsPerNode(GraphType* graph, Hashtable *ht,
float exponent)
{
	int i, j, k;
/*
	for ( i=0, j=0; i<graph->nvtxs; i++ )
	{
		int d = graph->xadj[i+1]-graph->xadj[i];
		if ( d == 0 )
		{
			continue;
		}

		int oldStart = graph->xadj[i];
		graph->xadj[i] = j;

		int numToRetain = 0;
		numToRetain = (int) ceil(log(d*1.0) + 1.0);
		if ( numToRetain > d )
			numToRetain	= d;
		int thresh;
		if ( numToRetain == d )
			thresh = 0;
		else
		{
			// randomselect also has the side-effect of placing
			// the topmost numToRetain elements in the right end
			// of the array 
			thresh = RandomSelectInts(graph->adjwgt, oldStart,
						graph->xadj[i+1]-1, d-numToRetain+1);
		}

		// side-effect of randomselect exploited here.
		for ( k=graph->xadj[i+1]-numToRetain; k<graph->xadj[i+1]; k++ )
		{
			if ( graph->adjwgt[k] >= thresh )
			{
				graph->adjncy[j] = graph->adjncy[k];
				graph->adjwgt[j++] = graph->adjwgt[k];
			}
			if ( j-graph->xadj[i] >= numToRetain )
				break;
		}
	}

	graph->xadj[graph->nvtxs] = j;
	printf("Done retaining top neighbors.\n");
	fflush(stdout);
*/

/*	idxtype *temp_adjncy, *temp_adjwgt;
	temp_adjncy = idxmalloc(graph->xadj[graph->nvtxs], "");
	temp_adjwgt = idxmalloc(graph->xadj[graph->nvtxs], "");
	for ( i=0; i<graph->xadj[graph->nvtxs]; i++ )
	{
		temp_adjncy[i] = graph->adjncy[i];
		temp_adjwgt[i] = graph->adjwgt[i];
	}
*/

	timer tmr;
/*	cleartimer(tmr);
	starttimer(tmr);
	for( i=0; i<graph->nvtxs; i++ )
	{
		ParallelQSortInts(graph->adjncy, graph->adjwgt,
		graph->xadj[i], graph->xadj[i+1]-1);
		if ( i % 5000 == 0 )
		{
			printf("%d..", i);
			fflush(stdout);
		}
	}
	stoptimer(tmr);
	printf("Time for sort:%.2f\n", gettimer(tmr));
	fflush(stdout);
*/
	cleartimer(tmr);
	starttimer(tmr);
	countingSortAdjListsByAdjwgt(graph,ht->numHashes); 
/*	for ( i=0; i<graph->nvtxs; i++ )
	{
		if ( i == 1 )
		{
			printf("Before sorting\n");
			for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
				printf("%d %d\n", graph->adjncy[j],	graph->adjwgt[j]);
		}
		ParallelQSortInts(graph->adjwgt, graph->adjncy, graph->xadj[i],
		graph->xadj[i+1] - 1);
		if ( i % 5000 == 0 )
		{
			printf("%d..", i);
			fflush(stdout);
		}
		if ( i == 1 )
		{
			printf("After sorting\n");
			for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
				printf("%d %d\n", graph->adjncy[j],	graph->adjwgt[j]);
			exit(0);
		}
	}
	*/
	stoptimer(tmr);
	printf("Time for counting sort:%.2f\n", gettimer(tmr));
	printf("Done sorting.\n");
	fflush(stdout);

	for ( i=0, j=0; i<graph->nvtxs; i++ )
	{
		k = graph->xadj[i];
		graph->xadj[i] = j;
		int d = graph->xadj[i+1] -k;
		int numToRetain = 0;
		if ( d > 0 )
		{
			if ( exponent == 0.5 )
				numToRetain = (int) ceil(sqrt(d*1.0));
			else
				numToRetain = (int) ceil(pow(d*1.0, exponent));
			if ( numToRetain <= 0 )
				numToRetain = 1;
	//		numToRetain = (int) ceil(log(d*1.0) + 1.0);
//			numToRetain *= (int) round(numToRetain * 1.25);
		}
		if ( numToRetain > d )
			numToRetain	= d;
//		fprintf(stderr, "numToRetain for %d is %d\n", (i+1),
//		numToRetain);
		int skip = d - numToRetain;
		for ( int l = 0; l < numToRetain; l++, j++ )
		{
		// top neighbors are at the end of the list, since
		// parallelqsort sorts in ascending order.
			graph->adjncy[j] = graph->adjncy[k+skip+l];
			graph->adjwgt[j] = graph->adjwgt[k+skip+l];
		}
	}
	graph->xadj[graph->nvtxs] = j;
	printf("Done retaining top neighbors.\n");
	fflush(stdout);

}

void exactPruneGraphCountingSort(GraphType* graph, int l, int
maxValue)
{
	idxtype* temp_a=idxmalloc(graph->nedges, "");
	idxtype* temp_w=idxmalloc(graph->nedges, "");
	int i;
	for ( i=0; i<graph->nedges; i++)
	{
		temp_a[i] = graph->adjncy[i];
		temp_w[i] = graph->adjwgt[i];
	}
	
	int nedges = graph->xadj[graph->nvtxs];
	int thresh = countingSortThresholdEntireGraphByAdjwgt(nedges,
	temp_a, temp_w, nedges-l-1, maxValue);

	free(temp_a);
	free(temp_w);
	exactPruneGraph(graph, 0, thresh);
}

void exactPruneGraph(GraphType* graph, int l, idxtype thresh )
{
/*	Matrix M;
	int i;
	M.nnz = graph->nedges;
	M.adjncy = graph->adjncy;
	M.adjwgt = (wgttype*) malloc( sizeof(wgttype)*M.nnz );
	for ( i=0; i<M.nnz; i++)
		M.adjwgt[i] = (wgttype)graph->adjwgt[i];

	exactPruneMatrix(&M, l);

	graph->nedges = M.nnz;
	graph->adjncy = M.adjncy;
	for ( i=0; i<M.nnz; i++)
		graph->adjwgt[i] = (idxtype) M.adjwgt[i];
	graph->adjwgt = (idxtype*) realloc (graph->adjwgt, sizeof(wgttype)*M.nnz);

	return;
*/
	int i, j,k;
	long nnz = graph->nedges;

	if ( thresh == 0 )
	{
		long s = sizeof(idxtype)*(nnz);
		idxtype* newadjwgt=(idxtype*)malloc(s);

		for(i = 0; i < nnz; i++)
			newadjwgt[i]=graph->adjwgt[i];
			
		thresh = RandomSelectInts(newadjwgt, 0, nnz-1, nnz-l+1);
		free(newadjwgt);
	}

//	printf("Threshold:%d\n",thresh);

	for(i=0,j=0; i<graph->nvtxs; i++)
	{
		k=graph->xadj[i];
//		if ( i > 0 )
//		{
			graph->xadj[i]=j;
//			printf("M->xadj[%d] set to %d\n", i, j);
//			printf("M->xadj[%d] is %d\n",i+1,M->xadj[i+1]);
//		}

		for(; k<graph->xadj[i+1]; k++)
		{
//			printf("M->adjwgt[%d]:%.2f\n",k,M->adjwgt[k]);
			if ( graph->adjwgt[k] >= thresh 
				|| graph->adjncy[k]	== i )
			{
				if ( j < k )
				{
					graph->adjncy[j]=graph->adjncy[k];
					graph->adjwgt[j]=graph->adjwgt[k];
				}
				j++;
			}
		}
//		printf("j:%d\n",j);
	}
	graph->xadj[graph->nvtxs]=j;
//	printf("M->nnz set to %d\n", j);
	
	if ( j < graph->nedges )
	{
		graph->nedges = j;
		graph->adjncy = (idxtype*) realloc(graph->adjncy,j*sizeof(idxtype));
		graph->adjwgt =	(idxtype*) realloc(graph->adjwgt,j*sizeof(idxtype));
	}

}

void convertMatrixToGraph(Matrix*M, idxtype **ret_xadj,
idxtype **ret_adjncy, idxtype **ret_adjwgt, wgttype scale)
{
	wgttype min, max, avg;

	int i, j, k;
	int minindex, minadj;
	min=max=avg=M->adjwgt[0]; 
	for(i=1;i<M->nnz;i++)
	{
//		M->adjwgt[i] /= 2;
		avg += M->adjwgt[i];
		if ( M->adjwgt[i] < min )
		{
			min = M->adjwgt[i];
			minindex = i;
			minadj = M->adjncy[i];
		}
		if ( M->adjwgt[i] > max )
			max = M->adjwgt[i];
	}
	avg /= M->nnz;

	if ( scale <= 0 )
	{
		if ( max/min > 1000 )
		{
			scale = 1.0/min;
			if ( min < 1.0e-07 )
				scale = 1.0e+05;
		}
		else
			scale = 1000/max; // scale = (1/(max/1000))
	}

//	printf("Statistics on conversion from matrix to graph:\n");
//	printf("min:%f, max:%f, avg:%f, nnz:%d, scale:%f\n", 
//				min, max, avg, M->nnz, scale );
//	printf("minindex:%d, minadj:%d\n", minindex, minadj);

	*ret_adjwgt=idxmalloc(M->nnz,
					"convertMatrixToGraph: new adjacency weights");
	*ret_adjncy=idxmalloc(M->nnz,"new adjacency lists");
	*ret_xadj=idxmalloc(M->nvtxs+1,"new xadj");
	(*ret_xadj)[0]=0;

//	printf("M->adjncy:%d,ret_adjncy:%d,M->xadj[0]:%d\n",
//					M->adjncy!=NULL,(*ret_adjncy)!=NULL,M->xadj[0]);

	for(k=0;k<M->nvtxs;k++)
	{
		for(i=M->xadj[k],j=(*ret_xadj)[k];i<M->xadj[k+1];i++)
		{
			if (M->adjncy[i]==k)
				continue;
			wgttype wt=scale*M->adjwgt[i];
			if ( wt >= 1 )
			{
//				printf("j:%d i:%d\n",j,i);
				(*ret_adjncy)[j]=M->adjncy[i];
				(*ret_adjwgt)[j]=(idxtype)round(wt);
				j++;
			}
		}
//		if ( k % 1000 == 0 )
//			printf("%d..",k);

		// if node has no neighbours, then add self-loop
		if ( (*ret_xadj)[k] == j )
		{
	//		(*ret_adjncy)[j]=k;
	//		(*ret_adjwgt)[j++]=1;
		}

		(*ret_xadj)[k+1]=j;
	}
}

// This one implements A + A'
void usual_symmetrize_directed(int nvtxs, int nedges, idxtype*
		xadj, idxtype* adjncy, idxtype* adjwgt, idxtype** ret_xadj,
		idxtype** ret_adjncy, idxtype** ret_adjwgt, int	addSelfLoops)
{
	Matrix *M, *M1, *M2;
	int i,j,k,doubleEdges, noOfSelfLoops;

	M1=(Matrix*)malloc(sizeof(Matrix));
	M1->nvtxs=nvtxs;M1->nnz=nedges;
	M1->xadj=xadj;M1->adjncy=adjncy;
	M1->adjwgt=(wgttype*)malloc(sizeof(wgttype)*nedges);
	if ( adjwgt != NULL )
	{
		for(i=0;i<nedges;i++)
			M1->adjwgt[i]=adjwgt[i];
	}
	else
	{
		for(i=0;i<nedges;i++)
			M1->adjwgt[i]=1.0;
	}
	// the sorting is so the addition later on will work.
	sortAdjLists(nvtxs, M1->xadj, M1->adjncy, M1->adjwgt);
	// getTranspose sorts M2, so we needn't worry about it.
	M2=getTranspose(M1);
	sortAdjLists(nvtxs, M2->xadj, M2->adjncy, M2->adjwgt);
	
	M=add(M1,M2);
	free(M1->adjwgt);
	free(M1);
	freeMatrix(M2);

	if ( addSelfLoops )
	{
		*ret_adjwgt=idxmalloc(M->nnz+M->nvtxs,
						"usual_symmetrize_directed: new adjwgt");
		*ret_adjncy=idxmalloc(M->nnz+M->nvtxs,
						"usual_symmetrize_directed:new adjlist");
	}
	else
	{
		*ret_adjwgt=idxmalloc(M->nnz,
						"usual_symmetrize_directed: new adjwgt");
		*ret_adjncy=idxmalloc(M->nnz,
						"usual_symmetrize_directed:new adjlist");
	}

	*ret_xadj=idxmalloc(M->nvtxs+1,"new xadj");
	(*ret_xadj)[0]=0;

	doubleEdges=0,noOfSelfLoops=0;
	for(k=0;k<M->nvtxs;k++)
	{
		int seenSelfloop=0;
		j=(*ret_xadj)[k];
		// self-loop
		if ( addSelfLoops )
		{
			(*ret_adjncy)[j]=k;
			(*ret_adjwgt)[j++]=1;
		}
		for(i=M->xadj[k];i<M->xadj[k+1];i++)
		{
	/*		if ( M->adjncy[i] > k 
				&& ( i==M->xadj[k] || M->adjncy[i-1] < k ) )
			{
				seenSelfloop=1;
				(*ret_adjncy)[j]=k;
				(*ret_adjwgt)[j++]=1;
			}
			if (M->adjncy[i]==k)
			{
				if ( seenSelfloop == 1 )
				{
					printf("Yikes, seeing self loop for the");
					printf(" second time!, nodeid:%d\n", k);
					abort();
				}
				noOfSelfLoops++;
				seenSelfloop=1;
			} */
			if (M->adjncy[i]==k)
				continue;
				
			wgttype wt=M->adjwgt[i];
			//if ( wt >= 1 )
			//{
				(*ret_adjncy)[j]=M->adjncy[i];
				(*ret_adjwgt)[j]=(idxtype)wt;
				j++;
				if ( wt > 1 )
					doubleEdges++;
			//}
		}

		(*ret_xadj)[k+1]=j;
	}
	doubleEdges/=2;
//	printf("2-way edges:%d, noOfSelfLoops:%d\n",doubleEdges, noOfSelfLoops);

	freeMatrix(M);
//	free(M1->adjwgt);
//	free(M1);
//	freeMatrix(M2);

}

int sizeOfSetIntersect(idxtype* a, int a_length, idxtype* b, int
b_length, int sort)
{
	if ( sort > 0 )
	{
		iidxsort(a_length, a);
		iidxsort(b_length, b);
	}

	int numMatches = 0;
	for ( int i=0, j=0; i<a_length && j < b_length; )
	{
		if ( a[i] == b[j] )
		{
			numMatches++;
			i++; j++;
		}
		else if ( a[i] > b[j] )
		{
			j++;
		}
		else 
			i++;
	}
	return numMatches;
}

/* This function can be used only when the adjacency lists for
 * each vertex are sorted! Beware!
 **/
Matrix* add(Matrix* M0, Matrix* M1)
{
	if ( checkSorted(M0->nvtxs, M0->xadj, M0->adjncy) == 0 || 
			checkSorted(M1->nvtxs, M1->xadj, M1->adjncy) == 0)
	{
		printf("Yikes! Input to matrix add is not sorted!\n");
		fflush(stdout);
		return NULL;
	}

	Matrix* ret;
	int i,j,k; 
	long ret_size=(M0->nnz > M1->nnz)? M0->nnz: M1->nnz ;
	int init_ret_size = ret_size;
	wgttype t;

	ret=allocMatrix(M0->nvtxs, init_ret_size ,0,0,0);
	ret->xadj[0]=0;
	for(i=0;i<M0->nvtxs;i++)
	{
		int l=ret->xadj[i];
		for(j=M0->xadj[i],k=M1->xadj[i];
			 j<M0->xadj[i+1] && k<M1->xadj[i+1];)
		{
			wgttype a=M0->adjwgt[j],b=M1->adjwgt[k];

			if ( l+2 >= ret_size )
			{
				ret_size+=ret_size;
				ret->adjncy=(idxtype*)realloc(ret->adjncy,
								(ret_size)*sizeof(idxtype));
				ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
								(ret_size)*sizeof(wgttype));
				if ( ret->adjncy == NULL || ret->adjwgt == NULL )
				{
					printf("Could not allocate");
					printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
					abort();
				}
			}

			if (M0->adjncy[j]==M1->adjncy[k])
			{
				ret->adjwgt[l]=a+b;
				ret->adjncy[l]=M0->adjncy[j];
				j++;k++;l++;
			}
			else 
			{
				if ( M0->adjncy[j] < M1->adjncy[k] )
				{
					ret->adjwgt[l]=a;
					ret->adjncy[l]=M0->adjncy[j];
					j++;l++;
				}
				else
				{
					ret->adjwgt[l]=b;
					ret->adjncy[l]=M1->adjncy[k];
					k++;l++;
				}
			}
		}
		if ( l+(M0->xadj[i+1]-j)+(M1->xadj[i+1]-k) >= ret_size )	
		{
			ret_size+=ret_size;
			ret->adjncy=(idxtype*)realloc(ret->adjncy,
							(ret_size)*sizeof(idxtype));
			ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
							(ret_size)*sizeof(wgttype));
			if ( ret->adjncy == NULL || ret->adjwgt == NULL )
			{
				printf("Could not allocate");
				printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
				abort();
			}
		}

		for(;j<M0->xadj[i+1];j++)
		{
			ret->adjncy[l]=M0->adjncy[j];
			ret->adjwgt[l++]=M0->adjwgt[j];
		}
		for(;k<M1->xadj[i+1];k++)
		{
			ret->adjncy[l]=M1->adjncy[k];
			ret->adjwgt[l++]=M1->adjwgt[k];
		}

		ret->xadj[i+1]=l;
	}

	ret->nnz=ret->xadj[ret->nvtxs];

	return ret;
}

Matrix* getTranspose2(Matrix* M)
{
	if ( M->nvtxs == 0 || M->xadj[M->nvtxs] == 0 )
	{
		printf("Matrix getTranspose invoked with empty matrix\n");
		return M;
	}

	int nvtxs = M->nvtxs;
	Matrix* ret=allocMatrix(M->nvtxs, M->xadj[nvtxs], 0, 0, 0);
	idxtype *inDegrees = idxmalloc(M->nvtxs,
	"getTranspose2:inDegrees");
	for ( int i=0; i<nvtxs; i++ )
		inDegrees[i] = 0;
	for ( int i=0; i<M->xadj[nvtxs]; i++ )
	{
/*		if ( M->adjncy[i] < 0 || M->adjncy[i] > M->nvtxs-1 )
		{
			printf("Yikes! Mistake in input, edge to %d\n",
			M->adjncy[i]);
			abort();
		}
*/		inDegrees[M->adjncy[i]]++;
	}

	ret->xadj[0]=0;
	for ( int i = 0; i < nvtxs; i++)
		ret->xadj[i+1] = ret->xadj[i] + inDegrees[i];
	
	idxtype *counters = inDegrees;
	for ( int i = 0; i < nvtxs; i++)
		counters[i] = ret->xadj[i];

	for ( int i = 0; i < M->nvtxs; i++ )
	{
		for ( int j = M->xadj[i]; j < M->xadj[i+1]; j++ )
		{
			int k = M->adjncy[j];
			ret->adjncy[counters[k]] = i;
			ret->adjwgt[counters[k]] = M->adjwgt[j];
			counters[k]++;
		}
	}
	return ret;
}

Matrix* getTranspose(Matrix* M)
{
	return getTranspose2(M);
}


/*
 * Does exact pruning of one column. n indicates number of column
 * entries, adjncy is an array of node ids, adjwgt is the array
 * of flow values and k is the number of entries that need to be
 * retained. The return value is the sum of the retained entries.
 */
wgttype exactPrune(int *n, idxtype* adjncy, wgttype* adjwgt,
						int k, wgttype initsum)
{
	int i,j,initn=*n;
	wgttype sum=0, *newadjwgt,thresh;
	// if column size less than k, don't need to do anything.
	if ( k >= *n )
		return initsum;

	// RandomSelect moves elements in the array around. Hence,
	// need a placeholder array.
	newadjwgt=(wgttype*)malloc(sizeof(wgttype)*(*n));
	for( i=0; i<*n; i++ )
	{
		newadjwgt[i]=adjwgt[i];
	}

/*	qsort((void*)newadjwgt, *n, sizeof(wgttype),
	comparewgttypes); */
//	thresh=newadjwgt[k-1];

	thresh=RandomSelect(newadjwgt,0,*n-1,*n-k+1);
//	printf("The %d biggest element selected by exactPrune",k);
//	printf(" is %f\n", thresh);
	free(newadjwgt);

	for(i=0,j=0;i<*n;i++)
	{
		if  ( adjwgt[i] >= thresh )
		{
			if ( j < i )
			{
				adjncy[j]=adjncy[i];
				adjwgt[j]=adjwgt[i];
			}
			sum += adjwgt[i];
			j++;
		}
	}
	*n=j;

/*	if (!( sum > 0)  || sum > 1 )
	{
		printf("sum:%f,n:%d,k:%d,initsum:%f,thresh:%f\n",
					sum,*n,k,initsum,thresh); 
		abort();
	}
*/
//	free(newadjwgt);
	return sum;
}
