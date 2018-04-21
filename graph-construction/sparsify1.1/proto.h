/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.43 2011-04-20 20:02:20 venu Exp $
 *
 */

/* bucketsort.c */
void BucketSortKeysInc(int, int, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void my_CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, int, idxtype *, idxtype *);
GraphType *my_SetUpCoarseGraph(GraphType *, int);
GraphType *SetUpCoarseGraph(GraphType *, int, int);
void ReAdjustMemory(GraphType *, GraphType *, int);
void CreateCoarseGraph_PowerLaw(CtrlType*, GraphType*, int,
		idxtype*, idxtype*, idxtype* );

/* Venu: My addition */
/* mis_coarsen.c */
Matrix** mis_Coarsen2Way(CtrlType *, GraphType*, int, int*);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* debug.c */
void fitToPiecewiseUniform(float*, int, float, float, float*,
float*);
double betaHalfToOneCdf(double, double, double);
float convertCosineForHashing(float);
void convertCosineForHashing(float *, int, float *);
void fitToBeta(float*, int, float*, float*, float lowerlimit=0,
float upperlimit=1);
double ahkEdges(GraphType*, float);
int numberOfMatchingEdges(GraphType*, GraphType*);
float* getOrderStatistics(float*, int, int);
float stdDeviation(idxtype*, int);
idxtype* getWeightsHistogram(GraphType*, int*, int);
idxtype* getDegreeHistogram(GraphType*, int*, int);
void ComputeAdjWgtSums(GraphType*);
int ComputeNumEdgesCut(GraphType*, idxtype*);
idxtype ComputeVolSubset(GraphType*, idxtype*);
int ComputeCut(GraphType *, idxtype *);
float ComputeConductance(GraphType *, idxtype *, int );
float ComputeNCut(GraphType *, idxtype *, int );
float Overlap_ComputeNcutVector(GraphType*, idxtype*, idxtype*,
idxtype, float*);
float ComputeNCutVector(GraphType*, idxtype*, int, float*);
int mapPartition(idxtype*, idxtype);
int CheckBnd(GraphType *);
int CheckBnd2(GraphType *);
int CheckNodeBnd(GraphType *, int);
int CheckRInfo(RInfoType *);
int CheckNodePartitionParams(GraphType *);
int IsSeparable(GraphType *);
idxtype* histogram(idxtype*, int, int);
//long* histogramForLongs(long*, long, long);
int checkValidUndirectedGraph(GraphType*);
int checkEqualityOfGraphs(GraphType*);

/* estmem.c */
void METIS_EstimateMemory(int *, idxtype *, idxtype *, int *, int *, int *);
void EstimateCFraction(int, idxtype *, idxtype *, float *, float *);
int ComputeCoarseGraphSize(int, idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *);

/* graph.c */
void my_SetUpGraph(GraphType*, int, idxtype*, idxtype*, idxtype*,
idxtype*, int, int);
void SetUpGraph(GraphType *, int, int, int, idxtype *, idxtype *, idxtype *, idxtype *, int);
void SetUpGraphKway(GraphType *, int, idxtype *, idxtype *);
void SetUpGraph2(GraphType *, int, int, idxtype *, idxtype *, float *, idxtype *);
void VolSetUpGraph(GraphType *, int, int, int, idxtype *, idxtype *, idxtype *, idxtype *, int);
void RandomizeGraph(GraphType *);
int IsConnectedSubdomain(CtrlType *, GraphType *, int, int);
int IsConnected(CtrlType *, GraphType *, int);
int IsConnected2(GraphType *, int);
int FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);
void my_Match_HEMN(CtrlType *, GraphType *);
void my_Match_SHEMN(CtrlType *, GraphType *);
void my_Match_RM(CtrlType *, GraphType *);
void my_Match_PowerLaw_FC(CtrlType *, GraphType *);


/* mwHash.c */
idxtype pickRandomElement(idxtype*, int);
idxtype getMinHashKey(idxtype*, int, int*, int, int);
Hashtable* buildHashtable(GraphType*, int, int*);
void freeHashtable(Hashtable *);
void generateRandoms(int, int*);
idxtype hash(idxtype, idxtype, idxtype, idxtype);

/* hash_match.c */
void generateRandoms(int, int*);

GraphType* getGlobalHashSparsifiedGraph(GraphType*, float, int);
GraphType* getExactGlobalSimSparsifiedGraph(GraphType*, float);
GraphType* getExactSimSparsifiedGraph(GraphType*, float);
GraphType* getRandomLocalSparsifiedGraph(GraphType*);
GraphType* getHashSparsifiedGraph(GraphType*, int, float);
void getHashGraph(int, int, idxtype*, idxtype*, idxtype*,
			idxtype**, idxtype**, int, int);
void match_hash(CtrlType*, GraphType*, int);
void freePrimes();

/* BetaCdf.c */
void my_gsl_error_handler(const char*, const char*, int, int);
double gsl_betaCdf(double, double, double);
double betaHalfToOneCdf(double, double, double);

/* convert_directed.c */
void addToMatrix(Matrix*, const int, const int, const wgttype);
Matrix* normBcMatrix(int, int*, int*, int*, int*, int*, float,
float);
void l2Normalize(Matrix*);
Matrix* reorderDimensionsAndVectors(Matrix*, idxtype*, wgttype*,
wgttype*, idxtype*);


Matrix* bayesianSimSearch(Matrix*, Beta, float, float, float,
float);
float *sampleUnprunedPairs_wrapper(GraphType *, int, float);
float* sampleSimilarities(GraphType*, int);
void getInAndOutDegrees(int, idxtype*, idxtype*, idxtype**,
idxtype**);
void convertMatrixToGraph(Matrix*,
idxtype**,idxtype**,idxtype**,wgttype scale=0);
void addGraphs(int, idxtype*,idxtype*,idxtype*,idxtype*,idxtype*,
idxtype*,idxtype** ,idxtype**,idxtype**);
void symmetrize_pagerank(GraphType*, idxtype**, idxtype**,
idxtype**, PageRankOptions);
void getReverseGraph(int, int, idxtype*, idxtype*, idxtype*,
idxtype**, idxtype**, idxtype**);
void convert_directed(int*, idxtype*, idxtype*,
idxtype*, idxtype*, idxtype*, idxtype*, idxtype**,
idxtype**, idxtype**);
void symmetrize_directed(const int, idxtype*, idxtype*,
idxtype*, idxtype*, idxtype*, idxtype*, idxtype**,
idxtype**, idxtype**, DirToUndirOptions);
void usual_symmetrize_directed(int, int, idxtype*, idxtype*,
		idxtype*, idxtype**, idxtype**, idxtype**, int);
void *bcWithThreshold_thread(void*);
Matrix* bcWithThreshold(int, idxtype*, idxtype*, idxtype*,
		float); 

/* mclutils.c */
int checkSorted(int, idxtype*, idxtype*);
void removeSelfLoops(GraphType*);
Matrix* removeSelfLoops(Matrix*, int);
void testCoarsening(int*, idxtype*, idxtype*, idxtype*, idxtype*,
int*, int);
void dump_graph(GraphType*);
void ncutifyWeights(Matrix*, int, int);
void normalizeColumns(Matrix*, int, int);
void printRow(Matrix*, int);
int numWithoutAttractors(Matrix*);
void getAttractorsForAll(Matrix*);
int isConverged(idxtype*, idxtype*, int);
int isConverged2(idxtype*, idxtype*, int);
int noOfUnique(idxtype*, int, idxtype*);
void freeMatrix(Matrix*);
Matrix* allocMatrix(int, int, int, int, int);
void dumpMatrix(Matrix*);
Matrix* setupCanonicalMatrix(int, int, idxtype*, idxtype*,
idxtype*, int );
void sortAdjLists(int, idxtype*, idxtype*, wgttype*);
void getPermutedGraph(idxtype*, idxtype*, int, int, idxtype*,
	idxtype*, idxtype*, idxtype**, idxtype**, idxtype** );
Matrix* permuteRowsAndColumns(const Matrix*, const idxtype*,
const idxtype*);

/* org_mcl.c */
void org_mcl(int*, idxtype*, idxtype*, idxtype*, idxtype*, 
int*, idxtype*,float,int,int, int, int);

/* rmcl.c */
/*void rmcl(int*, idxtype*, idxtype*, idxtype*, idxtype*, 
int*, idxtype*,float,int,int, int, int);
*/
Matrix* dprmcl(Matrix*, Matrix*, GraphType*, Options, int, int);
void dprmclWrapper(int*, idxtype*,idxtype*,idxtype*, 
idxtype*, int*, idxtype*, Options);
 

/* mlmcl.c */
void mlmcl(int*, idxtype*, idxtype*, idxtype*, idxtype*, int*,
				idxtype*, Options);

/* pagerank.c */
idxtype* maiyaSelectNodes(GraphType*, int*);
Matrix* commuteTimes(GraphType*, int);
GraphType* forestFireSparsify(GraphType*, float, float);
void initPageRankOptions(PageRankOptions *opt);
wgttype* pagerank(GraphType*, PageRankOptions, Matrix *M = NULL); 

/* mclbase.c */
void exactPruneGraphCountingSort(GraphType*, int, int);
int doTheyIntersect(idxtype*, int, idxtype*, int, int sort=0);
wgttype dotProduct(idxtype*, wgttype*, int, idxtype*, wgttype*,
				int, int sort =0);
int sizeOfSetIntersect(idxtype*, int, idxtype*, int, int sort=0);
void knnMatrix(Matrix*, int);
void checksums(Matrix*);
Matrix* getTranspose(Matrix*);
Matrix* getTranspose2(Matrix*);
Matrix* AtimesATransposeBlocking(Matrix* , wgttype*, wgttype*,
					int, wgttype);
Matrix* AtimesATranspose(Matrix*, wgttype*, wgttype*, int, int,
			int, int, wgttype, timer*);
Matrix* expand(Matrix*,Matrix*);
Matrix* expand_ht(Matrix*,Matrix*,idxtype*, int, wgttype
threshold = 0);
Matrix* allInOneStep(Matrix*, Matrix*, idxtype*,Options ,int, int);
Matrix* getDprAdjMatrix(Matrix*,Matrix*,idxtype*,wgttype);
Matrix* add(Matrix*, Matrix*);
void changeBetweenMatrices(Matrix*,Matrix*,wgttype*);
void changeBetweenMatrices2(Matrix*,Matrix*,wgttype*);
void inflate(Matrix*,float);
void pruneAndNormalize(Matrix*,int,int);
wgttype exactPrune(int*, idxtype*, wgttype*, int, wgttype);
void exactPruneMatrix(Matrix*, int);
void exactPruneGraph(GraphType*, int, idxtype thresh=0);
void retainTopNeighborsPerNode(GraphType*, Hashtable* , float );
int compareints(const void *, const void *);

/* memory.c */
void my_AllocateWorkSpace(CtrlType *ctrl, GraphType *graph);
void AllocateWorkSpace(CtrlType *, GraphType *, int);
void FreeWorkSpace(CtrlType *, GraphType *);
int WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, int);
void idxwspacefree(CtrlType *, int);
float *fwspacemalloc(CtrlType *, int);
void fwspacefree(CtrlType *, int);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* myqsort.c */
void iidxsort(int, idxtype *);
void iintsort(int, int *);
void ikeysort(int, KeyValueType *);
void ikeyvalsort(int, KeyValueType *);

/* subgraphs.c */
int selfLoopifyHubs(int, idxtype*, idxtype*, idxtype*, wgttype*,
int);
idxtype* removeHubs(GraphType*, int, int, GraphType**, int);
void getSubgraph(GraphType*, idxtype*, int, int, GraphType**);
GraphType* getCutGraph(GraphType*, idxtype*, int);
void globallySampleEdges(int, int, idxtype*, idxtype*, idxtype*,
		idxtype**, idxtype**, float);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* util.c */
void RandomPermuteFloatsInts(int, wgttype*, idxtype*);
void assignClustersToHubs(idxtype*, idxtype*, int, int,
GraphType*);
void dfTraversalMatrix(Matrix*, idxtype, idxtype*, int*, wgttype
minWgt);
void dfTraversal(GraphType*, idxtype, idxtype*, int*);
void initOptions(Options*);
void mapIndices(idxtype*, idxtype*, int, int);
idxtype* lookForSingletons(GraphType*, int*);
void errexit(char *,...);
#ifndef DMALLOC
int *imalloc(int, const char *);
idxtype *idxmalloc(int, const char *);
float *fmalloc(int, const char *);
int *ismalloc(int, int, const char *);
idxtype *idxsmalloc(int, idxtype, const char *);
//void *GKmalloc(int, char *);
void *GKmalloc(long, const char *);
#endif
void GKfree(void **,...); 
int *iset(int n, int val, int *x);
idxtype *idxset(int n, idxtype val, idxtype *x);
float *sset(int n, float val, float *x);
int iamax(int, int *);
int idxamax(int, idxtype *);
int idxamax_strd(int, idxtype *, int);
int samax(int, float *);
int samax2(int, float *);
int idxamin(int, idxtype *);
int samin(int, float *);
int idxsum(int, idxtype *);
int idxsum_strd(int, idxtype *, int);
void idxadd(int, idxtype *, idxtype *);
int charsum(int, const char *);
int isum(int, int *);
float ssum(int, float *);
float ssum_strd(int n, float *x, int);
void sscale(int n, float, float *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, int, float *, int);
void ParallelQSortFloatsInts(wgttype*, idxtype*, int, int);
void ParallelQSortIntsUsingScores(idxtype*, idxtype*, idxtype*,
			int, int);
void ParallelQSort(idxtype*,wgttype*,int,int);
void ParallelQSortInts(idxtype*,idxtype*,int,int);
void QSortIntsUsingInts(idxtype*, idxtype*, int, int);
void ParallelQSortLongs(long*,wgttype*,int,int);
int bsearch_insertPos(idxtype*, int, int, int);
void RandomPermute(int, idxtype *, int);
void permuteDegreeOrder(int, idxtype*, idxtype*);
wgttype RandomSelect(wgttype*, int, int, int);
idxtype RandomSelectInts(idxtype*, int, int, int);
double drand48();
void srand48(long);
int ispow2(int);
void InitRandom(int);
int log2(int);

/* io.c */
void writeFloats(float*, int, char*);
void ReadMatrix(Matrix *, char*, wgttype threshold=0);
void readMemberships(char*, int, idxtype*, idxtype**, idxtype**,
				idxtype*);
idxtype* getNodesToComponentMap(Matrix*, int*, wgttype);
idxtype* compSizeDistribution(GraphType*, int*);
int isGraphConnected(GraphType*);
void WriteRMap(const char*, idxtype*, int);
void printHistogram(idxtype*, int, FILE*);
int readClustering(const char *, int *, int);
void ReadGraph(GraphType *, const char *, int *, int, int);
//void ReadTxtGraph(GraphType *, char *, int *, int);
void WritePartition(const char *, idxtype *, int, int);
void my_WritePartition(const char *, idxtype *, int, float);
void my_WritePartitionAddOne(const char *, idxtype *, int);
void my_WriteMappedPartition(const char *, idxtype *, idxtype *, int);
void WriteMeshPartition(const char *, int, int, idxtype *, int, idxtype *);
void WritePermutation(const char *, idxtype *, int);
int CheckGraph(GraphType *);
idxtype *ReadMesh(const char *, int *, int *, int *);
void WriteGraph(const char *, int, idxtype *, idxtype *);
void WriteTxtGraph(const char *, int, idxtype *, idxtype *);
void WriteMatrix(const char *, int, idxtype*, idxtype*, wgttype*);
void WriteGraphWithWts(const char *, int, idxtype*, idxtype*,
idxtype*);
void WriteMappedTxtGraphWithWts(const char *, int,
		idxtype*,idxtype*,idxtype*,idxtype*, int);


/* merge.c */
ListGraph* createClusterGraph(const idxtype*, int, const GraphType* );
void mergeBestClusters(ListGraph*, idxtype*, int, int);


/* bayesSimSearch.c */
Matrix* bss_degreeDiscountedWrapper(GraphType*, int, float,
float, float, float);
