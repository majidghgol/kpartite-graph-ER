/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h,v 1.12 2011-04-20 20:02:20 venu Exp $
 */

/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT

#include <stdint.h>
const int largePrime = 2147483647;

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype;
#else
typedef short idxtype;
#endif

/* Venu : my addition below */
/* Undefine WGTTYPE_FLOAT in order to use double as the wgttype
 * */
#define WGTTYPE_FLOAT
/* wgttype is used to store weights in the matrices*/
#ifdef WGTTYPE_FLOAT
typedef float wgttype;
#else
typedef double wgttype;
#endif

#define MAXIDX	(1<<8*sizeof(idxtype)-2)


/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;


/*************************************************************************
* The following data structure will hold a node of a doubly-linked list.
**************************************************************************/
struct ListNodeType {
  int id;                       	/* The id value of the node */
  struct ListNodeType *prev, *next;     /* It's a doubly-linked list */
};

typedef struct ListNodeType ListNodeType;



/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct PQueueType {
  int type;                     /* The type of the representation used */
  int nnodes;
  int maxnodes;
  int mustfree;

  /* Linear array version of the data structures */
  int pgainspan, ngainspan;     /* plus and negative gain span */
  int maxgain;
  ListNodeType *nodes;
  ListNodeType **buckets;

  /* Heap version of the data structure */
  KeyValueType *heap;
  idxtype *locator;
};

typedef struct PQueueType PQueueType;


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edegreedef {
  idxtype pid;
  idxtype ed;
};
typedef struct edegreedef EDegreeType;


/*************************************************************************
* The following data structure stores an edge for vol
**************************************************************************/
struct vedegreedef {
  idxtype pid;
  idxtype ed, ned;
  idxtype gv;
};
typedef struct vedegreedef VEDegreeType;

/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/*************************************************************************
* This data structure holds various working space data
**************************************************************************/
struct workspacedef {
  idxtype *core;			/* Where pairs, indices, and degrees are coming from */
  int maxcore, ccore;

  EDegreeType *edegrees;
  VEDegreeType *vedegrees;
  int cdegree;

  idxtype *auxcore;			/* This points to the memory of the edegrees */

  idxtype *pmat;			/* An array of k^2 used for eliminating domain 
                                           connectivity in k-way refinement */
};

typedef struct workspacedef WorkSpaceType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 int id, ed;            	/* ID/ED of nodes */
 int ndegrees;          	/* The number of different ext-degrees */
 EDegreeType *edegrees;     	/* List of edges */
};

typedef struct rinfodef RInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* vol-based partition
**************************************************************************/
struct vrinfodef {
 int id, ed, nid;            	/* ID/ED of nodes */
 int gv;            		/* IV/EV of nodes */
 int ndegrees;          	/* The number of different ext-degrees */
 VEDegreeType *edegrees;     	/* List of edges */
};

typedef struct vrinfodef VRInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 idxtype edegrees[2];  
};

typedef struct nrinfodef NRInfoType;

/* Venu: my addition */
/****************
* This data structures holds a (square) sparse matrix
*/
struct matrixdef{
	int nvtxs, nnz; /* nnz stands for number of non-zero
	entries*/
	idxtype* xadj; /* xadj[i+1]-xadj[i] gives the number of
	non-zero entries in column i */
	idxtype* adjncy;
	wgttype* adjwgt; /* array that stores weights of the
	adjacency lists */
	wgttype* adjwgtsum; /* sum of adjacency weights of each node,
	or the sum of a column, basically. assigned only during the
	expandAndInflate phase */
	wgttype* maxwgt; /* max wgt of each column, assigned only
	during the expandAndInflate phase */
	idxtype* attractors; /* the row with the maximum weight in a
	column, and which has weight > 0.50 */

	idxtype* rmap; /* This is for mis_coarsen */
	int currentSize;
	int sizeIncrement;
};

typedef struct matrixdef Matrix;

/* hashStoreType defines granularity with which we will add new
 * hashes. If we want to add more hashes at each increment, use
 * long; otherwise use int. */
// typedef long hashStoreType;
/* stores the random vectors with which dot products will be
 * computed to estimate cosine similarity between vectors. */
struct coshashvectors
{
	int nDimensions;
	uint8_t* hashVectors; // the first nHashes/8 bytes give the bits
	// of the first dimension in each of the nHashes vectors. So,
	// the hash vectors are stored with a stride of nHashes bits. 
	
	int nHashes; 
	// must be a multiple of 8.
};

typedef struct coshashvectors CosHashVectors;

class BoolVector
{
	public:
		uint8_t *a;
		int length; // must be a multiple of 8.

		BoolVector();
		~BoolVector();
};

class CosineSketches
{
	public:

	int nPoints;
	BoolVector *extraSketches;
	uint8_t *minSketches; 
	int minSketchSize; // must be a multiple of 8.

		CosineSketches(int, int);
		void buildMinSketches(const Matrix*, const
		CosHashVectors);
		~CosineSketches();
};

struct coshashsketches
{
	int nPoints;
	int nHashes; // must be a multiple of 8.
	void *sketches; // first nHashes/8 bytes contain sketch of
	//first point, and so on.
};

typedef struct coshashsketches CosHashSketches;

struct beta
{
	float alpha, beta;
};
typedef struct beta Beta;

class PwUnif{
	public:
	float p1, p2;
	float lb, ub;
	float changePt;
};

const int LnGamma_nudger = 200;
class LnGamma{
	float init;
	int max;
	float *lg;
	int lgIndex;

	public:
		LnGamma(float, int);
		void fillLnGammas(int);
		float lnGamma(int, int fill=1);
		~LnGamma();
};

const int BetaCdf_nudger = 10000;
class BetaCdf{
	float init_a, init_b, x;
	float min_a_b;
	int maxHashes;
	LnGamma *lgA, *lgB, *lgAB;
	float *betaIncs; 
	int biIndex_a, biIndex_b;
	float lnx, lnminusx;

	public:
		BetaCdf(float , float , float , int );
		~BetaCdf();
		void fillLnGammas(int);
		void fillBetaIncs(int, int);
		float lnGamma(int, int fill=0);
		float betaInc(int, int, int fill=0);
		void clear();
		
};

struct bssOptionsDef{
	float cosThreshold;
	float angleThreshold;
	float cosDelta;
	float eps1;
	float eps2;
	int numSamples;
};

typedef struct bssOptionsDef BssOptions;

class BayesSimSearch{

	int nPoints;
	int nDimensions;

	idxtype *xadj, *xadj_ends, *adjncy, *invIndexAdjncy,
	*invIndex_starts, *invIndex_ends;
	wgttype *adjwgt, *invIndexWgts;

	long numPruned;
	long numFalsePositives;
	long numDotProducts;
	int minMatchesCalcs;

	Matrix* simPairs();
	void addPairsForVertex(const int, Matrix*, uint8_t*);
	float posteriorCdfUsingCaches(int nMatches, int nTrials,
	float x);
	int computeMinMatches(int);
	void computeMinMatches();

	public:

	Matrix *fvs, *reordered ;
	Beta betaPrior;
	PwUnif puPrior;
	int isBetaPrior;
	
	timer priorTimer, skTimer, simTimer;
	BssOptions opts;

	wgttype* maxInWeights, *maxOutWeights;
	idxtype* inDegrees;
	CosHashVectors chv;
	CosineSketches *cs;

	BetaCdf *oneHalfBetaCache, *changePtBetaCache;
	float *thresholdCdfCache;
	int *minMatches, numMinMatches;

	BayesSimSearch(Matrix *fvs, BssOptions opts, CosHashVectors
	chv):fvs(fvs), opts(opts), chv(chv){}

	Matrix* run();
};


struct hashtabledef
{
	int numHashes;
	int numNodes;

	// sortedNodeIds is the list of nodeIds sorted according to
	// their fingerprints.
	idxtype *sortedNodeIds; 

	idxtype *hashes;
};

typedef struct hashtabledef Hashtable;

typedef unsigned int UInt32;
//typedef int UInt32;

struct threadData {
	int nvtxs;
	idxtype *xadj, *adjncy, *adjwgt, *r_xadj, *r_adjncy,
	*r_adjwgt, *outDegrees, *inDegrees;
	float threshold, error;
	float alpha, beta;
	idxtype **ret_xadj, **ret_adjncy, **ret_adjwgt;
	Matrix **ret;
};

/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  idxtype *gdata, *rdata;	/* Memory pools for graph and refinement data.
                                   This is where memory is allocated and used
                                   the rest of the fields in this structure */

  int nvtxs, nedges;		/* The # of vertices and edges in the graph */
  idxtype *xadj;		/* Pointers to the locally stored vertices */
  idxtype *vwgt;		/* Vertex weights */
  idxtype *vsize;		/* Vertex sizes for min-volume formulation */
  idxtype *adjncy;		/* Array that stores the adjacency lists of nvtxs */
  idxtype *adjwgt;		/* Array that stores the weights of the adjacency lists */

  idxtype *adjwgtsum;		/* The sum of the adjacency weight of each vertex */

  idxtype *label;

  idxtype *cmap;
  
  /* Venu: my addition
  * Refine map: maps vertices in coarse graph to vertices in
   * refined graph. Two maps needed as 2 vertices are mapped to
   * one vertex in the coarse graph */
  idxtype *rmap1;
  idxtype *rmap2;
  idxtype *numDescendants;
  /* indicates if this graph is the original
  (i.e. the most refined) graph*/
  int isOrgGraph; 
  int isDirected; 
  wgttype *pagerank;

  /* Partition parameters */
  int mincut, minvol;
  idxtype *where, *pwgts;
  int nbnd;
  idxtype *bndptr, *bndind;

  /* Bisection refinement parameters */
  idxtype *id, *ed;

  /* K-way refinement parameters */
  RInfoType *rinfo;

  /* K-way volume refinement parameters */
  VRInfoType *vrinfo;

  /* Node refinement information */
  NRInfoType *nrinfo;


  /* Additional info needed by the MOC routines */
  int ncon;			/* The # of constrains */ 
  float *nvwgt;			/* Normalized vertex weights */
  float *npwgts;		/* The normalized partition weights */

  struct graphdef *coarser, *finer;
};

typedef struct graphdef GraphType;

struct floatintDef{
	float f;
	int i;
};
typedef struct floatintDef FloatInt;

struct intlistdef{
	idxtype* l;
	int length;
	int allocSize;
	int increment;
};
typedef struct intlistdef ListInt;

struct wgtlistdef{
	wgttype* l;
	int length;
	int allocSize;
	int increment;
};
typedef struct wgtlistdef ListWgt;

struct listgraphdef{
	int nvtxs;
	int nedges;
	ListInt* adjLists;
	ListWgt* wgtLists;
	wgttype* vols; // volumes, i.e. sum of edge weights
};
typedef struct listgraphdef ListGraph;


/*************************************************************************
* The following structure stores information used by Metis
**************************************************************************/
struct controldef {
  int CoarsenTo;		/* The # of vertices in the coarsest graph */
  int dbglvl;			/* Controls the debuging output of the program */
  int CType;			/* The type of coarsening */
  int IType;			/* The type of initial partitioning */
  int RType;			/* The type of refinement */
  int maxvwgt;			/* The maximum allowed weight for a vertex */
  float nmaxvwgt;		/* The maximum allowed weight for a vertex for each constrain */
  int optype;			/* Type of operation */
  int pfactor;			/* .1*prunning factor */
  int nseps;			/* The number of separators to be found during multiple bisections */
  int oflags;

  WorkSpaceType wspace;		/* Work Space Informations */

  /* Various Timers */
  timer TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, UncoarsenTmr, 
        SepTmr, RefTmr, ProjectTmr, SplitTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;

};

typedef struct controldef CtrlType;

struct optionsdef{
	int coarsenTo;
	float gamma;
	int iter_per_level;
	int num_last_iter;
	int ncutify;
	int hubRemoval;
	float hubPct;
	int mis_coarsenType;
	int matchType;
	int exact;
	int k;
	float dpr_threshold;
	int transformAdj;
	float penalty_power; // used in transformAdj
};

typedef struct optionsdef Options;

struct pagerankOptionsDef{
	wgttype alpha; // random jump probability.
	int max_iters;
	wgttype convergeThreshold;
};

typedef struct pagerankOptionsDef PageRankOptions;

struct dirToUndirOptionsDef{
	int conversionMethod;
	wgttype threshold, threshold2;
	int blockSize;
	wgttype bcWeight;
	wgttype prunePercent;
	wgttype scale;
	int invDegreeType;
};

typedef struct dirToUndirOptionsDef DirToUndirOptions;

/*************************************************************************
* The following data structure stores max-partition weight info for 
* Vertical MOC k-way refinement
**************************************************************************/
struct vpwgtdef {
  float max[2][MAXNCON];
  int imax[2][MAXNCON];
};

typedef struct vpwgtdef VPInfoType;




