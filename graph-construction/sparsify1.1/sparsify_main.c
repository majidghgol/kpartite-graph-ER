/*
 *
 * $Id: sparsify.c,v 1.7 2010-11-30 23:17:20 venu Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/

void print_help(const char * program_name)
{
	fprintf(stderr, "Usage: %s <GraphFile> ",program_name);
	fprintf(stderr, "[-b sparsifyType] [-m numHashes] ");
	fprintf(stderr, "[-e sparsifyExponent] [-r globalFraction] [-o output file]");
	fprintf(stderr, "\nsparsifyType: 0 (default) - fast,");
	fprintf(stderr, "approximate, local sparsify\n ");
	fprintf(stderr, "\t\t1 - slow, exact local sparsify\n ");
	fprintf(stderr, "\t\t2 - fast, approximate global sparsify\n"); 
	fprintf(stderr, "\t\t3 - slow, exact global sim sparsify\n");
	fprintf(stderr, "numHashes: number of hashes for approx.");
	fprintf(stderr, " similarity (default: 30) (only for approx. similarity)\n");
	fprintf(stderr, "sparsifyExponent (0 < e < 1): each node with degree d, d^e edges will be retained (default: 0.5) (only for local sparsification)\n");
	fprintf(stderr, "globalFraction (0 < r < 1): fraction of edges to retain in sparsified graph (default: 0.25) (only for global sparsification)\n");

}

main(int argc, char *argv[])
{

	GraphType graph;
	char filename[256],outputFile[256];
	int outputFileGiven = 0;
	int sparsifyType = 0, numHashes = 30; 
	float globalSamplingRatio = 0.25, degreeExponent = 0.5;
	timer TOTALTmr, METISTmr, IOTmr;
	int wgtflag = 0, addSelfLoop = 1, txtFormat = 0,
	minEdgeWeight = 0;

	if ( argc < 2 )
	{
		print_help(argv[0]);
		exit(0);
	}

	for (argv++; *argv != NULL; argv++)
	{
		if ((*argv)[0] == '-')
		{
			int temp;
			switch ((*argv)[1])
			{
			case 'o':
			case 'O':
				outputFileGiven = 1;
				strcpy( outputFile, *(++argv));
				break;
			case 'b':
			case 'B':
				sparsifyType = atoi(*(++argv));
				break;
			case 'm':
			case 'M':
				numHashes = atoi(*(++argv));
				break;
			case 'e':
			case 'E':
				degreeExponent = atof(*(++argv));
				break;
			case 'w':
			case 'W':
				minEdgeWeight = atoi(*(++argv));
				break;
			case 'r':
			case 'R':
				globalSamplingRatio = atof(*(++argv));
				break;
			}
		}
	    else
		{
	      strcpy(filename, *argv);
	    }
	}


	if ( !outputFileGiven )
	{
		fprintf(stderr, 
			"Please specify output file using -o option\n");
		exit(0);
	}

	cleartimer(TOTALTmr);
	cleartimer(METISTmr);
	cleartimer(IOTmr);

	starttimer(TOTALTmr);
	starttimer(IOTmr);

	ReadGraph(&graph, filename, &wgtflag, addSelfLoop, txtFormat);

	if (graph.nvtxs <= 0) {
	  printf("Empty graph. Nothing to do.\n");
	  exit(0);
	}
	stoptimer(IOTmr);

	printf("Graph Information ---------------------------------------------------\n");
	printf("  Name: %s, #Vertices: %d, #Edges: %d \n", filename,
	graph.nvtxs, graph.nedges/2 );
	fflush(stdout);

	starttimer(METISTmr);

	GraphType *inputg;
	inputg = &graph;

	idxtype *new_xadj, *new_adjncy, *new_adjwgt;
	if ( sparsifyType == 0 ) // Estimate similarity using minhash
	//and then retain similar edges.
	{
		GraphType *rg;
		rg = getHashSparsifiedGraph(inputg, numHashes,
		degreeExponent);
	//	new_xadj = rg->xadj;
	//	new_adjncy = rg->adjncy;

		new_xadj = rg->xadj;
		new_adjncy = rg->adjncy;
		usual_symmetrize_directed(graph.nvtxs,rg->xadj[graph.nvtxs],
		rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy,
		&new_adjwgt, 0); 
	}
	else if (sparsifyType == 2 )
	{
		GraphType *rg;
		rg =
		getGlobalHashSparsifiedGraph(inputg,globalSamplingRatio,
		numHashes);
	//	new_xadj = rg->xadj;
	//	new_adjncy = rg->adjncy;
	
		usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs],
		rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy,
		&new_adjwgt, 0); 
		
	}
	else if ( sparsifyType == 1 )
	{
		GraphType *rg;
		rg = getExactSimSparsifiedGraph(inputg, degreeExponent);

		usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs],
		rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy,
		&new_adjwgt, 0); 
	}
	else if ( sparsifyType == 3 )
	{
		GraphType *rg;
		rg = getExactGlobalSimSparsifiedGraph(inputg,
		globalSamplingRatio);

		usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs],
		rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy,
		&new_adjwgt, 0); 
	}

	stoptimer(METISTmr);

	printf("No. of edges in sparsified graph:%d\n",
			(new_xadj[inputg->nvtxs]/2));

	starttimer(IOTmr);
	if ( 0 /*sparsifyType == 7*/ )
		WriteGraphWithWts(outputFile, inputg->nvtxs, new_xadj,
		new_adjncy, new_adjwgt);
	else
		WriteGraph(outputFile, inputg->nvtxs, new_xadj, new_adjncy);
	stoptimer(IOTmr);
	stoptimer(TOTALTmr);

	printf("\nTiming Information --------------------------------------------------\n");
	printf(" I/O:\t\t\t\t\t\t\t%7.3f\n", gettimer(IOTmr));
	printf(" Time for sparsifying:\t\t\t\t\t%7.3f\n",
					gettimer(METISTmr));
	printf(" Total:\t\t\t\t\t\t\t%7.3f\n", gettimer(TOTALTmr));
	printf("**********************************************************************\n");


	

}
