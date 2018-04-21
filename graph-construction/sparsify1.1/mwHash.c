#include <metis.h>

int isPrime(idxtype n)
{
	int i,end;
	//rule out evens
	if (n <= 2) {return n == 2;}
	if (n % 2 == 0){return 0;}
	// rule out odds up to the square root of n
	for (i = 3, end = (idxtype)sqrt((float)n);i <= end; i += 2) 
	{
		if (n % i == 0) {return 0;}
	}
	return 1;
}

void generatePrimes(idxtype* tprimes, int numHashes)
{
  idxtype start=1000000001;
  int i;
  for(i=0;i<(numHashes*2);i++)
    {
      while (!isPrime(start)) 
	  {
	  	start++;
	  }
      tprimes[i]=start;
      start+=100000000;
    }
  //shuffle them
  for(i=0;i<1000;i++)
    {
      int val1  = ((int)rand()) % (numHashes*2); 
      int val2  = ((int)rand()) % (numHashes*2);
      idxtype temp = tprimes[val1];
      tprimes[val1]=tprimes[val2];
      tprimes[val2]=temp;
    }
  
}

idxtype hash(idxtype a, idxtype b, idxtype c, idxtype val)
{
	UInt32 val1 = (a*val);
	UInt32 val2 = val1+b;
	UInt32 val4 = val2 % ((UInt32) c);
	idxtype val3=(idxtype)val4;
	if(val3<=0)
	{
		printf("%d hashed to %d! a=%d, b=%d, c=%d\n", val, val3, a, b,
		c);
		  //cout <<a << "," << b << "," << c << "," << val << "," << val1 << "," << val2 << "," << val3 << endl;
	}
	return val3;
}

void generateRandoms(int numToGenerate, int* randoms )
{
	int i;
	int max = largePrime;
	int minRandomNumber = largePrime/100;
	for ( i=0; i<numToGenerate; i++ )
	{
		randoms[i] = (rand()) % max;
	//	if ( primes[i] == 0 )
		if ( randoms[i] < minRandomNumber )
			i--;
//		primes[i] = primes[i] % max;
	}
}



idxtype pickRandomElement(idxtype* set, int length)
{
	return set[RandomInRangeFast(length)];	
}

idxtype getMinHashKey(idxtype* set, int length, int *randoms, int
thisK, int numHashes)
{
	idxtype minHash=largePrime, minVal=-1;
	int i=0;
    for( i=0; i<length; i++)
    {
      idxtype x =  hash(randoms[thisK],randoms[thisK+numHashes],
	  						largePrime,set[i]);
      //  cout << "val=" << (*buffer)[i] << ",hash=" << x << endl;
    	if(x<minHash)
	    {
 	  		minHash=x;
			minVal=set[i];
		}
    }
    //cout << "Min hash is " << minHash << "," << minVal << endl;
    //cout << minHash << endl;
    //getchar();
    return minVal;
}

void getMinHashKeys( const idxtype * const set, int length, const
int * const randoms, int
totalNumHashes, int * const out , int * const mins, int numHashes)
{
	int t;
	int startHash = 0;
	for ( int i=0; i<length; i++ )
	{
		int feature = set[i];
		for ( int k=startHash, l=0; k<startHash+numHashes; k++,
		l++ )
		{
			t = hash(randoms[k], randoms[k+totalNumHashes],
						largePrime, feature);
			if ( t < mins[l] || i==0 )
			{
				out[l] = feature;
				mins[l] = t;
			}
		}
	}
}

Hashtable* buildHashtable(GraphType* graph, int numHashes, int*
randoms)
{
	int i;
	Hashtable* ht = (Hashtable*)malloc(sizeof(Hashtable));
	idxtype* base;
	
	ht->numNodes = graph->nvtxs;
	ht->numHashes = numHashes;
	ht->sortedNodeIds = idxmalloc(graph->nvtxs,
								"buildHashtable:ht->sortedNodeIds");
	ht->hashes = idxmalloc(graph->nvtxs*numHashes,
									"buildHashtable:ht->hashes");

	int *mins = (int*)malloc(sizeof(int)*numHashes);

	base = ht->hashes;
	idxtype* adjlist = graph->adjncy;
	for ( int j=0; j<graph->nvtxs; j++ )
	{
		int adjlistLength = graph->xadj[j+1]-graph->xadj[j];
		if ( adjlistLength == 1 )
		{
			for ( int i=0; i<numHashes; i++ )
				base[i] = adjlist[0];
		}
		else
		{
			getMinHashKeys(adjlist, adjlistLength, randoms,
			numHashes, base, mins, numHashes);
		}

		base += numHashes;
		adjlist += adjlistLength;
	}

	free(mins);
/*	for ( i=0; i<graph->nvtxs; i++ )
	{
		ht->sortedNodeIds[i] = i;
	}
*/
	return ht;
}

void freeHashtable(Hashtable *ht)
{
	if ( ht != NULL )
	{
		if ( ht->hashes != NULL )
			free(ht->hashes);
		if ( ht->sortedNodeIds != NULL )
			free(ht->sortedNodeIds);
		free(ht);
		ht = NULL;
	}
}

