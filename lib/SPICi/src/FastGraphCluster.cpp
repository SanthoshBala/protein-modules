#include "FastGraphCluster.h"
#include "FibonacciHeap.h"
#include "DegreeArray.h"
#define UNEXPLORED -1
#define SECONDSEED
#define MAX_NAME_LENGTH 1000

FastGraphCluster::FastGraphCluster(string file,double density,size_t lowersize,double lowerincrease,char mode)
:	file(file),m_nLowerSize(lowersize),m_dLowerDensity(density),m_dLowerIncrease(lowerincrease),
	m_pNeighbor(NULL),m_pWeight(NULL),graphMode(mode)
{
	// pass through the file for the first time
	size_t i,j;
	string sa,sb;
	float weight;
	ifstream fin(file.c_str());
	size_t ia,ib;

	vector<size_t> vecneighborcnt;

	// temporary list for id and weights, may not be used
	list<size_t> *idalist,*idblist;
	list<float> *weightlist;
	list<string> *namelist;

	float *array;

	idalist = new list<size_t>;
	idblist = new list<size_t>;
	weightlist = new list<float>;
	namelist = new list<string>;

	m_nEdge = m_nVertex = 0;

	map<string,size_t>::iterator miter;
	
	char tempchar[MAX_NAME_LENGTH];
	float weightsum = 0;

	while(!fin.eof())
	{
		fin.getline(tempchar,MAX_NAME_LENGTH,'\t');
		sa = tempchar;

		fin.getline(tempchar,MAX_NAME_LENGTH,'\t');
		sb = tempchar;

		fin >> weight;
		fin.getline(tempchar,MAX_NAME_LENGTH,'\n');	// jump the last space

		if (sa.length() ==0 || sb.length()==0 || fin.eof())break;	// reach the end of file

		if ( sa == sb ){
			cerr << "Warning: cannot include self interaction:\t"<< sa << '\t' << sb << endl;
			continue;
		}

		if (weight <=0 || weight > 1){
			cerr << "Warning: illegal interaction weight:\t" << sa << '\t' << sb << '\t' << weight << endl;
			continue;
		}

		weightsum += weight;
		m_nEdge++;

		// get the vertex index map first
		miter = vertexNameMap.find(sa);
		if( miter == vertexNameMap.end())
		{
			ia = vertexNameMap[sa] = m_nVertex;
			m_nVertex++;
			vecneighborcnt.push_back(0);
			namelist->push_back(sa);

		}else{
			ia = miter->second;
		}

		miter = vertexNameMap.find(sb);
		if( miter == vertexNameMap.end())
		{
			ib = vertexNameMap[sb] = m_nVertex;
			m_nVertex++;
			vecneighborcnt.push_back(0);
			namelist->push_back(sb);
		}else{
			ib = miter->second;
		}

		vecneighborcnt[ia]++;
		vecneighborcnt[ib]++;

		if (mode == '0')
		{	// sparse graph
			idalist->push_back(ia);
			idblist->push_back(ib);
			weightlist->push_back(weight);

		}else if (mode == '1')
		{	// nearly complete graph
			if ( ia < ib ){swap(ia,ib);}	// ia > ib

			// add up matrix size if necessary
			while ( m_pMatrix.size() <= ia )
			{
				array = new float[m_pMatrix.size()];
				for (i=0;i<m_pMatrix.size();i++){array[i] = 0;}
				m_pMatrix.push_back(array);
			}

			// write down lower triangle
			m_pMatrix[ia][ib] = weight;
		}
		//else{} the rest mode, wait for second read
	}

	weightsum /= m_nEdge;

	if (weightsum < m_dLowerDensity )
	{
		cerr << "Warning: low average edge confidence " << weightsum << endl;
	}

	// transfer the vertex name
	vertexName = new string[m_nVertex];

	list<string>::iterator siter = namelist->begin();
	for (i=0;i<m_nVertex;i++,siter++){
		vertexName[i] = *siter;
	}
	delete namelist;

	// these statistics should be useful for any mode
	neighborCnt = new size_t[m_nVertex];
	neighborWeightCnt = new double[m_nVertex];
	changedFlag = new bool[m_nVertex];

	for (i=0;i<m_nVertex;i++)
	{
		neighborCnt[i] = 0;
		neighborWeightCnt[i] = 0;
		changedFlag[i] = false;
	}

	if (mode!='1')
	{	// not nearly complete graph, use linked list
		m_pNeighbor = new size_t*[m_nVertex];
		m_pWeight = new float*[m_nVertex];

		for (i=0;i<m_nVertex;i++){
			m_pNeighbor[i] = new size_t[vecneighborcnt[i]];
			m_pWeight[i] = new float[vecneighborcnt[i]];
		}
	}

	fin.close();

	ifstream secfin;

	if (mode=='2')
	{	// double read for non-sparse, non-nearly complete graph
		secfin.open(file.c_str());
	}

	list<size_t>::iterator aiter,biter;
	list<float>::iterator witer;
	aiter = idalist->begin();
	biter = idblist->begin();
	witer = weightlist->begin();

	if (mode == '1'){
		for (i=0;i<m_nVertex;i++)
		{
			for (j=0;j<i;j++)
			{
				weight = m_pMatrix[i][j];
				if ( weight > 0 ){
					neighborCnt[i]++;
					neighborCnt[j]++;

					neighborWeightCnt[i] += weight;
					neighborWeightCnt[j] += weight;
				}
			}
		}

	}else{
		for (i=0;i<m_nEdge;)
		{
			if (mode =='2')
			{
				secfin >> sa >> sb >> weight;
				if ( sa == sb )continue;	// jump all self size_teractions
				ia = vertexNameMap[sa];
				ib = vertexNameMap[sb];
			}else { // mode == '0'
				ia = *aiter;
				ib = *biter;
				weight = *witer;
			}

			m_pNeighbor[ia][neighborCnt[ia]] = ib;
			m_pNeighbor[ib][neighborCnt[ib]] = ia;

			m_pWeight[ia][neighborCnt[ia]] = weight;
			m_pWeight[ib][neighborCnt[ib]] = weight;

			neighborCnt[ia]++;
			neighborCnt[ib]++;

			neighborWeightCnt[ia] += weight;
			neighborWeightCnt[ib] += weight;
			i++;

			if (mode =='0')
			{
				aiter++;
				biter++;
				witer++;
			}
		}

		if (mode == '2')
		{	// double read
			secfin.close();
		}
	}

	delete idalist;
	delete idblist;
	delete weightlist;

	double maxWeightDegree = 0;

	// open the space for Fibonacci heap nodes
	m_pHeapNode = new FiboNode*[m_nVertex];
	for (i=0;i<m_nVertex;i++){
		m_pHeapNode[i] = new FiboNode(HeapNode(UNEXPLORED,i));
		if( neighborWeightCnt[i] > maxWeightDegree ) maxWeightDegree = neighborWeightCnt[i];
	}

	// prepare seeds list
	seedArray = new DegreeArray(neighborWeightCnt,m_nVertex,maxWeightDegree);
}

FastGraphCluster::~FastGraphCluster(void)
{
	for (size_t i=0;i<m_nVertex;i++){
		if (m_pHeapNode[i]!=NULL){
			delete m_pHeapNode[i];
		}
	}

	if ( m_pWeight != NULL )
	{
		for (size_t i=0;i<m_nVertex;i++){
			delete[] m_pWeight[i];
			delete[] m_pNeighbor[i];
		}

		delete[] m_pWeight;
		delete[] m_pNeighbor;
	}

	delete[] changedFlag;
	delete[] m_pHeapNode;
	delete[] vertexName;
	delete[] neighborCnt;
	delete[] neighborWeightCnt;
	delete seedArray;

	// free potential matrix graph, even if it might be empty
	vector<float*>::iterator iter = m_pMatrix.begin();
	for (;iter!=m_pMatrix.end();iter++){
		delete[] *iter;
	}
}

size_t FastGraphCluster::fastCluster(string output)
{
	size_t i=0,size,total = 0;
	size_t clusterCnt = 0,maxCluster = 0;

	vector<size_t> result;
	vector<size_t>::iterator iter;
	ofstream fout(output.c_str());
	fout.clear();

	while (!seedArray->empty() && seedArray->top >= m_nLowerSize)
	{
		i = seedArray->getMax();

		expand(i,result);

		size = (size_t)result.size();
		total += size;

		if ( size >= m_nLowerSize ){
			for ( i = 0 ; i < size ; i++ ){
				fout<< vertexName[result[i]] <<(i==size-1?'\n':'\t');
			}

			clusterCnt++;
			maxCluster += size;
		}
	}

	fout.close();

	if(clusterCnt > 0){
	  cout <<"average cluster size: "<< (float)maxCluster/clusterCnt <<endl;
	}

	return clusterCnt;
}

double FastGraphCluster::getWeight(double edgeWeight,double vertexWeight){
  return m_nVertex*(size_t)(5*edgeWeight) + vertexWeight;
}

double FastGraphCluster::expand(size_t index,vector<size_t> &result)
{
	size_t i,j,nVertex = 0,*neighborlist,neighborcnt;
	double density = 0,increase,w,totalWeight = 0,oldWeight;
	float *weightlist;

	FibonacciHeap heap;	// local expanding heap
	vector<size_t> changed;

	FiboNode *ptr = m_pHeapNode[index];
	ptr->key.wsum = 0;
	heap.insert(ptr);
	result.clear();

#ifdef SECONDSEED
	// start second level heuristic value selection
	size_t maxindex;
	float maxWeight;
	neighborcnt = neighborCnt[index];
	maxWeight = 0;
	maxindex = 0;

	if ( graphMode != '1' )
	{	// not complete graph
		neighborlist = m_pNeighbor[index];
		weightlist = m_pWeight[index];

		for (i=0;i<neighborcnt;i++)
		{
			j = neighborlist[i];

			if (m_pHeapNode[j] == NULL){	// it is searched before
				neighborcnt ++;	// recover the neighbor count
				continue;
			}

			//w = 1 + sqrt(weightlist[i]*neighborWeightCnt[j]);	// edge weight + vertex weight sum
			//w = weightlist[i];
			w = getWeight(weightlist[i],neighborWeightCnt[j]);

			if (w > maxWeight ) {
				maxindex = i;
				maxWeight = (float)w;
			}
		}

		// now maxindex is what we want here
		// we change this edge to be really big, so for the expanding procedure, we definitely select this edge first
		weightlist[maxindex] += maxWeight;

	}else{
		for (i=0;i<m_nVertex;i++)
		{
			if ( i == index || m_pHeapNode[i] == NULL){
				continue;
			}

			if (i>index){
				w = m_pMatrix[i][index];
			}else{
				w = m_pMatrix[index][i];
			}

			if (w==0){
				continue;
			}

			w = getWeight(w,neighborWeightCnt[i]);

			if (w > maxWeight){
				maxindex = i;
				maxWeight = (float)w;
			}
		}

		i = maxindex;
		j = index;
		if (i<j)swap(i,j);
		m_pMatrix[i][j] += maxWeight;
	}

#endif

	while (!heap.empty())
	{
		ptr = heap.extractMin();
		i = ptr->key.index;
		w = ptr->key.wsum;

		// increasing the size
		totalWeight += w;
		nVertex ++;

		if (nVertex > 1){
		        density = 2.0*totalWeight/(nVertex*(nVertex-1));
			increase = w/(density*nVertex);

			if ( density < m_dLowerDensity || increase < m_dLowerIncrease ){
				// recover the old density and exit
				nVertex--;
				totalWeight -= w;
				density = 2.0*totalWeight/(nVertex*(nVertex-1));
				break;	// exceeding the lowest density
			}
		}

		// free relevant fields
		delete ptr;
		m_pHeapNode[i] = NULL;

		// remove merged vertex from seed list
		seedArray->remove(i,neighborWeightCnt[i]);

		result.push_back(i);	// write back result

		if ( graphMode != '1' )
		{
			neighborcnt = neighborCnt[i];
			neighborlist = m_pNeighbor[i];
			weightlist = m_pWeight[i];

			// expanding vertex neighbor and updated relevant values
			for (i=0;i<neighborcnt;i++)
			{
				j = neighborlist[i];
				w = weightlist[i];
				ptr = m_pHeapNode[j];

				if (ptr!=NULL)
				{	// not expand
					if (ptr->key.wsum == UNEXPLORED)
					{	// a new expanded vertex
						ptr->key.wsum = w;
						heap.insert(ptr);
					}else{	// previously expanded.
						heap.decrease(ptr,HeapNode( ptr->key.wsum + w , j ) );	// increasing
					}

					// decrease its weight degree no matter whether it is deleted.
					oldWeight = neighborWeightCnt[j];

#ifdef SECONDSEED
					if ( ( nVertex == 1 ) && (i == maxindex) ){
						// since we changed the edge weight to "cheat", we need to be careful when deducing the first level seed value
						seedArray->decrease(j,oldWeight,oldWeight + maxWeight - w);
						neighborWeightCnt[j] -= (w-maxWeight);
					}else{
						seedArray->decrease(j,oldWeight,oldWeight - w);
						neighborWeightCnt[j] -= w;
					}
#else
					seedArray->decrease(j,oldWeight,oldWeight - w);
					neighborWeightCnt[j] -= w;
#endif

					neighborCnt[j]-- ;

					if (changedFlag[j] == false){
						changedFlag[j] = true;
						changed.push_back(j);
					}

				}else{
					neighborcnt++;	// put back deduction on neighborCnt
				}
			}
		}else{	// mode == '1' complete graph
			// i is the new input node
			// expanding vertex neighbor and updated relevant values
			for (j=0;j<m_nVertex;j++)
			{
				ptr = m_pHeapNode[j];
				if ( i == j || ptr == NULL ){continue;}

				if (i<j) w = m_pMatrix[j][i];
				else w = m_pMatrix[i][j];

				if (w==0){continue;}

				if (ptr->key.wsum == UNEXPLORED)
				{	// a new expanded vertex
					ptr->key.wsum = w;
					heap.insert(ptr);
				}else{	// previously expanded.
					heap.decrease(ptr,HeapNode( ptr->key.wsum + w , j ) );	// increasing
				}

				// decrease its weight degree no matter whether it is deleted.
				oldWeight = neighborWeightCnt[j];

#ifdef SECONDSEED
				if ( ( nVertex == 1 ) && (j == maxindex) ){
					// since we changed the edge weight to "cheat", we need to be careful when deducing the first level seed value
					seedArray->decrease(j,oldWeight,oldWeight + maxWeight - w);
					neighborWeightCnt[j] -= (w-maxWeight);
				}else{
					seedArray->decrease(j,oldWeight,oldWeight - w);
					neighborWeightCnt[j] -= w;
				}
#else
				seedArray->decrease(j,oldWeight,oldWeight - w);
				neighborWeightCnt[j] -= w;
#endif

				neighborCnt[j]-- ;

				if (changedFlag[j] == false){
					changedFlag[j] = true;
					changed.push_back(j);
				}
			}
		}

#ifdef SECONDSEED
		if (nVertex == 1 ){
			totalWeight -= maxWeight;	// deduce back the extra weight we put for seeding pair
		}
#endif
	}

	// for all changed posize_ts, future is a new beginning for them.
	for ( vector<size_t>::iterator iter = changed.begin() ; iter != changed.end() ; iter++ )
	{
		i = *iter;
		changedFlag[i] = false;
		if ( m_pHeapNode[i] != NULL ){
			m_pHeapNode[i]->key.wsum = UNEXPLORED;
		}
	}

	return density;
}

bool HeapNode::operator >(const HeapNode &right)
{
	return wsum<right.wsum;
}

bool HeapNode::operator <(const HeapNode &right)
{
	return wsum>right.wsum;
}

bool HeapNode::operator <=(const HeapNode &right)
{
	return wsum >= right.wsum;
}

bool HeapNode::operator >=(const HeapNode &right)
{
	return wsum <= right.wsum;
}
