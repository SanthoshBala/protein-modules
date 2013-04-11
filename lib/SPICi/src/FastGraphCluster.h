#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
using namespace std;

class FiboNode;
class BinaryHeap;
class DegreeArray;

// node used in heuristic vertex adding
class HeapNode{
public:
	double wsum;	// how many edges weight are connected to current module
	size_t index;	// index of the vertex

	HeapNode( double wsum = 0 , size_t index = 0 ):wsum(wsum),index(index){}

	// the order here are totally reversed. we will get a max heap instead
	bool operator >(const HeapNode &right);
	bool operator >=(const HeapNode &right);
	bool operator <(const HeapNode &right);
	bool operator <=(const HeapNode &right);
};

class FastGraphCluster
{
public:
	FastGraphCluster(string file,double density,size_t lowersize,double lowerincrease,char mode);
	~FastGraphCluster(void);

	// a fast neighbor searching method: depth: neighbor radius to grasp
	size_t fastCluster(string output);

private:
	// starting from a local seed posize_t and heuristically get the local dense module
	// index: starting posize_t
	double expand(size_t index,vector<size_t> &result);
	
	// return a weight decided by edge and node neighbors
	double getWeight(double edgeweight,double vertexweight);

	size_t *neighborCnt;
	double *neighborWeightCnt;	// total sum of the weight for its neighbor, seed value heuristic

	string *vertexName;
	FiboNode **m_pHeapNode;	// we will use Fibonacci heap
	size_t m_nVertex,m_nEdge;
	string file;
	bool *changedFlag;

	DegreeArray *seedArray;	// extract min + decrease key

	// parameters for the algorithm searching
	size_t m_nLowerSize;	// lower size of the modules
	double m_dLowerDensity,m_dLowerIncrease;	// lower density of the modules

	// map every vertex name to its index
	map<string,size_t> vertexNameMap;

	// Graph representation, choose one according to mode
	char graphMode;
	// linked list graph
	size_t **m_pNeighbor;	// neighbor index and its
	float **m_pWeight;	// weight value

	// adjacent matrix
	vector<float*> m_pMatrix;
};
