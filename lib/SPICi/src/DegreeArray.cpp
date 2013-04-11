#include "DegreeArray.h"
#include <cstdlib>
using namespace std;

#define ROUND(x)((size_t)(x+0.5))

DegreeArray::DegreeArray(const double *degreeWeight,size_t V,double maxWeight)
:m_nVertex(V){
	size_t i,d;
	top = ROUND(maxWeight);

	indexmap = new list<size_t>::iterator[V];	// map from vertex index to round off degree bucket

	m_pDegree = new list<size_t>[top+1];

	for (i=0;i<V;i++)
	{
		d = ROUND(degreeWeight[i]);
		list<size_t>& l = m_pDegree[d];
		l.push_front(i);
		indexmap[i] = l.begin();
	}
}

DegreeArray::~DegreeArray(void)
{
	delete[] m_pDegree;
	delete[] indexmap;
}

size_t DegreeArray::getMax(){
	return m_pDegree[top].front();
}

size_t DegreeArray::extractMax()
{
	if (empty())
	{
		cerr << "cannot extract elements from an empty degree array."<<endl;
		exit(1);
	}

	size_t head = m_pDegree[top].front();
	m_pDegree[top].pop_front();
	m_nVertex--;

	while (m_pDegree[top].size() == 0 && top > 0 ){
		top --;
	}

	return head;
}

void DegreeArray::decrease(size_t inx,double oValue,double nValue)
{
	size_t oldvalue = ROUND(oValue),
		newvalue = ROUND(nValue);

	if (oldvalue == newvalue)return;

	list<size_t>::iterator ptr = indexmap[inx];
	m_pDegree[oldvalue].erase(ptr);
	m_pDegree[newvalue].push_front(inx);
	indexmap[inx] = m_pDegree[newvalue].begin();

	while (m_pDegree[top].size() == 0 && top > 0 ){
		top --;
	}
}

void DegreeArray::remove(size_t inx,double value)
{
	size_t v = ROUND(value);
	list<size_t>::iterator ptr = indexmap[inx];
	m_pDegree[v].erase(ptr);
	m_nVertex --;

	while (m_pDegree[top].size() == 0 && top > 0 ){
		top --;
	}
}
