#pragma once
#include <iostream>
#include <list>
using namespace std;

class DegreeArray
{
public:
	DegreeArray(const double *degreeWeight,size_t V,double maxWeight);
	~DegreeArray(void);
	size_t extractMax();	// return the vertex index
	size_t getMax();
	bool empty(){return m_nVertex ==0;}
	void decrease(size_t inx,double oValue,double nValue);
	void remove(size_t inx,double value);
	size_t m_nVertex,top;

private:
	list<size_t>::iterator *indexmap;
	list<size_t> *m_pDegree;	// round off the vertex weight to integer
};
