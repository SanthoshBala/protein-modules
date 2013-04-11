/************************************************************************/
/* Implemented according to <<The introduction to algorithm>>
	operator > and < of Key type must be overloaded!
	Peng Jiang		Email: flute.peng.jiang@gmail.com					*/
/************************************************************************/
#pragma once
#include <iostream>
#include "FastGraphCluster.h"
using namespace std;

// first version for programming convenience
typedef HeapNode KEYTYPE;

class FiboNode{
public:
	FiboNode(KEYTYPE key,FiboNode *parent = NULL,FiboNode *child = NULL,
		FiboNode *left = NULL ,FiboNode *right = NULL,
		size_t degree = 0,bool mark = false);
	~FiboNode(){}

	KEYTYPE key;	// key value
	FiboNode *parent,*child;	//	child points to any one of its child
	FiboNode *left,*right;	//	link this node into the double linked list of Sibling
	size_t degree;	// the number of children in the child list
	bool mark;	// whether lost a child since the last time x was made the child of another node
};

// in order to get an inverted version. simply change the operate >

class FibonacciHeap
{
public:
	FibonacciHeap(void);
	~FibonacciHeap(void);

	FiboNode* insert(KEYTYPE x);	// insert a new node and return the relevant pointer
	void insert(FiboNode *node);	// make the node yourself first,don't forget to clear the node

	KEYTYPE getMin();	// get the minimum key value
	FiboNode* extractMin();	// extract the min heap

	void decrease(FiboNode *node,KEYTYPE newKey);	// decrease the key
//	void remove(FiboNode *x);
	bool empty(){return minheap == NULL;}
	void combine(FibonacciHeap* rheap);	// combine rheap with ourself
private:
	FiboNode *minheap;
	bool *mark;
	size_t m_nNode;	// the current number of nodes
	void consolidate();
	void cut(FiboNode *x,FiboNode *y);
	void cascadingCut(FiboNode *y);
	void fibHeapLink(FiboNode* y,FiboNode* x);
	void dispose(FiboNode *y);	// dispose the subtree rooted at y
};
