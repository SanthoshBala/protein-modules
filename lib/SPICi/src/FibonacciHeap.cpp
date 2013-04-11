#include "FibonacciHeap.h"
#include <cstdlib>
#include <cmath>
using namespace std;

#define EXDEBUG

FibonacciHeap::FibonacciHeap(void)
{
	minheap = NULL;
	m_nNode = 0;
}

FibonacciHeap::~FibonacciHeap(void)
{
	//dispose(minheap);
}

FiboNode::FiboNode(KEYTYPE key,
	FiboNode *parent /* = NULL */,FiboNode *child /* = NULL */,
	FiboNode *left /* = NULL  */,FiboNode *right /* = NULL */,
	size_t degree /* = 0 */,bool mark /* = false */
):key(key),parent(parent),child(child),left(left),right(right),degree(degree),mark(mark){}


// no consolidations here
FiboNode* FibonacciHeap::insert(KEYTYPE x)
{
	FiboNode *node = new FiboNode(x);
	node->left = node->right = node;
	if (minheap==NULL){
		minheap = node;
	}else{
		// concatenate the root list containing x with root list H
		node->left = minheap->left;
		node->right = minheap;
		minheap->left->right = node;
		minheap->left = node;

		// update the minheap
		if (x<minheap->key){
			minheap = node;
		}
	}
	m_nNode++;

	return node;
}

// no consolidatation here
void FibonacciHeap::insert(FiboNode *node)
{
	node->left = node->right = node;
	node->parent = node->child = NULL;
	node->degree = 0;
	node->mark = false;

	if (minheap==NULL){
		minheap = node;
	}else{
		// concatenate the root list containing x with root list H
		node->left = minheap->left;
		node->right = minheap;
		minheap->left->right = node;
		minheap->left = node;

		// update the minheap
		if ( node->key < minheap->key ){
			minheap = node;
		}
	}
	m_nNode++;
}

KEYTYPE FibonacciHeap::getMin()
{
	if (minheap==NULL){
		cout<<"There is no element in the heap!"<<endl;
		exit(1);
	}
	return minheap->key;
}

void FibonacciHeap::combine(FibonacciHeap* rheap)
{
	// concatenate the root list of H2 with the root list of H
	FiboNode *right0 = minheap->right;
	FiboNode *right1 = rheap->minheap->right;

	minheap->right = rheap->minheap;
	rheap->minheap->right = minheap;

	right0->left = right1;
	right1->left = right0;

	if (minheap==NULL || (rheap->minheap!=NULL && rheap->minheap->key < minheap->key ))
	{
		minheap = rheap->minheap;
	}
	m_nNode += rheap->m_nNode;
}

FiboNode* FibonacciHeap::extractMin()
{
	FiboNode *z = minheap,*ptr,*child;

	if (minheap!=NULL)
	{
		child = minheap->child;

		// for each child x of z add x to the root list of H
		if (child!=NULL)
		{
			ptr = child;
			// set all the child list parent pointer to be NULL
			do {
				ptr->parent = NULL;
				ptr = ptr->right;
			} while(ptr!=child);

			// link it into the root list
			ptr = child->right;
			minheap->right->left = child;
			child->right->left = minheap;
			child->right = minheap->right;
			minheap->right = ptr;
		}

		// remove z from the root list of H (done in the next step)


		if (z == z->right){	// only node in the root list
			minheap = NULL;
		}else{
			minheap = z->right;

			// unlink z from the root list
			z->left->right = minheap;
			minheap->left = z->left;

			consolidate();
		}
		m_nNode--;

	//	KEYTYPE tk = z->key;
	//	delete z;

		return z;
	}else{
		cout<<"Cannot extract element from empty heap!"<<endl;
		return NULL;
	//	exit(1);
	}
}

void FibonacciHeap::consolidate()
{
	size_t i,Dn = 1+(size_t)(log((double)m_nNode)/log(0.5*(1.0+sqrt(5.0))));
	FiboNode **A = new FiboNode*[Dn];
	FiboNode *x,*y;
	size_t d;
	FiboNode *pos = minheap;
	FiboNode *end = minheap->left;
	bool tflag,endflag = false;

	for (i=0;i<Dn;i++)A[i] = NULL;

	while(endflag == false)
	{
		if (pos == end){
			endflag = true;
		}

		tflag = false;
		x = pos;
		d = pos->degree;

		while (A[d]!=NULL)
		{
			y = A[d];
			if (x->key > y->key)
				swap(x,y);

			if (y==pos)
			{	// we will convert a root node to a child
				pos = pos->right;
				tflag = true;
			}

			fibHeapLink(y,x);
			A[d] = NULL;
			d++;
		}
		A[d] = x;
		if (tflag ==false){
			pos = pos->right;
		}
	}

	minheap = NULL;

	for (i=0;i<Dn;i++)
	{
		if (A[i]!=NULL)
		{
			// add A[i] to the root list of H
			if (minheap==NULL){
				minheap = A[i];
				minheap->right = minheap->left = minheap;
			}else{
				A[i]->right = minheap->right;
				A[i]->left = minheap;
				minheap->right->left = A[i];
				minheap->right = A[i];
			}

			// update minheap
			if ( minheap==NULL || A[i]->key < minheap->key ){
				minheap = A[i];
			}
		}
	}

	delete[] A;
}

void FibonacciHeap::fibHeapLink(FiboNode* y,FiboNode* x)
{
#ifdef EXDEBUG
	if (y->parent !=NULL || y->right == y)
	{	// y is not a root node or y is the only root
		cout<<"Illegal operation in fibHeapLink."<<endl;
		exit(1);
	}
#endif

	// remove y from the root list of H
	y->left->right = y->right;
	y->right->left = y->left;

	// make y a child of x
	y->parent = x;
	x->degree++;

	if (x->child==NULL)
	{	// x has no child before
		x->child = y;
		y->right = y->left = y;
	}else{
		y->right = x->child->right;
		y->left = x->child;
		x->child->right->left = y;
		x->child->right = y;
	}

	// unmark y
	y->mark = false;
}

void FibonacciHeap::decrease(FiboNode *x,KEYTYPE k)
{
#ifdef EXDEBUG
	if (x==NULL){
		cout<<"Cannot decrease Null node value!"<<endl;
		exit(1);
	}
#endif
	if ( k >= x->key ){
		cout<<"There is no need to increase the key!"<<endl;
		return;
	}

	x->key = k;
	FiboNode *y = x->parent;

	if( y!=NULL && x->key < y->key )
	{
		cut(x,y);
		cascadingCut(y);
	}

	if (x->key < minheap->key){
		minheap = x;
	}
}

void FibonacciHeap::cut(FiboNode *x,FiboNode *y)
{
#ifdef EXDEBUG
	if (x->parent!=y)
	{
		cout<<"Error link in cut"<<endl;
		exit(1);
	}
#endif

	// remove x from the child list of y,
	if (y->child==x)
	{// need to update the child pointer
		if (x->right ==x)
		{// only one child in the children list
			y->child = NULL;
		}else{
			y->child = x->right;
			x->right->left = x->left;
			x->left->right = x->right;
		}
	}else{	// just unlink x
		x->right->left = x->left;
		x->left->right = x->right;
	}

	y->degree--;

	// add x to the root list of H
	minheap->left->right = x;
	x->left = minheap->left;
	x->right = minheap;
	minheap->left = x;

	x->parent = NULL;
	x->mark = false;
}

void FibonacciHeap::cascadingCut(FiboNode *y)
{
	FiboNode *z = y->parent;
	if (z!=NULL)
	{
		if(y->mark == false)
		{// lost a child
			y->mark = true;
		}else{
			cut(y,z);
			cascadingCut(z);
		}
	}
}
/*
void FibonacciHeap::remove(FiboNode *x)
{
	KEYTYPE k = getMin();
	decrease(x,k-1);	// speciall Note. If Key Type is not integer. This step would fail
	extractMin();
}
*/

void FibonacciHeap::dispose(FiboNode *y)
{
	if (y==NULL)return;
	FiboNode *pos = y->child;
	delete y;
	if (pos==NULL)return;
	FiboNode *end = pos->left;
	while (pos!=end)
	{
		y = pos;
		pos = pos->right;
		delete y;
	}
	delete end;
}
