/**************************************************************************

    begin                : June 30 2004
	version		 : 1.2 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au

 *   This file contains several classes: support_vector, SVSetNode and     *
 *   Forest.                                                               *
 *                                                                         *
 *      support_vector is a vector, label and a value, as used in the      *
 *      cutting angle method.                                              *
 *      SVSetNode    represents a combination of n support vectors         *
 *      when organiser in a tree (ie it's a tree node)                     *
 *      Forest is a set of trees of SVSetNodes                             *
 *
	    Forest takes care of maintaining the tree structures, in which 
		parent nodes have references to children nodes, and allows queries
		starting from the root(s)

        SVSetNode allows to perform certain tests on nodes and does some
		housekeeping

		These classes are not to be used directly but from within
		Interpolant class. These are workers which perform all computations
		required by Interpolant. 

	See documentation about the methods used for further information

 *                                                                         *
 *  © Gleb Beliakov, 2004								   *
 * This program is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful, but     *
 * WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * General Public License for more details.                                *
 *                                                                         *
 * You should have received a copy of the GNU General Public License       *
 * along with this program; if not, write to the Free Software Foundation, *
 * Inc., 59 Temple Place Suite 330, Boston, MA 02111-1307 USA.             *
 ***************************************************************************/


#ifdef _MSC_VER  
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINT; //this type myst be 8 bytes long
#define ULLMASK 0xFFFFFFFF00000000UL 
//#define NOMINMAX 
#else
#define ULLMASK 0xFFFFFFFF00000000ULL 
typedef unsigned long long int ULINT; //this type myst be 8 bytes long
#endif

// this is compiler dependent
// these macros are to cast double to 8 byte integer and back
#define d2ulint(a)  (*(ULINT*) &(a))
#define ulint2d(a) (*(double*) &(a))


#define TNT_NO_BOUNDS_CHECK
#define real double // 8 bytes: choose one of these

typedef unsigned int SVINDEX ;

// to save RAM for functions with < 15 variables. 
#define SMALL_DIM  // up to 15 variables, for more than 15 vars, undefine it

#define TRIANGULATION1

#include <cstdio>
#include <cstdlib>
#include <list>
#include <cmath>
#include <ctime>
#include <deque>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;

#include "tnt/tnt.h"
#include "memblock.h"

//using namespace TNT;
typedef TNT::Vector<float>		fVec;
typedef TNT::Vector<real>		dVec;
typedef TNT::Vector<int>		iVec;
typedef TNT::Matrix<real>		dMat;
typedef TNT::Matrix<int>		iMat;
typedef TNT::Vector<SVINDEX>	siVec; 

typedef vector<char> shortindexvector;
typedef vector<char>::iterator shortindexvectoriter;
typedef set<UINT> indexset;
typedef set<UINT>::iterator indexsetiter;

const real Infinity  = 1.0e16;
const real SInfinity = 1.0e7; // "small infinity" for boundary points (to avoid loss of precision when taking 1-x[i])

#define sqr_(a) ((a)*(a))
#define min_(a,b) ((a)<(b)?(a):(b))
#define max_(a,b) ((a)>(b)?(a):(b))
double ElapsedTime();
void   ResetTime() ;


#ifndef SMALL_DIM

#define POS_VEC(a) (((a)>>24) & 0xFF)   // only the last 6 bits, the first 2 bits reserved

#else
// pos 3F would mean root (1,2,3,...,n), no more than 255 variables
#define POS_NUMCHLF(a) (((a)>>28))   // + 1 as 0 children does not make sense, so 0 means 1, and to ensure F is NOT used (special for root)
#define POS_VEC(a) (((a)>>24) & 0x0F)   // only 15 possible positions

#endif

#define POSVEC_VEC(a,b) ((a)<<24 | (b))
#define VEC_VEC(a) ((a) & 0x00FFFFFF)
 
#define IS_ROOT(a) (((a) &0xFFFFFF) == 0xFFFFFF) 
#define SET_ROOT(a) ((a) |= 0xFFFFFF) 




/* ---------- Aux classes------------------------*/

#define SVSetNodePtr UINT

// just to pack these into 24 bytes or less
#define vecnumber SVSetNodeData[0]
#define children__ SVSetNodeData[1]
#define numchildren__ SVSetNodeData[2]

#ifndef SMALL_DIM
#define parent__ SVSetNodeData[3]
#else
#define parent__ SVSetNodeData[2]
#endif




class support_vector {
public:
	unsigned int label;
	dVec		 vec;
	real		 funvalue;

// create a support vector from x and val
	void	SVForm(dVec& x, real val);
// create a support vector from x and val (different syntaxix)
	void	SVForm(real* x, real val);

// same as above, byt x and val are already stored in vec and funvalue
	void	SVForm(int Label);

// increment funvalue by meps to break the ties.
	void	Increment();

// update the components if the function value changes
	short int	ChangeF(real newval); // returns 1 if old value < newval

// returns the coordinates of the point x
	void ReturnX(dVec& x);

	support_vector* This() {return this;};
};

// to store the list of support vectors
typedef deque <support_vector>	SVDeque;


// One local minimum of  function H -- as a node of the tree

class SVSetNode {
public:

#ifdef SMALL_DIM
	UINT	SVSetNodeData[2]; // packs children and parent  8 bytes only!!!
#else
	UINT	SVSetNodeData[3]; // packs children and parent
#endif
	float Dval;				 // value of the max of vertices, use float to save space 
//  real  Dval;  

// end data members--------------------


	SVSetNode(); // constructor, assigns NULL to pointers
	~SVSetNode(); 
	void		Init();
	SVSetNode*	This() {return this;}
	int			IsValid();
	SVSetNodePtr GetParent() ;

// how many children this node has
#ifdef SMALL_DIM
	int	GetNumChildren() { 
		UINT a=POS_NUMCHLF(vecnumber);
		if((a - 15) <= 0) return -1; 
		else return a; };
	void SetNumChildren(int ncld) { 
		if(ncld==-1) vecnumber|=0xF0000000; else {
			vecnumber &= 0x0FFFFFFF; vecnumber |= (ncld) << 28; }};
#else
	inline int	GetNumChildren() {if(numchildren__ < 0xFFFFFFFF) return numchildren__; else return -1; };
	inline void SetNumChildren(int ncld) {if(ncld<0) numchildren__=0xFFFFFFFF; else numchildren__ = ncld; };
#endif

	// attaches a child "node" to this, at position pos
	void AddChild(SVSetNodePtr thisnode, SVSetNodePtr node, int pos);

	// deletes all children. Used to clear memory when destroying the tree
	void Clear();

	// removes just the reference to the child, not destroys the child
	void RemoveChild(SVSetNodePtr child, int pos);//	{	children[pos]=NULL; }

	// these two methods test cond (2) for SVector v
	// the first version is to test nodes other than root (index is not important)
	// the second version is to test root, in which case index should be the
	// list of indices of SV comprising this node
	// returns 0 if passed, 1 if failed (dominance), 2 if nonstrict dominance, and 3 if below best function value,
	int TestVector(dVec& v, siVec* index);
	int TestVectorIndex(dVec& v, siVec* index);
	int TestVectorIndexQ(dVec& v, siVec* index);
	int TestVectorQ(dVec& v, siVec* index);


	// for the root node generates the indices of SVectors. for ROOT returns 1,2,3,,,.n
	// otherwise returns the acural indices, stored in VectorPos
	void GenerateInitVector(siVec* initvec, SVSetNodePtr thisnode);

	// tests cond (1) with SV at position pos. Assumes that the parent
	// satisfies this condition, and hence tests only column pos
	// index contains the actual SV indices. Also returns the olddiag, the value
	// of the element on diagonal to be replaced. It will be used in updating DVal
	int TryNewVectorIndex(support_vector* SV, int pos, siVec* index, real &olddiag);

	void CopyTo(SVSetNode* copy);

// computes the maximum of funvalues of the participating support vectors
	real	ComputeMaximumF(siVec* index);

// computes the value of the local minimum
	real	ComputeFunValue(dVec& X, siVec* index);

};



/*-----------------------------------------------------------------------------
	This class implements a tree (rather forest). Leaves are the local minima of saw-tooth
	cover. 


	There are 2 types of methods, the routine insert/delete and problem-specific queries
	see documentation about the methods used


  The root keeps its participating support vectors in full

------------------------------------------------------------------------*/
struct HeadStruc {
	SVSetNodePtr Head;
	siVec*	p_index;
};

class Forest {
public:
	int size, sizevirtual, sizemem; // aux. members for testing
	int sizepacked;

	siVec m_initvec,	m_index, temp_index; // just not to create it in all functions
	// provide temp. storage passed to through pointer

	support_vector		SVT;

	deque<HeadStruc>	Heads; // here we keep the roots of the trees

	SVSetNodePtr		m_TempChildren;

	int					Initiated;  // flag to indicate the forest has a root

	indexset			m_indexset;

public:
// constructor
	Forest() {Initiated=0; m_indexset.clear();};

	void Init(); // to create Heap and aux. storage
// destructor
	void	EraseAll();

// Routine methods
	int GetSizeMem()	{return sizemem; }; // returns the size of the forest
	int GetSize()		{return size; };	// returns the size of the forest
	size_t SizeRoot()		{return Heads.size();}; // how many roots
	size_t Size()			{return (SizeRoot() <<24); };  //not used

	int ComputeSize(SVSetNodePtr node);  // size of this branch
	int ComputeSize(SVSetNodePtr node, int not_this_child);  // same but excluding this child's branch

										
	void	AddRootNode(  SVSetNodePtr node);	// starter: called in InitPopulate
	void	AddTree(SVSetNodePtr root);			// add a branch
	siVec*	GetVecAddress() {return &m_initvec;}; // provides working memory

	void	AddLeaf(SVSetNodePtr node); // called recursively to find the leafs and insert into the heap
private:	
	void	ClearBranch(SVSetNodePtr branch); //like EraseBranch, but not removed from heap

public:
	void	EraseBranch(SVSetNodePtr branch, int processparent=-1); 
	void	EraseRootEntry(SVSetNodePtr branch); // like erase branch, but processes roots

// these are problem-specifis methods
private:
// called internally from ProcessAll. This is the working horse
// given new SV, and the root index vector *initvec (calculated before the first call)
// returns  1 if test (2) fails (needs to split this node). If there are children,
//	processes them recursively
// returns 0 if not affected by SV. In this case processing stops (children not processed)
	int		ProcessNode(SVSetNodePtr node, support_vector* SV, siVec* initvec);


public:
// called outside. Starts at roots and processes all trees in the forest. Splits
// and updates the tree automatically
	void	ProcessAll(support_vector* SV);

// as above, by the SV changes dynamically 
	void	ProcessAllDyn(support_vector* SV);
	int		ProcessNodeDyn(SVSetNodePtr node, support_vector* SV, siVec* initvec);


// these methods are to transfer branches between processors
// PackBranchStart and UnPackBranchStart should be called for specified branch
//
	void	PackBranchStart(SVSetNodePtr branch, char** buffer, int* pos);
	void	UnPackBranchStart(char** buffer, SVSetNodePtr* branch);
	void	PackBranch(SVSetNodePtr branch, char* buffer, int& pos);
	void	UnPackBranch(char* buffer, int& pos, SVSetNodePtr branch);


};


