/**************************************************************************

    begin                : June 30 2004
	version				 : 2.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au

 *   This file contains implementation of  support_vector, SVSetNode and   *
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

#include "interpol.h"   


extern Interpolant*		Parent;   // points to the parent Interpolant class
extern MemoryBlock<SVSetNode>	MBSV;  // memory pool for SVSetNodes

// some static global variables 
	int GlobPos;
	SVSetNodePtr newnode;
	real olddiag;
	SVSetNode *Newnode,*tempnode;
	int AtLeastOneFound;

// global (not in recursion)
	int Globvecnumber_;

	real Meps=1.0e-5;


#define IS_VALID_PTR(a) ((a)!=NULL && (a) !=0xFFFFFFFF)

void CopyNode(SVSetNode* source , SVSetNode* dest) {	memcpy(dest, source, sizeof(SVSetNode)); }

clock_t clockS,clockF;
double TotalTime;

void   ResetTime() {TotalTime=0; clockS=clock();}
double ElapsedTime()
{
	clockF=clock();
	double duration=(double)(clockF - clockS) / CLOCKS_PER_SEC;
	TotalTime += duration;
	clockS=clockF;
	return TotalTime;
}


/**************  support_vector  ********************************************/

void support_vector::SVForm(dVec& x, real val)
// forms the support vector from x and value
{
	vec.newsize(x.size());
	funvalue=val;
	for(int i=0;i<vec.size();i++) 
		vec[i]=funvalue/Parent->m_Constants[i] - x[i];
}

void support_vector::SVForm(real* x, real val)
// forms the support vector from x and value
{
	vec.newsize(Parent->Dim);
	funvalue=val;
	for(int i=0;i<vec.size();i++) 
		vec[i]=funvalue/Parent->m_Constants[i] - x[i];
}

void	support_vector::SVForm(int Label)
{
	label=Label;
	for(int i=0;i<vec.size();i++) 
		vec[i]=funvalue/Parent->m_Constants[i] - vec[i];
}

void support_vector::Increment()
{	ChangeF(funvalue + max_(funvalue,1)*Meps);  }

short int	support_vector::ChangeF(real newval)
{
	short int i;
	for( i=0;i<vec.size();i++) 
		vec[i]+= (-funvalue + newval)/Parent->m_Constants[i];
	i= (newval>=funvalue);
	funvalue=newval;
	return i;
}

void support_vector::ReturnX(dVec& x)
{
	for(int i=0;i<vec.size();i++)
		x[i]=funvalue/Parent->m_Constants[i] - vec[i];
}


/**************  svsetnode  ********************************************/



SVSetNode::SVSetNode()   { 	Init(); }

void SVSetNode::Init() { 
		children__=MB_BADINDEX;
		vecnumber=MB_BADINDEX; 
		Dval=(float) -Infinity;	// for triangulation test
}

SVSetNode::~SVSetNode() { }

void SVSetNode::CopyTo(SVSetNode* copy) { 
	copy->Dval=Dval;
	copy->vecnumber= vecnumber;
}


void SVSetNode::Clear()
{
	int i;
	if(!IsValid()) return;
	if(children__ != MB_BADINDEX)
	for(i=0;i< GetNumChildren();i++) 
		MBSV.FreeBlock(children__ + i);  // may be process all children for killers? ??????
	children__=MB_BADINDEX; // shall we care?
}


void SVSetNode::RemoveChild(SVSetNodePtr child, int pos)	
{		
		SVSetNode* n=MBSV.GetAt(children__ + pos);
		n->Clear();
}

inline int SVSetNode::IsValid()
{
//	if(POS_VEC(vecnumber) != 0xFF) return 1; else return 0;
	return 1;
}

inline SVSetNodePtr SVSetNode::GetParent() { 	return 0; }


int SVSetNode::TestVectorIndex(dVec& v, siVec* index)
{
	int i,flag=0;
	real u;				// make them global?
	if(index==NULL) 
	for(i=0; i<v.size(); i++) {
		u=Parent->SVectors[ i ].vec[i] - v[i];
		if( u > 0) {return 0;}
		if(u==0) 
		{ flag +=1;
		} 
	}
	else
	for(i=0; i<v.size(); i++) {
		u=Parent->SVectors[ (*index)[i] ].vec[i] - v[i];
		if( u > 0) {return 0;}
		if(u==0) { 
			flag+=1; // indicates nonstrict dominance
		}	
	}
	if(flag==1) return 2;  // nonstrict dominance, test not failed but an extra minimum (copy) is needed
	else 
		return 1;  // strict dominance, test failed
}

int SVSetNode::TestVector(dVec& v, siVec* index)
{
		real u;
		int pos = POS_VEC(vecnumber);
		int vecnumber_ = VEC_VEC(vecnumber);

		if(IS_ROOT(vecnumber)) { // root, use index explicitly
			return TestVectorIndex(v,index);
		} else { // test only one diagonal element, assumes parent has passed the test
			u = Parent->SVectors[ vecnumber_ ].vec[pos] -  v[pos];
			if( u >  0) { return 0;}  
			if( u == 0) 
			{	return 2; }      // nonstrict dominance
			return 1; // strict dominance
		}
	//	return 1; //should not get here
}

// as above but strict inequality
int SVSetNode::TestVectorIndexQ(dVec& v, siVec* index)
{	int i;
	for(i=0; i<v.size(); i++) {
		if( Parent->SVectors[ (*index)[i] ].vec[i] > v[i]) {return 0;}
	}
	return 1;
}

int SVSetNode::TestVectorQ(dVec& v, siVec* index)
{
		int pos = POS_VEC(vecnumber);
		int vecnumber_ = VEC_VEC(vecnumber);

		if(IS_ROOT(vecnumber)) { // root, use index explicitly
			return TestVectorIndex(v,index);
		} else { // test only one diagonal element, assumes parent has passed the test
			if( Parent->SVectors[ vecnumber_ ].vec[pos] > v[pos]) {return 0;} //<=
			return 1;
		}
//		return 1; //should not get here
}

//returns 0 if insuccessful, 1 if cond 1 is satisfied
int  SVSetNode::TryNewVectorIndex(support_vector* SV, int pos, siVec* index, real &olddiag )
{
	// test of cond (1). Cond (2) supposedly passed. Only one column is tested. 
	// olddiag is returned for another method (to calculate Dval)
	olddiag= Parent->SVectors[ (*index)[pos]].vec[pos]; // old diagonal
	real val=SV->vec[pos];
	for(int i=0;i<Parent->Dim;i++) { // if diag not smaller than this column, exit 0
		if((pos!=i) && (val >  Parent->SVectors[ (*index)[i]].vec[pos]) ) return 0;   
	}
	return 1;
}

void SVSetNode::GenerateInitVector(siVec* initvec, SVSetNodePtr thisnode)
{
	// this can be called only for one of the roots, to generate its
	// index set, for subsequent processing of children

	if(initvec==NULL) return;

// will retrieve init vector from the killers list
	if(IS_ROOT(vecnumber) ) // means root
	{		// attempt to find this root in the list of heads
	deque<HeadStruc>::iterator iter;

	for(iter=Parent->HeapPossibleMin.Heads.begin(); iter!=Parent->HeapPossibleMin.Heads.end(); iter++)
		if((*iter).Head == thisnode) {
			(*initvec) = *((*iter).p_index);
			return;
		}
		 for(int i=0;i<Parent->Dim; i++) (*initvec)[i]=i;
	}
}



// it is possible to keep this value instead of calculating it every time

real	SVSetNode::ComputeMaximumF(siVec* index)
{
	if(Dval>-Infinity) return Dval;

	real a=-Infinity;
	real b=Infinity;
	real c;
	for(int i=0;i<Parent->Dim;i++) {
		c=Parent->SVectors[ (*index)[i]].funvalue;
		if(a<c) a=c; // max
		if(b>c) b=c; // min
	}
		Dval= (float) (a + 0.001   * (a - b  + 0.0000001)); // max - min
		return (double) Dval;
}

real	SVSetNode::ComputeFunValue(dVec& X, siVec* index)
{
	real b,c=0,a=-Infinity;
	for(int i=0;i<Parent->Dim;i++) {
		b = Infinity;
		for(int j=0;j<Parent->Dim;j++){
			c= Parent->m_Constants[j] * (Parent->SVectors[ (*index)[i]].vec[j] + X[j]);
			if(c < b)  b = c;
			if(b<a) break; // too small
		}
		if(a<c) a=b;
	}
	return a;
}




///////////////////////////////Forest////////////////////////

void Forest::Init() 
{
	size=sizemem=0; sizevirtual=0;

	m_initvec.newsize(Parent->Dim); 
	temp_index.newsize(Parent->Dim); 
	m_index.newsize(Parent->Dim);
	SVT.vec.newsize(Parent->Dim);

	m_TempChildren=MBSV.GetNextFree(Parent->Dim);
	Initiated=1;
};

/*
int Forest::ComputeSize(SVSetNodePtr node)
{
	int sz=1;
	int i;
	SVSetNode* n=MBSV.GetAt(node);
	if(n->IsValid())
	if( n->children__!=MB_BADINDEX) {
		for(i=0; i<n->GetNumChildren(); i++) {
			sz += ComputeSize(n->children__ + i);
		}
	}
	return sz;
}

int Forest::ComputeSize(SVSetNodePtr node, int not_this_child)
{
	int sz=1;
	int i;
	SVSetNode* n=MBSV.GetAt(node);
	if(n->IsValid())
	if( n->children__!=MB_BADINDEX) 
		for(i=0;i<n->GetNumChildren();i++) {
			if(i!=not_this_child) sz += ComputeSize(n->children__ + i);
	}
	return sz;
}

*/


void	Forest::AddRootNode( SVSetNodePtr node)
{
// this method called only once in the initpopulate
	SVSetNode* n=MBSV.GetAt(node);

	SET_ROOT(n->vecnumber);
	n->SetNumChildren(-1);
	HeadStruc Head;
	Head.Head=node;
	Head.p_index = new(siVec);
	*(Head.p_index) = temp_index;

	Heads.push_back(Head);
	size++;
}

void	Forest::AddLeaf(SVSetNodePtr node)
{
	// inserts the leaf into the heap, otherwise calls recursively
	int i;
	SVSetNode* n= MBSV.GetAt(node);

	if(n->GetNumChildren()<=0 /*==MB_BADINDEX*/) // means this is a leaf
	{
		size++; // size of tree
		return;
	}
	else { // process children if any
		for(i=0;i<n->GetNumChildren();i++) {
			// the child can be empty
			 AddLeaf(n->children__ +i);   // where is memory allocation?
		}
		size++; 
	}
}

void	Forest::AddTree(SVSetNodePtr root)
{
	if(root==MB_BADINDEX) return;

	HeadStruc Head;
	Head.Head=root;
	Head.p_index = new(siVec);
	*(Head.p_index) = temp_index;
	Heads.push_back(Head);

//	SVSetNode* R= MBSV.GetAt(root);

	AddLeaf(root);
}



void	Forest::EraseBranch(SVSetNodePtr branch, int processparent)
{
	SVSetNode* Branch= MBSV.GetAt(branch);

	if(!Branch->IsValid()) return;

	int i;
	if(Branch->GetNumChildren()>0/*!=MB_BADINDEX*/) {
		for(i=0;i<Branch->GetNumChildren();i++)
			EraseBranch(Branch->children__ +i);
	}
	
	MBSV.FreeBlock(branch);	
}

void	Forest::EraseRootEntry(SVSetNodePtr branch)
{
	// find the root in the list of roots and erase it.
	deque<HeadStruc>::iterator iter;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++)
		if((*iter).Head == branch) {
			delete (*iter).p_index;
			iter=Heads.erase(iter); 
			return;
		}
}

void	Forest::ClearBranch(SVSetNodePtr branch)
{
// this method is called from destructor. It differes from EraseBranch in that
// the nodes are not deleted from the heap (to save time)
	SVSetNode* Branch= MBSV.GetAt(branch);
	if(!Branch->IsValid()) return;

	int i;
	if(Branch->GetNumChildren()>0/*!=MB_BADINDEX*/) {
		for(i=0;i<Branch->GetNumChildren();i++)
			ClearBranch(Branch->children__ +i);
	}

// destructor, keeps in the heap invalid reference
	size--;
	MBSV.FreeBlock(branch);	
}




void Forest::EraseAll()
{
	deque<HeadStruc>::iterator iter;

	if(MBSV.IsValid() )
		for(iter=Heads.begin(); iter!=Heads.end(); iter++) { delete (*iter).p_index; ClearBranch((*iter).Head);  }

//	MBCL.ClearAll();
//	MBSV.ClearAll();
/**/
	Heads.clear();
	size=0;
	Initiated=0;
}





int Forest::ProcessNode(SVSetNodePtr node, support_vector* SV, siVec* initvec)
{
	// recursive calls
	SVSetNode* Node=MBSV.GetAt(node);

	int P,i,numchld;


/*------------------- for indices --*/
#ifdef TRIANGULATION1
	UINT	idx=0;
	int Globpos;
	Globpos = POS_VEC(Node->vecnumber);

	int Globvecnumber_ = VEC_VEC(Node->vecnumber);
	if(!IS_ROOT(Node->vecnumber) ) {
		idx = (*initvec)[Globpos];
		(*initvec)[Globpos]=Globvecnumber_; // update the svset
	}
#endif
/*------------------- for indices --*/


	// here I can keep indexvector (starting from the top), so no computevectors is necessary.
	P=Node->TestVector(SV->vec, initvec);


	switch(P) {
	case 2: 
			// remember position to restore after return
//			Pos=GlobPos;
		
	case 1: // dominance, split this node
		if(Node->GetNumChildren()<0 ) { // means leaf
			AtLeastOneFound++;

			if(P==2 ) 
				SV->Increment();  // data not in general position, perturb the data

			// must be no children__ at this stage
			// create children__ if any
			numchld=0;

			for(i=0;i<Parent->Dim;i++)
				if(Node->TryNewVectorIndex(SV, i, initvec, olddiag)) { // cond (1)

					// add new node
					size++;
					newnode = m_TempChildren + numchld;
					numchld++;

					Newnode = MBSV.GetAt(newnode);
					Newnode->Init();
					Newnode->vecnumber=POSVEC_VEC(i, SV->label);
					Newnode->SetNumChildren(-1);

				} 

				Node->SetNumChildren(numchld);
				Node->children__= MBSV.GetNextFree(numchld);
				sizemem+=numchld;

				for(i=0;i<numchld;i++) { // copy to the actual nodes
					newnode=Node->children__ + i;
					tempnode= MBSV.GetAt(m_TempChildren + i);
					Newnode = MBSV.GetAt(newnode);
					CopyNode(tempnode, Newnode);
				}
			
		} else { // this was a branch, recursively process children
			for(i=0;i<Node->GetNumChildren();i++) 
			{
				 ProcessNode(Node->children__ + i, SV, initvec); // everything is done here
			}
		}
		// on exit undo IDX
		break;

	case 0: //test (2) passed
		// do nothing
		break;
	}

/*------------------- for indices --*/
#ifdef TRIANGULATION1
 	if(!IS_ROOT(Node->vecnumber) ) {
		(*initvec)[Globpos]=idx;
	}
#endif

	return 0;
}

void Forest::ProcessAll(support_vector* SV)
{
	deque<HeadStruc>::iterator iter;
	SVSetNode *node;

	AtLeastOneFound=0;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++) {
		node=MBSV.GetAt( (*iter).Head );
		m_initvec = *((*iter).p_index);				// root is not necessarily {1,2,3,...} on multiprocessor system
		ProcessNode( (*iter).Head, SV, &m_initvec);
	}
	if(AtLeastOneFound==0)
	{
		Parent->Match=1;
//		cout << SV->label <<" "<<SV->funvalue<<endl;
	}
}




int Forest::ProcessNodeDyn(SVSetNodePtr node, support_vector* SV, siVec* initvec)
{
	// recursive calls
	SVSetNode* Node=MBSV.GetAt(node);

	short int i;
	UINT	idx=0;
	int Globpos;

	Globpos = POS_VEC(Node->vecnumber);

	int Globvecnumber_ = VEC_VEC(Node->vecnumber);
	if(!IS_ROOT(Node->vecnumber) ) {
		idx = (*initvec)[Globpos];
		(*initvec)[Globpos]=Globvecnumber_; // update the svset
	}

	// need to know index in advance
	real avr	= Node->ComputeMaximumF(initvec);

	// adapt SV to this node 
	i=SV->ChangeF(avr);

// need to do index , as the SV has changed!!
	if(i) i=Node->TestVectorQ(SV->vec, initvec)==1; else i=Node->TestVectorIndexQ(SV->vec, initvec)==1; 
	if(i) 
// equivalent to the next line, but may not work on every compiler
//	if( (i && Node->TestVectorQ(SV->vec, initvec)==1) || (Node->TestVectorIndexQ(SV->vec, initvec)==1) )  // skips full test if avr did not decrease

	{  // dominance, split this node
		if(Node->GetNumChildren()<0) { // means leaf 
			for(i=0;i<Parent->Dim;i++)
				m_indexset.insert((*initvec)[i]);

		} else { // this was a branch, recursively process children
			for(i=0;i<Node->GetNumChildren();i++) 
			  { // save index
				 ProcessNodeDyn(Node->children__ + i , SV, initvec);
				 // restore index
			  }
		}

	}// do nothing otherwise

 	if(!IS_ROOT(Node->vecnumber) ) {
		(*initvec)[Globpos]=idx;
	}

	return 0;
}
	
void Forest::ProcessAllDyn(support_vector* SV)
{
	m_indexset.clear();
	if(!Initiated) return;

	deque<HeadStruc>::iterator iter;
	SVSetNode *node;

	for(iter=Heads.begin(); iter!=Heads.end(); iter++) {
		node=MBSV.GetAt( (*iter).Head );
		m_initvec = *((*iter).p_index);
		ProcessNodeDyn( (*iter).Head, SV, &m_initvec);
	}
}




/*-----------------------------------------------------------
	Packing routines: to transfer branches between processors
	Only to be used on multiprocessor system under MPI
	 not finished yet...
*/

#define CODE_EL 1
#define CODE_EL_ROOT 4
#define CODE_CHILDREN 2
#define CODE_KILLERS 3
#define CODE_INDEX 5
#define CODE_END 6
#define CODE_CONSTRAINED_MIN 7

void	Forest::PackBranch(SVSetNodePtr branch, char* buffer, int& pos)
{
// pack itself
	int i,j;
	sizepacked++;
	SVSetNode* Branch=MBSV.GetAt(branch);


	buffer[pos++]=CODE_EL;

	memcpy(buffer+pos,&(Branch->vecnumber), sizeof(int)); pos +=sizeof(int);
	if(!Branch->IsValid()) { buffer[pos++]=CODE_END; return;}

	memcpy(buffer+pos,&(Branch->Dval), sizeof(Branch->Dval)); pos +=sizeof(Branch->Dval);

	if(pos>0x7FFFFF) {pos=0; return; }
// packs children if any
	if(Branch->GetNumChildren()>0 /*Branch->children__!=MB_BADINDEX*/) {
		buffer[pos++]=CODE_CHILDREN; // recursion
		j=Branch->GetNumChildren();
		memcpy(buffer+pos,&(j),sizeof(int)); pos +=sizeof(int);
		for(i=0;i<Branch->GetNumChildren();i++) 
			{
				PackBranch(Branch->children__ + i,buffer,pos);
			} 
	}

	buffer[pos++]=CODE_END; // terminate this element
}

void	Forest::PackBranchStart(SVSetNodePtr branch, char** buffer, int* pos)
{
// starts packing, calls pack recursively

	*buffer=(char*) calloc(0x7FFFFF,1); // Size???
	*pos=0;
	sizepacked=0;


	if(branch==MB_BADINDEX) {(*buffer)[(*pos)++] = CODE_END; return;}

	// ensure this branch has no possible killers ???
	SVSetNode* Branch=MBSV.GetAt(branch);

	siVec index(Parent->Dim);
//	Branch->ComputeVectors(index,&m_initvec,branch);

	// now index
	(*buffer)[(*pos)++]=CODE_INDEX;
	memcpy(*buffer+*pos, index.begin(), index.size()*sizeof(SVINDEX));
		*pos += index.size()*sizeof(SVINDEX);

	int savevec=Branch->vecnumber;
	SET_ROOT(Branch->vecnumber); // forget the vecnumber
	PackBranch(branch, *buffer, *pos);
	if(*pos==0) { // means too many nodes
		return;
	} 

	(*buffer)[(*pos)++]=CODE_END;

	Branch->vecnumber=savevec;

	EraseBranch(branch,1);
}

void	Forest::UnPackBranch(char* buffer, int& pos, SVSetNodePtr branch)
{
	int i,j,cont=1;
//	unsigned int k;
	SVSetNode* Child;
	SVSetNodePtr child;
	SVSetNode* Branch=MBSV.GetAt(branch);


	while(cont) {
		i=buffer[pos++];
		switch(i) {
		case CODE_END: cont=0; break;

		case CODE_INDEX:
			SET_ROOT(Branch->vecnumber);
			Branch->parent__=MB_BADINDEX;
			memcpy(temp_index.begin(),buffer+pos,temp_index.size()*sizeof(SVINDEX));
			pos += temp_index.size()*sizeof(SVINDEX);
				break;


		case CODE_EL:
			memcpy(&(Branch->vecnumber),buffer+pos, sizeof(int)); pos +=sizeof(int);
			if((Branch->vecnumber & 0xFFFFFF) != 0xFFFFFF)  // valid
			{memcpy(&(Branch->Dval),buffer+pos, sizeof(Branch->Dval)); pos +=sizeof(Branch->Dval);}
			else EraseBranch(branch);

			break;

		case CODE_CHILDREN:
			memcpy(&(j),buffer+pos, sizeof(int)); pos +=sizeof(int);
			Branch->SetNumChildren(j);
			Branch->children__=MBSV.GetNextFree(j);

			for(i=0;i<j ;i++) {
				if(buffer[pos] != 0) { // if next one is OK
					child=Branch->children__+i;
					Child=MBSV.GetAt(child);
					Child->Init();
					Child->parent__=branch;
					Child->SetNumChildren(-1);

					//	new(SVSetNode);
					UnPackBranch(buffer,pos,child);
				}
				else pos++;
			}
				break;

		}
	}
}

void	Forest::UnPackBranchStart(char** buffer, SVSetNodePtr* branch )
{
	SVSetNodePtr parentnode=MBSV.GetNextFree();
	SVSetNode* Parent=MBSV.GetAt(parentnode);  //new(SVSetNode);
	Parent->Init();

	int pos=0;
	UnPackBranch(*buffer,pos, parentnode);

	if(pos==1) // empty branch
	{
		MBSV.FreeBlock(parentnode);
	}

	free(*buffer);
	*buffer = NULL;
	*branch=parentnode; // for return
}

