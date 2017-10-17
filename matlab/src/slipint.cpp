/**************************************************************************

    begin                : Apr 19 2005
	version				 : 2.0 
    copyright            : (C) 2005 by Gleb Beliakov
    email                : gleb@deakin.edu.au


  SLipInt.cpp: declaration of the Simple Lipschitz interpolant class.

  SLipInt class implements the method of Lipschitz interpolation
  and smoothing. The interpolant is computed as
  g(x)= 0.5(H_upper(x)  + H_lower(x))

  with 
  H_upper(x)= min_k (y^k + LipConst d(x,x^k))
  H_lower(x)= max_k (y^k - LipConst d(x,x^k))

  where the input data is (x^k,y^k), k=1,...npts.

  This is the best interpolant in the worst case scenario, if the interpolated
  function is known to be Lipschitz with the Lipschitz constant LipConst.

  There are no restrictions on the distribution of data x^k in R^dim

  The enhancements in version 2 include smoothing, monotone approximation,
  automatic calculation of the Lipschitz constant using sample splitting and
  cross-validation.

  See documentation for more details.

 *                                                                         *
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

#include "slipint.h"




#ifdef _MSC_VER  
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINT; //this type myst be 8 bytes long
#else
typedef unsigned long long int ULINT; //this type myst be 8 bytes long
#endif


double sqr__(double a) {return a*a; }
double max__(double a, double b) { return((a>b)?a:b); }
double min__(double a, double b) { return((a<b)?a:b); }

// various distances

// Euclidean
double SLipInt::dist(int dim,  double* x, double* xk, double* param)
{
	double d=0;
	dim--;
	for( ; dim>=0; dim--) d+=sqr__(x[dim] - xk[dim]);
	return sqrt(d);
}


// l_infty
double SLipIntInf::dist(int dim, double* x, double* xk, double* param)
{
	double dk,d=-1;
	int i = dim-1;
	for( ; i>=0; i--) {
		dk=fabs(x[i] - xk[i]);
		if(dk>d) {
			d=dk;
		}
	}
	return (d);
}
// l_infty, also returns the direction
double SLipIntInf::distInfDir(int dim,  double* x, double* xk, int* dir)
{
	double dk,d=-1;
	int i=dim-1;
	for( ; i>=0; i--) {
		dk=(x[i] - xk[i]);
		if(fabs(dk)>d) {
			d=fabs(dk);
			if(dk>=0) *dir=i; else *dir=dim+i;
		}
	}
	return (d);
}

// l_infty, returns the direction
double SLipIntInf::distDir(int dim,  double* x, double* xk, int* dir)
{
	double dk,d=-1;
	int i = dim-1;
	for( ; i>=0; i--) {
		dk=fabs(x[i] - xk[i]);
		if(dk>d) {
			d=dk;
			*dir=i; 
		}
	}
	return (d);
}




// Constrained distances

double SLipInt::dist(int dim,  double* x, double* xk, int* Cons, double* param)
{ //||(x-xk)V+||
	double d=0;
	dim--;
//	for( ; dim>=0; dim--) d+=sqr__(  max__(  Cons[dim] * (x[dim] - xk[dim]), 0) );
	for( ; dim>=0; dim--) d+=sqr__(  Cons[dim]!=0 ? max__(  Cons[dim] * (x[dim] - xk[dim]), 0) : x[dim] - xk[dim] );
	return sqrt(d);
}

double SLipIntInf::dist(int dim,  double* x, double* xk,int* Cons, double* param)
{
	double dk,d=-1;
	int i = dim-1;
	for( ; i>=0; i--) {
		dk =  Cons[i]!=0 ? max__(  Cons[i] * (x[i] - xk[i]), 0) : fabs(x[i] - xk[i]) ;
		//dk=fabs(x[i] - xk[i]);
		if(dk>d) {
			d=dk;
		//	*dir=i; 
		}
	}
	return (d);
}



double SLipInt::distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param)
{ //||(x-xk)V+,alpha||
	double d=0;
	dim--;
    for( ; dim>=0; dim--) {
		if( Cons[dim]==0) d+=sqr__(x[dim] - xk[dim]);
		else if(Cons[dim]>0) d+= sqr__ (max__(  (x[dim] - xk[dim]),  min__(0,  (LeftRegion[dim]-xk[dim])))); 
		else d+=sqr__(max__(  (xk[dim] - x[dim]),  min__(0,  1 *(LeftRegion[dim]-x[dim])) ) );
	} 
	return sqrt(d);
}
double SLipInt::distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param)
{ //||(x-xk)V+,beta||
	double d=0;
	dim--;
	for( ; dim>=0; dim--) 
		if( Cons[dim]==0) d+=sqr__(x[dim] - xk[dim]);
		else if(Cons[dim]>0) d+= sqr__ (max__((x[dim] - xk[dim]),  min__(0, (x[dim]-RightRegion[dim]) ) ) );
		else d+= sqr__ (max__(xk[dim] - x[dim],  min__(0,xk[dim]-RightRegion[dim] )));
	return sqrt(d);
}



double SLipInt::distAll(int dim,  int type, double* x, double* xk, int* Cons, double* Region, double* param)
{
	switch(type) {
		case 0: return dist(dim,x,xk,param);
		case 1: return dist(dim,x,xk,Cons,param);
		case 2: return distLeftRegion(dim,x,xk,Cons,Region,param);
		case 3: return distRightRegion(dim,x,xk,Cons,Region,param);
	}
	return 0;
}



double SLipIntInf::distLeftRegion(int dim, double* x, double* xk, int* Cons, double* LeftRegion, double* param)
{ //||(x-xk)V+,alpha||
	double d=0;
	dim--;
    for( ; dim>=0; dim--) {
		if( Cons[dim]==0) d =max__(d,fabs(x[dim] - xk[dim]));
		else if(Cons[dim]>0) d+= max__(d, max__(  (x[dim] - xk[dim]),  min__(0,  (LeftRegion[dim]-xk[dim])))); 
		else d+=max__(d,max__(  (xk[dim] - x[dim]),  min__(0,  1 *(LeftRegion[dim]-x[dim])) ) );
	} 

	return d;
}

double SLipIntInf::distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param)
{ //||(x-xk)V+,beta||
	double d=0;
	dim--;
	for( ; dim>=0; dim--) 
		if( Cons[dim]==0) d+=max__(d,fabs(x[dim] - xk[dim]));
		else if(Cons[dim]>0) d+= max__ (d,max__((x[dim] - xk[dim]),  min__(0, (x[dim]-RightRegion[dim]) ) ) );
		else d+= max__ (d,max__(xk[dim] - x[dim],  min__(0,xk[dim]-RightRegion[dim] )));
	return sqrt(d);
}

double SLipIntInf::distAll(int dim,  int type, double* x, double* xk, int* Cons, double* Region, double* param)
{
	switch(type) {
		case 0: return dist(dim,x,xk,param);
		case 1: return dist(dim,x,xk,Cons,param);
		case 2: return distLeftRegion(dim,x,xk,Cons,Region,param);
		case 3: return distRightRegion(dim,x,xk,Cons,Region,param);
	}
	return 0;
}



// simplicial distance
double SLipIntInf::distSimp(int dim, double* x, double* xk, int* dir)
// int I, int J)
{	
	double d,t;
	double u1,u2;
	double uu,vv;
	// simplicial distance
	int i;
	d=0;
	u1=u2=0;
	for(i=0;i<dim;i++) {     
		uu=x[i];
		vv=xk[i];
		t=uu-vv;
		u1 += uu;
		u2 += vv;
		if(d<t) {d=t; *dir=i;}      // max(x-y)
	}
	uu=1-u1; // slack variables
	vv=1-u2;
	t=uu-vv;
	if(d<t) { d=t; *dir=dim;}      // max(x-y)

// direction is the coordinate where max is achieved
	return d;
}

double SLipInt::distInf1(int dim, double* x, double* xk, int* dir)
{
	double dk,d=-1;
	int i = dim-1;
	for( ; i>=0; i--) {
		dk=fabs(x[i] - xk[i]);
		if(dk>d) {
			d=dk;
			*dir=i; 
		}
	}
	return (d);
}
// just procedural versions of the methods above
double distInf2(int dim, double* x, double* xk, int* dir)
{
	double dk,dka,d=-1;
	int i = dim-1;
	for( ; i>=0; i--) {
		dk = (x[i] - xk[i]);
		dka=fabs(dk);
		if(dka>d) {
			d=dka;
			if(dk>=0)
				*dir=i;
			else *dir=i+dim;
		}
	}
	return (d);
}



 double SLipIntLp::dist(int dim,  double* x, double* xk, double* param)
{
	double d=0;
	dim--;
	for( ; dim>=0; dim--) d+=pow(x[dim] - xk[dim], m_P);
	return pow(d, m_P1);
}
// Euclidean
 double SLipIntLp::dist(int dim,  double* x, double* xk, int* Cons, double* param) // constrained
 { //||(x-xk)V+||
	double d=0;
	dim--;
	for( ; dim>=0; dim--) d+=pow(  Cons[dim]!=0 ? max__(  Cons[dim] * (x[dim] - xk[dim]), 0) : x[dim] - xk[dim] ,     m_P);
	return pow(d,m_P1);
}

 double SLipIntLp::distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param)
{ //||(x-xk)V+,alpha||
	double d=0;
	dim--;
    for( ; dim>=0; dim--) {
		if( Cons[dim]==0) d+=pow(x[dim] - xk[dim], m_P);
		else if(Cons[dim]>0) d+= pow (max__(  (x[dim] - xk[dim]),  min__(0,  (LeftRegion[dim]-xk[dim]))), m_P); 
		else d+=pow(max__(  (xk[dim] - x[dim]),  min__(0,  1 *(LeftRegion[dim]-x[dim])) ), m_P);
	} 
	return pow(d, m_P1);
}

 double SLipIntLp::distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param)
{ //||(x-xk)V+,beta||
	double d=0;
	dim--;
	for( ; dim>=0; dim--) 
		if( Cons[dim]==0) d+=pow(x[dim] - xk[dim], m_P1);
		else if(Cons[dim]>0) d+= pow (max__((x[dim] - xk[dim]),  min__(0, (x[dim]-RightRegion[dim]) ) ) , m_P1);
		else d+= pow (max__(xk[dim] - x[dim],  min__(0,xk[dim]-RightRegion[dim] )), m_P1);
	return pow(d, m_P1);
}

 double SLipIntLp::distAll(int dim,  int type, double* x, double* xk, int* Cons, double* Region, double* param)
{
	switch(type) {
		case 0: return dist(dim,x,xk,param);
		case 1: return dist(dim,x,xk,Cons,param);
		case 2: return distLeftRegion(dim,x,xk,Cons,Region,param);
		case 3: return distRightRegion(dim,x,xk,Cons,Region,param);
	}
	return 0;
}








/* auxiliary class CLargeSet
This class represents a set of integers up to K (excluding K). K can be large 
(of order of 10^6. The presence or absence of an element is indicated by the
corresponding bit in a mask, represented by an array of 64 bit integers.
Set operations are implemented though bit masks. Should be very fast.
*/
class CLargeSet {
public:
	int	m_size, m_sizearray;
	ULINT* m_els;

	CLargeSet(int K) {m_size=K; m_sizearray=m_size/sizeof(ULINT)/8 + 1;  // els needed
		m_els=(ULINT*) malloc(m_sizearray*sizeof(ULINT));
		Clear();
	}; // how many bytes we need

	~CLargeSet(){free(m_els);};

// sets the bit for element i (adds i to the set)
	void	Set(int i);
// sets the bit for element i (adds i to the set), and returns 1 if i was already there 0 otherwise
	int		SetChanged(int i);

// removes i from the set
	void	Remove(int i);
// clears the set = empty set
	void	Clear();
// does i belong to the set?
	int		IsPresent(int i);
// what is the next element in this set
	int		NextElement(int i);
// removes from this set elements that belong to OtherSet
	void	Remove(CLargeSet* OtherSet);
// copies elements from the OtherSet to this one (OtherSet must be at least as big as this one)
	void	Copy(CLargeSet* OtherSet);

private:
	int		pos;
	ULINT	mask;
// compute the position of the bit for element i in the mask
	void	ComputeBit(int i);
};



typedef struct s_neighbor { 
public:
	double  dist;     // distance from the point datum
	int		datum;
}t_neighbor;
// for sorting elements wrt dist
struct Less_than {  
  bool operator()(const t_neighbor& a, const t_neighbor& b) {
      return a.dist < b.dist;  // based on last names only
  }
};

typedef struct s_neighborEx { 
public:
	double  dist;     // distance from the point datum
	int		datum;
	double	diff;
}t_neighborEx;
// for sorting elements wrt dist
struct Less_thanEx {  
  bool operator()(const t_neighborEx& a, const t_neighborEx& b) {
      return a.dist < b.dist;  // based on last names only
  }
};

/* auxiliary class OneRow
This class represents one row of the pairwise distance matrix. Besides the distances,
it also holds the information about the direction of the other data points wrt to this one.
As we use simplicial distance, there are Dim+1 possible directions.
The distance matrix is used to sort data points wrt to any of these directions.
It also keeps information about which data points have bigger j-th coordinate.
*/

class	OneRow {
public:
	CLargeSet**  m_sets;     // an array of sets of size Dim+1, each set holding up to Ndata members
	t_neighbor* m_neighbors; // the distances to this data point, and their indices live here
	int		Dim, Ndata;      //  dimension and the total number of data points
// used internally, but needs to be public for the Pack method
	int*	m_last, *m_first;// indices in the array m_neighbors indicating starting and ending indices for
							 // each direction

	OneRow(int dim, int ndata);
	~OneRow();

// add a distance d to the point J, in the direction dir (computed together with d)
	void	AddDistance(double d, int dir, int J);
// sort distances in the increasing order (in each direction)
	void	SortAll();

// returns the next closest data point, and its direction J. Automatically increments the index,
// so that the next call returns the subsequent data point. Returns a negative value if the list finished
	int		GetNextJ(double& d, int& J);
// should be called before the fist call to GetNextJ, resets the counter
	void	ResetCounter();

// removes all members of the set m_sets[dir] which are also members of the same set in the OtherRow
// does not affect the list of distances m_neighbors, just the elements of the set
	void	RemoveReferences(OneRow* OtherRow, int dir);

// Recomputes the sets m_sets from the list of distances m_neighbors
	void	ComputeSets();

// Packs a long array of distances in temprow into a shoter array m_neighbors of this instance
	void	Pack(OneRow* temprow);

// sets all the counters/sets to their initial values
	void	Reset();
private:
	int		counter,dir;
};

/*----------class OneRow-------------------------------------------------*/

OneRow::OneRow(int dim, int ndata)
{
	Dim=dim;
	Ndata=ndata;
	m_sets=(CLargeSet**) malloc(dim*sizeof(CLargeSet*));
	int i;
	for(i=0;i<dim;i++) m_sets[i]=new CLargeSet(ndata);

	m_neighbors=(t_neighbor*) malloc(/*dim*/ndata*sizeof(t_neighbor));
	m_last=(int*)malloc(sizeof(int)*dim);
	m_first=(int*)malloc(sizeof(int)*dim);
	for(i=0;i<dim;i++) m_last[i]=i*ndata;	
	for(i=0;i<dim;i++) m_first[i]=i*ndata;	
}

OneRow::~OneRow()
{
	int i;
	for(i=0;i<Dim;i++) delete(m_sets[i]);
	free(m_sets);
	if(m_neighbors!=NULL) free(m_neighbors);
	free(m_last);
	free(m_first);
}

void OneRow::AddDistance(double d, int dir, int J)
{
	m_neighbors[m_last[dir]].dist=d;
	m_neighbors[m_last[dir]].datum=J;
	m_sets[dir]->Set(J);
	m_last[dir]++;
}

// used only by one instance tempneighbor
void OneRow::Reset()
{	int i;
	for(i=0;i<Dim;i++) m_last[i]=i*Ndata/Dim;	 // this is for tempneigbor, careful for others
	for(i=0;i<Dim;i++) m_first[i]=i*Ndata/Dim;
	for(i=0;i<Dim;i++) m_sets[i]->Clear();
}

void OneRow::Pack(OneRow* temprow)
{
	int i,j,k;
	for(i=0;i<Dim;i++) m_sets[i]->Copy(temprow->m_sets[i]);
	k=0;
	for(i=0;i<Dim;i++) {
		m_first[i]=k;
		for(j=temprow->m_first[i]; j<temprow->m_last[i]; j++) {
			m_neighbors[k]=temprow->m_neighbors[j];
			k++;
		}
		m_last[i]=k;
	}
}

Less_than less_than;  /* declare a comparison function object, to
                          pass to sort and search algorithms */
Less_thanEx less_thanEx;  /* declare a comparison function object, to
                          pass to sort and search algorithms */
void OneRow::SortAll()
{
	for(int i=0;i<Dim;i++) 
		sort(&(m_neighbors[m_first[i]]),&(m_neighbors[m_last[i]]),less_than);
}

void	OneRow::ResetCounter()  {	counter=0; dir=0; }

int		OneRow::GetNextJ(double& d, int& J) // next after start
{
	int idx;
L1:	while(counter + m_first[dir] < m_last[dir]) {
		idx=m_neighbors[m_first[dir] + counter].datum;
		if(m_sets[dir]->IsPresent(idx)) 
		{
			d=m_neighbors[m_first[dir] + counter].dist;
			J=dir;
			counter++; // for the next time
			return idx;
		}
		else counter++;
	}

	if(dir < Dim-1) {
		dir++;
		counter=0;
		goto L1;
	}
	d=0; J=Dim;
	return -1;
}

//remove all indices from m_sets[dir] that are also in OtherRows's m_sets
void	OneRow::RemoveReferences(OneRow* OtherRow, int dir)
{ 
	m_sets[dir]->Remove(OtherRow->m_sets[dir]);
}

void	OneRow::ComputeSets()
{
	int i,j;
	for(i=0;i<Dim;i++) // direction
		for(j=m_first[i]; j < m_last[i]; j++) {
			m_sets[i]->Set( m_neighbors[j].datum );
		}
// I no longer need these, free memory
	if(m_neighbors!=NULL) free(m_neighbors);
	m_neighbors=NULL;
}


#define BMASK 0x3F;
#define SHIFT 6
void	CLargeSet::Set(int i)
{
	ComputeBit(i);
	m_els[pos] |= mask;
}

int	CLargeSet::SetChanged(int i)
{
	ComputeBit(i);
	i=((m_els[pos] & mask) != 0); // what it was
	m_els[pos] |= mask; // set new
	return i;
}

void	CLargeSet::Remove(int i)
{
	ComputeBit(i);
	m_els[pos] &= ~mask;
}

void	CLargeSet::Clear()
{
	for(int i=0;i<m_sizearray;i++) m_els[i]=0;
}

int		CLargeSet::IsPresent(int i)
{
	ComputeBit(i);
	return ((m_els[pos] & mask) != 0);
}

void	CLargeSet::ComputeBit(int i)
{
	pos = i >> SHIFT;
	mask = i & BMASK;
	mask = 1 << mask;
}

int		CLargeSet::NextElement(int i)
{ // this can be accelerated using bit masks
	i++;
	while(i< m_size)
		if(IsPresent(i)) return i;
		else i++;
	return -2; // finished, no more elements
}
	
// must be of the same size
void	CLargeSet::Remove(CLargeSet* OtherSet)
{
	for(int i=0;i<m_sizearray;i++) {
		m_els[i] &= ~(OtherSet->m_els[i]);
	}
}
// other set must be of larger or equal capacity!
void	CLargeSet::Copy(CLargeSet* OtherSet)
{
	for(int i=0;i<m_sizearray;i++) {
		m_els[i] = (OtherSet->m_els[i]);
	}
}


/*---------------------------end auxiliary classes---------------------------------------*/

// returns 1 if x >> y wrt Cons
int		SLipIntBasic::Dominates(int dim, double* x, double * y, int* Cons)
{
	int i;
	for(i=0;i<dim;i++) {
		if(Cons[i]!=0 && x[i]<y[i]) return 0;
	}
	return 1;
}


void SWAP(double* a, double* b, double *temp)
{
	*temp=*a;
	*a=*b;
	*b=*temp;
}

// Assumes the data in XData are stored columnwise, as in Fortran.
// Used for compatibility of this library with other packages, e.g.,Matlab, R, which
// may use column format. 
// Converting the data is faster than accessing it in F-format and converting the indices every time,
// as coordinates of one datum are not stored consequtively in memory
void	SLipIntBasic::ConvertXData(int dim, int npts,  double* XData)
{
	double temp;
	for(i=0;i<dim;i++)
		for(j=0;j<npts;j++) {
			SWAP(&(XData[i*dim+j]),&(XData[j*npts+i]),&temp);
		}		
}
void	SLipIntBasic::ConvertXData(int dim, int npts,  double* XData, double* auxStorage)
{
	for(i=0;i<dim;i++)
		for(j=0;j<npts;j++) {
			auxStorage[i*dim+j]=XData[j*npts+i];
		}		
}

double	SLipIntBasic::Value(int dim, int npts,   double* x, double* XData, double* YData, int* index)
{	return Value(dim,npts,x, XData, YData, MaxLipConst, index); }

double	SLipIntBasic::Value(int dim, int npts,   double* x, double* XData, double* YData,   double Lipconst, int* index)
{
	g1=-10e20;
	g2=-g1;

	if(index==NULL) {
		for(i=0;i<npts;i++) {
			j=i*dim;
			d2 = dist(dim,  x, &(XData[j]));
			d1= YData[i] - Lipconst * d2;
			d2= YData[i] + Lipconst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} else {
		for(i1=0;i1<npts;i1++) {
			i=index[i1];
			j=i*dim;
			d2 = dist(dim,  x, &(XData[j]));
			d1= YData[i] - Lipconst * d2;
			d2= YData[i] + Lipconst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	}
	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&Lipconst));
		g2=min__(g2,ExtraUpperBound(dim,x,&Lipconst));
	}
	return 0.5*(g1+g2);
}

int	SLipIntBasic::FindVoronoi(int dim, int npts,   double* x, double* XData, double &d)
{ // using DT
	int m,k;
	k=i=npts / 2;  // just some random starting point, not efficient in 1d 
	j=i*dim;
	d3 = d = dist(dim,  x, &(XData[j]));
	while(1) {
		for(i1=pneighbors[i];i1<pneighbors[i+1];i1++) {
			j=neighbors[i1];
			m=j*dim;
			d2 = dist(dim,  x, &(XData[m]));
			if(d3>d2) {k=j; d3=d2;}
		}
		if(d==d3) { // no improvement, finish
			return i;
		}
		i=k; // and continue
		d=d3;
	}
}

// Returns the value of the interpolant , with the Lipschitz constant
// computed from the data. Can be used after ComputeLocalLipschitz
double	SLipIntBasic::ValueLocal(int dim, int npts,   double* x, double* XData, double* YData)
{
	g1=-10e20;
	g2=-g1;

	int j1,k12;
	int k1;
	int lim;

	double dt;

	for(i=0;i<npts;i++) {
		k1=i*dim;
//		k12 = i*(npts-1);
		k12=GridLim[i+npts];

		d2 = dist(dim,  x, &(XData[k1]));
		lim=GridLim[i];	
		if(d2<GridR[k12]) {d1 =  d2/GridR[k12]*GridVal[k12];}
		else if(d2>=GridR[k12+lim-1]) d1 = GridVal[k12 + lim-1] + (GridVal[k12 + lim-1] - ((lim>1)?GridVal[k12 + lim -2]:0) ) * 
			(d2-GridR[k12 + lim -1]) / (GridR[k12 + lim -1]-((lim>1)?GridR[k12 + lim -2]:0));
		else {
			j1=BinSearch(d2, &(GridR[k12]), 0, GridLim[i]);
			d1=(d2-GridR[k12+j1])/(GridR[k12+j1+1]-GridR[k12+j1]);
			d1=d1*GridVal[k12+j1+1] + (1-d1)*GridVal[k12+j1];
		}

		dt= YData[i] - d1;		
		d1= YData[i] + d1;

		if(g1<dt) g1=dt;
		if(g2>d1) g2=d1;
	}
	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&MaxLipConst));
		g2=min__(g2,ExtraUpperBound(dim,x,&MaxLipConst));
	}

//	return g1;

	return 0.5*(g1+g2);



}



/********************************** Monotone approximation ****************************/


double	SLipIntBasic::ValueCons(int dim, int npts,  int* Cons, double* x, double* XData, double* YData)
{	return ValueCons(dim,npts, Cons, x, XData, YData, MaxLipConst); }

double	SLipIntBasic::ValueCons(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   double Lipconst, int* index)
{
	g1=-10e20;
	g2=-g1;

	if(index==NULL) {
		for(i=0;i<npts;i++) {
			j=i*dim;
			d2 = dist(dim,   &(XData[j]), x, Cons);
			d1= YData[i] - Lipconst * d2;

			d2 = dist(dim,  x, &(XData[j]), Cons);
			d2= YData[i] + Lipconst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} else {
		for(i1=0;i1<npts;i1++) {
			i=index[i1];
			j=i*dim;
			d2 = dist(dim,   &(XData[j]), x, Cons);
			d1= YData[i] - Lipconst * d2;

			d2 = dist(dim,  x, &(XData[j]), Cons);
			d2= YData[i] + Lipconst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	}

	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&Lipconst));
		g2=min__(g2,ExtraUpperBound(dim,x,&Lipconst));
	}

	return 0.5*(g1+g2);
}


int SLipIntBasic::BinSearch(double r, float* Arr, int le, int ri)
{
	int  l,u,mid;
	l=le; u=ri-1;
l2:
	if(u-l <=1) {  return l;}
	mid=(l+u)/2;
	if(r<(Arr)[mid]) u=mid; else l=mid;
	goto l2;

}

// Returns the value of the interpolant, with the Lipschitz constant
// computed from the data. Can be used after ComputeLocalLipschitz2
double	SLipIntBasic::ValueLocalCons(int dim, int npts, int* Cons,  double* x, double* XData, double* YData)
{
	return ValueLocal2Consinternal(dim,npts,Cons,x,XData, YData,0,NULL);
}

double	SLipIntBasic::ValueLocalConsLeftRegion(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, double* Region)
{	return ValueLocal2Consinternal(dim,npts,Cons,x,XData, YData,1,Region); }
double	SLipIntBasic::ValueLocalConsRightRegion(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, double* Region)
{	return ValueLocal2Consinternal(dim,npts,Cons,x,XData, YData,2,Region); }

// Returns the value of the interpolant in l_2 norm, with the Lipschitz constant
// computed from the data. Can be used after ComputeLocalLipschitz2
double	SLipIntBasic::ValueLocal2Consinternal(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, int reg, double* Region)
{
	g1=-10e20;
	g2=-g1;

	int j1,k12;
	int k1;
	int lim;

	double dt;

	for(i=0;i<npts;i++) {
		k1=i*dim;
//		k12 = i*(npts-1);
		k12=GridLim[i+npts];

		switch(reg) {
		case 1: d2=distLeftRegion(Dim,  &(XData[k1]), x, Cons, Region); break;
		case 2: d2=distRightRegion(Dim,  &(XData[k1]), x, Cons, Region);break;
		default:
		case 0: d2 = dist(dim,  &(XData[k1]),x, Cons);
			break;
		}
		
		lim=GridLim[i];	
		if(d2<GridR[k12]) {d1 =  d2/GridR[k12]*GridVal[k12];}

		else if(d2>=GridR[k12+lim-1]) d1 = GridVal[k12 + lim-1] + (GridVal[k12 + lim-1] - ((lim>1)?GridVal[k12 + lim -2]:0) ) * 
			(d2-GridR[k12 + lim -1]) / (GridR[k12 + lim -1]-((lim>1)?GridR[k12 + lim -2]:0));

//		else if(d2>=GridR[k12+GridLim[i]-1]) d1= GridVal[k12 + GridLim[i]-1];
		else {
			j1=BinSearch(d2, &(GridR[k12]), 0, GridLim[i]);
			d1=(d2-GridR[k12+j1])/(GridR[k12+j1+1]-GridR[k12+j1]);
			d1=d1*GridVal[k12+j1+1] + (1-d1)*GridVal[k12+j1];
		}

		dt= YData[i] - d1;

		switch(reg) {
		case 1: d2=distLeftRegion(Dim,  x, &(XData[k1]),  Cons, Region); break;
		case 2: d2=distRightRegion(Dim, x, &(XData[k1]),  Cons, Region);break;
		default:
		case 0: d2 = dist(dim, x, &(XData[k1]),Cons);
			break;
		}

		if(d2<GridR[k12]) {d1 =  d2/GridR[k12]*GridVal[k12];}
		else if(d2>=GridR[k12+lim-1]) d1 = GridVal[k12 + lim-1] + (GridVal[k12 + lim-1] - ((lim>1)?GridVal[k12 + lim -2]:0) ) * 
			(d2-GridR[k12 + lim -1]) / (GridR[k12 + lim -1]-((lim>1)?GridR[k12 + lim -2]:0));
//		if(d2<GridR[k12]) {d1 =  d2/GridR[k12]*GridVal[k12];}
//		else if(d2>=GridR[k12+GridLim[i]-1]) d1= GridVal[k12 + GridLim[i]-1];
		else {
			j1=BinSearch(d2, &(GridR[k12]), 0, GridLim[i]);
			d1=(d2-GridR[k12+j1])/(GridR[k12+j1+1]-GridR[k12+j1]);
			d1=d1*GridVal[k12+j1+1] + (1-d1)*GridVal[k12+j1];
		}
		
		d1= YData[i] + d1;

		if(g1<dt) g1=dt;
		if(g2>d1) g2=d1;
	}
	

	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&MaxLipConst));
		g2=min__(g2,ExtraUpperBound(dim,x,&MaxLipConst));
	}

	return 0.5*(g1+g2);
}


// Returns the value of the interpolant in l_2 norm, assuming it is monotone for x<< LeftRegion
double	SLipIntBasic::ValueConsLeftRegion(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* LeftRegion, int* index)
{
	g1=-10e20;
	g2=-g1;

	if(index==NULL) {
		for(i=0;i<npts;i++) {
			j=i*dim;
			d2 = distLeftRegion(dim,  &(XData[j]), x, Cons, LeftRegion);
			d1= YData[i] - LipConst * d2;

			d2 = distLeftRegion(dim, x, &(XData[j]), Cons, LeftRegion);
			d2= YData[i] + LipConst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} else {
		for(i1=0;i1<npts;i1++) {
			i=index[i1];
			j=i*dim;
			d2 = distLeftRegion(dim,  &(XData[j]), x, Cons, LeftRegion);
			d1= YData[i] - LipConst * d2;

			d2 = distLeftRegion(dim, x, &(XData[j]), Cons, LeftRegion);
			d2= YData[i] + LipConst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	}

	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&LipConst));
		g2=min__(g2,ExtraUpperBound(dim,x,&LipConst));
	}

	return 0.5*(g1+g2);

}

// Returns the value of the interpolant in l_2 norm,  assuming it is monotone for x>> RightRegion
double	SLipIntBasic::ValueConsRightRegion(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* RightRegion, int* index)
{
	g1=-10e20;
	g2=-g1;

	if(index==NULL) {

		for(i=0;i<npts;i++) {
			j=i*dim;
			d2 = distRightRegion(dim,  &(XData[j]), x, Cons, RightRegion);
			d1= YData[i] - LipConst * d2;

			d2 = distRightRegion(dim, x, &(XData[j]), Cons, RightRegion);
			d2= YData[i] + LipConst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} else {
		for(i1=0;i1<npts;i1++) {
			i=index[i1];
			d2 = distRightRegion(dim,  &(XData[j]), x, Cons, RightRegion);
			d1= YData[i] - LipConst * d2;

			d2 = distRightRegion(dim, x, &(XData[j]), Cons, RightRegion);
			d2= YData[i] + LipConst * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} 
	
	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,&LipConst));
		g2=min__(g2,ExtraUpperBound(dim,x,&LipConst));
	}

	return 0.5*(g1+g2);
}


// Verifies the data is monotone wrt specified variables
int		SLipIntBasic::VerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC, double eps)
{
	int i,j;
	for(i=0;i<npts;i++) {
		j=dim*i; 
		if(fabs( ValueCons( dim,  npts,   Cons, &(XData[j]),  XData,  YData,  LC) - YData[i]) > eps ) 
			{	return 0; }
	}
	return 1;
}

// Verifies the data is monotone wrt specified variables in the region On x<< LeftRegion 
int		SLipIntBasic::VerifyMonotonicityLeftRegion (int dim, int npts, int* Cons,  double* XData, double* YData, 
															 double* LeftRegion,double LC, double eps)
{
	int i,j;
	for(i=0;i<npts;i++) {
		j=dim*i; 
		if( Dominates(dim, LeftRegion, &(XData[j]) , Cons))  // only if the point is in the right region
		if(fabs( ValueConsLeftRegion( dim,  npts,   Cons, &(XData[j]),  XData,  YData, LC, LeftRegion) - YData[i]) > eps ) return 0;
	}
	return 1;
}
// Verifies the data is monotone wrt specified variables in the region On x >> RightRegion 
int		SLipIntBasic::VerifyMonotonicityRightRegion (int dim, int npts, int* Cons,  double* XData, double* YData, 
															  double* RightRegion, double LC, double eps)
{
	int i,j;
	for(i=0;i<npts;i++) {
		j=dim*i; 
		if( Dominates(dim, &(XData[j]) , RightRegion, Cons)) // only if the point is in the right region
		if(fabs( ValueConsRightRegion( dim,  npts,   Cons, &(XData[j]),  XData,  YData,  LC, RightRegion) - YData[i]) > eps ) return 0;
	}
	return 1;
}



double SLipIntBasic::value(int dim, int npts,   double* x, double* XData, double* YData,   double LipConst, int* index,
		 int type, int* Cons, double* Region)
// various parameters
{
	switch(type) {
	case 1: // monotonicity constraints
		return ValueCons(dim,npts,Cons,x,XData,YData,LipConst, index);
	case 2: // + left region
		ValueConsLeftRegion(dim,npts,Cons,x,XData,YData,LipConst,Region,index);
	case 3: // + Right region
		ValueConsRightRegion(dim,npts,Cons,x,XData,YData,LipConst,Region,index);
	case 0:
	default:
		return Value(dim,npts,x,XData,YData,LipConst,index);
	}
}

double SLipIntBasic::valuelocal(int dim, int npts,   double* x, double* XData, double* YData,
		 int type, int* Cons, double* Region)
{
	switch(type) {
	case 1: // monotonicity constraints
		return ValueLocalCons(dim,npts,Cons,x,XData,YData);
	case 2: // + left region
		ValueLocalConsLeftRegion(dim,npts,Cons,x,XData,YData,Region);
	case 3: // + Right region
		ValueLocalConsRightRegion(dim,npts,Cons,x,XData,YData,Region);
	case 0:
	default:
		return ValueLocal(dim,npts,x,XData,YData);
	}
}



/***** Implementations of the derived classes *******/


void	SLipInt::ComputeLipschitz(int dim, int npts, double* XData, double* YData)
{
	int k1,k2;
	MaxLipConst=0;
// for all pairs
	for(i=0;i<npts;i++) {
		k1=dim*i; 
		for(j=i+1;j<npts;j++) {
				k2=dim*j;
				d1=dist(dim,&(XData[k2]),&(XData[k1]));
				g1=fabs(YData[j] - YData[i]);

				if(d1>0) 
					MaxLipConst = max__(MaxLipConst, g1/d1);		
		}
	}
}


// Computes the local Lipschitz constant in l_2 norm, compatible with the data
void	SLipIntBasic::ComputeLocalLipschitz(int dim, int npts, double* XData, double* YData)
{
	if(GridR!=NULL) {
		if(npts!=NPTS) {
			free(GridR); free(GridVal); free(GridLim);
			GridLim=(int*)malloc(sizeof(int)*(npts*2+1));
		}
	} else {	GridLim=(int*)malloc(sizeof(int)*(npts*2+1));	}


	t_neighbor* m_neighbors = (t_neighbor*) malloc(sizeof(t_neighbor)*npts);

	vector<float>GridRT;
	vector<float>GridValT;

	GridRT.reserve(npts*3);
	GridValT.reserve(npts*3);

	MaxLipConst=0;

	NPTS=npts;

	int i,j,k, idx;
	int k1,k2,klim=0;
	double t1=0,t2=0;
	double r0,M0,Mt,w;
	double prevR, prevM;
	double prevprevR, prevprevVal, prevVal;

	GridLim[npts]=0;

	for(i=0;i<npts;i++) {
		// the i-th datum
		k1=i*dim;
		k=0;
		M0=0; Mt=0;
// start with the average Lip const for this datum
		for(j=0;j<npts;j++) if(i!=j) {
			k2=dim*j;
			t1 = dist(dim,&(XData[k1]),&(XData[k2]));
			m_neighbors[k].dist= t1;
			m_neighbors[k].datum= j;

			g1=fabs(YData[j] - YData[i]);
			w=1./t1;
			Mt += w;
			M0 += g1/t1 * w; //* weight?
			k++;
		}
		sort(&(m_neighbors[0]),&(m_neighbors[k]),less_than);
		//M0 /= k; 
		M0 /= Mt;
		r0=0;
// now find the grid and vals
		
		for(j=0;j<k;j++) {
			d1=m_neighbors[j].dist;
			g1=fabs(YData[i] - YData[m_neighbors[j].datum]);
			t2=g1/d1;
			if(t2>=M0) {
				if(r0==0) {r0=d1; t2 = g1; klim=j; goto L1;} 
			}
		}
L1:		if(r0==0)
		{
			r0=d1; // max dist
			idx=GridLim[npts + i ];
			GridRT.push_back((float)d1);
			GridValT.push_back((float)g1);
			GridLim[i]=1;
			goto L2; // end of the loop
		}


		idx=GridLim[npts + i ];
		GridRT.push_back((float)r0);
		GridValT.push_back((float)t2);
		M0=t2/r0;

		GridLim[i]=1;
		MaxLipConst=max__(M0,MaxLipConst);
		prevprevR=0; prevprevVal=0; prevR=r0; prevM=M0;
		prevVal=t2;

	// now all the rest
		for(j=klim+1;j<k;j++) {
			d1=m_neighbors[j].dist;
			t1 = fabs(YData[i] - YData[m_neighbors[j].datum]);
			if(d1>prevR) { 
				if(t1> (d1-prevR)*prevM + t2 ) {
					GridRT.push_back((float)d1);
					GridValT.push_back((float)t1);
					GridLim[i]++;
					prevprevR=prevR;
					prevprevVal=prevVal;
					 prevM = (t1-t2) / (d1-prevR); prevR=d1; prevVal=t1;

					t2=t1;
					MaxLipConst=max__(MaxLipConst, prevM);
				}} else { //d1==prevR, update is LipConst >
					if(t1>t2) {
						GridValT[GridValT.size()-1]=(float)t1;
						prevM = (t1-prevprevVal)/(d1-prevprevR);
						prevprevVal=prevVal;
						prevVal=t1;
						t2=t1;
						MaxLipConst=max__(MaxLipConst, prevM);
					}
				}

		} //j loop

L2:
		GridLim[npts+i+1]=GridLim[i]+GridLim[npts+i];

	}
// reallocate GridR and GridVal
	idx=GridLim[npts+npts];
	GridR=(float*) malloc(sizeof(float)*idx);
	GridVal=(float*) malloc(sizeof(float)*idx);
	for(i=0;i<idx;i++) {
		GridR[i]=GridRT[i];
		GridVal[i]=GridValT[i];
	}

	GridRT.clear(); GridValT.clear();
	free(m_neighbors);

}



// Computes the local Lipschitz constant in l_2 norm, compatible with the data
void	SLipIntBasic::ComputeLocalLipschitzCons(int dim, int npts, int _type, int* Cons, double* XData, double* YData, double* Region)
{
	if(GridR!=NULL) {
		if(npts!=NPTS) {
			free(GridR); free(GridVal); free(GridLim);
			GridLim=(int*)malloc(sizeof(int)*(npts*2+1));
		}
	} else {
			GridLim=(int*)malloc(sizeof(int)*(npts*2+1));
	}


	t_neighborEx* m_neighbors = (t_neighborEx*) malloc(sizeof(t_neighborEx)*npts);
//	t_neighbor* m_neighbors1 = (t_neighbor*) malloc(sizeof(t_neighbor)*npts);

	vector<float>GridRT;
	vector<float>GridValT;

	GridRT.reserve(npts*3);
	GridValT.reserve(npts*3);


	type=_type;
	LocalCons=Cons;
	LocalRegion=Region;
	MaxLipConst=0;

	NPTS=npts;

	int i,j,k, idx;
	int k1,k2,klim;
	double t1,t2;
	double r0,M0,Mt,w;
	double prevR, prevM;
	double prevprevR, prevprevVal, prevVal;


	GridLim[npts]=0;

	for(i=0;i<npts;i++) {
		// the i-th datum
		k1=i*dim;
		k=0;
		M0=0; Mt=0;
		for(j=0;j<npts;j++) if(i!=j) {
			k2=dim*j;
			t1 = distAll(dim,type,&(XData[k1]),&(XData[k2]), Cons,Region);
			t2 = distAll(dim,type,&(XData[k2]),&(XData[k1]), Cons,Region);
			g1=(YData[i] - YData[j]);

			if(g1>0 || t2==0) {
				m_neighbors[k].dist= t1;
				m_neighbors[k].datum= j;
				m_neighbors[k].diff = g1;
				w=1./t1;
				Mt += w;
				M0 += g1/t1 * w; //* weight?
			} else {
				m_neighbors[k].dist= t2;
				m_neighbors[k].datum= j;
				m_neighbors[k].diff = -g1;
				w=1./t2;
				Mt += w;
				M0 += -g1/t2 * w; //* weight?
			}

			k++;
		}
		sort(&(m_neighbors[0]),&(m_neighbors[k]),less_thanEx);
//		sort(&(m_neighbors1[0]),&(m_neighbors1[k]),less_than);

// now find the grid and vals
		M0 /= Mt;
		r0=0;
		for(j=0;j<k;j++) {
			d1=m_neighbors[j].dist;
			g1=m_neighbors[j].diff;
			if(g1>0) {
				t2=g1/d1;
				if(t2>=M0) {
					if(r0==0) {r0=d1; t2 = g1; klim=j; goto L1;} 
				}
			}
		}
L1:		// save them
		if(r0==0)
		{
			r0=d1; // max dist
			idx=GridLim[npts + i ];
			GridRT.push_back((float)d1);
			GridValT.push_back((float)g1);
			GridLim[i]=1;
			goto L3; // end of the loop
		}

		idx=GridLim[npts + i ];
		GridRT.push_back((float)r0);
		GridValT.push_back((float)t2);
		M0=t2/r0;

		GridLim[i]=1;
		MaxLipConst=max__(M0,MaxLipConst);
		prevprevR=0; prevprevVal=0; prevR=r0; prevM=M0;
		prevVal=t2;


	// now all the rest
		for(j=klim+1;j<k;j++) {
			d1 = m_neighbors[j].dist;
			t1 = m_neighbors[j].diff;
			if(d1>prevR) { 
				if(t1> (d1-prevR)*prevM + t2 ) {
					GridRT.push_back((float)d1);
					GridValT.push_back((float)t1);
					GridLim[i]++;
					prevprevR=prevR;
					prevprevVal=prevVal;
					 prevM = (t1-t2) / (d1-prevR); prevR=d1; prevVal=t1;

					t2=t1;
					MaxLipConst=max__(MaxLipConst, prevM);
				}} else { //d1==prevR, update is LipConst >
					if(t1>t2) {
						GridValT[GridValT.size()-1]=(float)t1;
						prevM = (t1-prevprevVal)/(d1-prevprevR);
						prevprevVal=prevVal;
						prevVal=t1;
						t2=t1;
						MaxLipConst=max__(MaxLipConst, prevM);
					}
				}

		} //j loop




L3:
		GridLim[npts+i+1]=GridLim[i]+GridLim[npts+i];

	}
// reallocate GridR and GridVal
	idx=GridLim[npts+npts];
	GridR=(float*) malloc(sizeof(float)*idx);
	GridVal=(float*) malloc(sizeof(float)*idx);
	for(i=0;i<idx;i++) {
		GridR[i]=GridRT[i];
		GridVal[i]=GridValT[i];
	}

	GridRT.clear(); GridValT.clear();
	free(m_neighbors);

}



// additions to version 2
int		SLipIntBasic::ComputeScaling(int dim, int npts, double* XData, double* YData)
{
//		ComputeLipschitzInf( dim,  npts,  XData,  YData);

		if(Scaling!=NULL) {
		if(dim!=Dim) {
			free(Scaling);
			Scaling=(double*) malloc(sizeof(double)*dim);
		}
		} else
			Scaling=(double*) malloc(sizeof(double)*dim);

		double* STD=(double*) malloc(sizeof(double)*dim);

		Dim=dim;
		int m;

		for(i=0;i<dim;i++) Scaling[i]=STD[i]=0;

		for(i=0;i<npts;i++) {
			m=i*dim;
			for(j=0;j<dim;j++) {
				Scaling[j]+=sqr__(XData[m+j]);
				STD[j]+=XData[m+j];
			}
		}

		for(j=0;j<dim;j++) Scaling[j] =sqrt( (Scaling[j]-sqr__(STD[j])/npts) / (npts-1));				

		for(i=0;i<dim;i++)
			Scaling[i] = 1.0/Scaling[i];

		free(STD);
		return 0;
}





/**************************************************************************
// methods below refer to Lipschitz smoothing

 ***************************************************************************/


// Smooth the data subject to given Lipschitz constant in Euclidean norm
void	SLipInt::SmoothLipschitz(int dim, int npts,  double* XData, double* YData, double* TData, double LC)
{	SmoothLipschitz2internal(dim,npts, XData,  YData,  TData, 0,0,0, &LC,  &LC, &Dim);}

// Smooth the data subject to given Lipschitz constant in Euclidean norm
void	SLipInt::SmoothLipschitzW(int dim, int npts,  double* XData, double* YData, double* TData, double LC, double* W)
{	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,1,0, &LC,  W,  &Dim);}



// same, subject to monotonicity constraints
void	SLipInt::SmoothLipschitzCons(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC)
{ 	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,0, 1, &LC,  &LC, Cons); }

void	SLipInt::SmoothLipschitzConsLeftRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* LeftRegion)
{ 	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,0, 1, &LC,  &LC, Cons, 1, LeftRegion); }

void	SLipInt::SmoothLipschitzConsRightRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* RightRegion)
{ 	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,0, 1, &LC,  &LC, Cons, 2, RightRegion); }

void	SLipInt::SmoothLipschitzWConsLeftRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC,double* W, double* LeftRegion)
{ 	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,1, 1, &LC,  W, Cons, 1, LeftRegion); }

void	SLipInt::SmoothLipschitzWConsRightRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* W, double* RightRegion)
{ 	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,1, 1, &LC,  W, Cons, 2, RightRegion); }


// same, subject to monotonicity constraints
// Smooth the data subject to given Lipschitz constant in Euclidean norm
void	SLipInt::SmoothLipschitzWCons(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double LC, double* W)
{	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,1,1, &LC,  W, Cons);}






void	SLipInt::ComputeLipschitzInf(int dim, int npts, double* XData, double* YData)
{
	if(LipConst!=NULL) {
		if(dim!=Dim) {
			free(LipConst);
			LipConst=(double*) malloc(sizeof(double)*dim);
		}
	} else
		LipConst=(double*) malloc(sizeof(double)*dim);

	for(j=0;j<dim;j++) LipConst[j]=0;

	Dim=dim;

	int dir;
	int k1,k2;

// for all pairs

	for(i=0;i<npts;i++) {
		k1=dim*i; 
		for(j=i+1;j<npts;j++) {
				k2=dim*j;
				d1=distInf1(dim,&(XData[k2]),&(XData[k1]),&dir);
				g1=fabs(YData[j] - YData[i]);

				if(d1>0) 
					LipConst[dir] = max__(LipConst[dir], g1/d1);
		}
	}

	MaxLipConst=0;
	for(i=0;i<dim;i++) MaxLipConst=max__(MaxLipConst,LipConst[i]);

}





/* ******************************** SLipIntInf class **********************************************/

// overwrites, uses the direction for directional Lip Constants
double	SLipIntInf::ValueDir(int dim, int npts,  double* x, double* XData, double* YData)
{	return ValueDir(dim,npts,x, XData, YData, LipConst); }

double	SLipIntInf::ValueDir(int dim, int npts, double* x, double* XData, double* YData,   double* Lipconst, int* index)
{ // overwrites, uses the direction for directional Lip Constants
	g1=-10e20;
	g2=-g1;

	if(index==NULL) {
		for(i=0;i<npts;i++) {
			j=i*dim;
			d2 = distDir(dim, x, &(XData[j]),&j);
			d1= YData[i] - Lipconst[j] * d2;
			d2= YData[i] + Lipconst[j] * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	} else {
		for(i1=0;i1<npts;i1++) {
			i=index[i1];
			j=i*dim;
			d2 = distDir(dim, x, &(XData[j]),&j);
			d1= YData[i] - Lipconst[j] * d2;
			d2= YData[i] + Lipconst[j] * d2;
			if(g1<d1) g1=d1;
			if(g2>d2) g2=d2;
		}
	}
	if(UseOtherBounds) {
		g1=max__(g1,ExtraLowerBound(dim,x,Lipconst));
		g2=min__(g2,ExtraUpperBound(dim,x,Lipconst));
	}
	return 0.5*(g1+g2);
}




void	SLipIntInf::ComputeLipschitz(int dim, int npts, double* XData, double* YData)
{
	if(LipConst!=NULL) {
		if(dim!=Dim) {
			free(LipConst);
			LipConst=(double*) malloc(sizeof(double)*dim);
		}
	} else
		LipConst=(double*) malloc(sizeof(double)*dim);

	for(j=0;j<dim;j++) LipConst[j]=0;

	Dim=dim;

	int dir;
	int k1,k2;

// for all pairs

	for(i=0;i<npts;i++) {
		k1=dim*i; 
		for(j=i+1;j<npts;j++) {
				k2=dim*j;
				d1=distDir(dim, &(XData[k2]),&(XData[k1]),&dir);
				g1=fabs(YData[j] - YData[i]);

				if(d1>0) 
					LipConst[dir] = max__(LipConst[dir], g1/d1);
		}
	}

	MaxLipConst=0;
	for(i=0;i<dim;i++) MaxLipConst=max__(MaxLipConst,LipConst[i]);

}







/**************************************************************************
// methods below refer to Lipschitz smoothing

 ***************************************************************************/


// Smooth the data subject to given Lipschitz constant in l-infy norm
void	SLipIntInf::SmoothLipschitz(int dim, int npts,  double* XData, double* YData, double* TData, double LC)
{	SmoothLipschitzInfinternal(dim,npts,XData,  YData,  TData, 0,0, &LC,  &LC); }


void	SLipIntInf::SmoothLipschitzW(int dim, int npts,  double* XData, double* YData, double* TData, double LC, double* W)
{	SmoothLipschitzInfinternal(dim,npts,XData,  YData,  TData, 0,1, &LC,  W); }


// same, subject to monotonicity constraints
void	SLipIntInf::SmoothLipschitzCons(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC)
{	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,0, 1, &LC,  &LC, Cons);}

// same, subject to monotonicity constraints
// Smooth the data subject to given Lipschitz constant 
void	SLipIntInf::SmoothLipschitzWCons(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double LC, double* W)
{	SmoothLipschitz2internal(dim,npts,XData,  YData,  TData, 0,1,1, &LC,  W, Cons);}


// Same in simplicial distance
void	SLipIntInf::SmoothLipschitzSimp(int dim, int npts,  double* XData, double* YData, double* TData,  double LC)
{	SmoothLipschitzSimpinternal(dim,npts,XData,  YData,  TData, 0, LC,  &LC); }
// Same in simplicial distance
void	SLipIntInf::SmoothLipschitzSimpW(int dim, int npts,  double* XData, double* YData, double* TData,  double LC, double* W)
{	SmoothLipschitzSimpinternal(dim,npts,XData,  YData,  TData, 1, LC,  W); }


//#define IPT

void	SLipIntInf::SmoothLipschitzInfinternal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,
											double* LocalLC, double* W, int* index)
{
	int i,j,k, iind, jind, direction ,it, ii,jj;
	int Ndata=npts;
	int Dim=dim;
	double m_Lip= *LocalLC;

	OneRow**	m_theneighbors;

	vector<int> Row,Col;
	vector<double> Vals, Obj;

	Vals.push_back(0);
	Row.push_back(0);
	Col.push_back(0);
	Obj.push_back(0);

	double d;
	m_theneighbors=(OneRow**) malloc(sizeof(OneRow*)*Ndata);
	for(i=0;i<Ndata;i++) m_theneighbors[i]=new OneRow(Dim*2,Ndata);  // memory needed is Ndata*Ndata*(1)

	// tempneighbor is used as a temp storage for all the distances, which are then
	// transferred to m_theneighbors in a compact form
	OneRow* tempneigbor=new OneRow(Dim*2,Ndata*(Dim*2));
	tempneigbor->Reset(); // reset all the indices

	for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
		i=(index==NULL) ? ii : index[ii];
		iind=i*dim;

		for(jj=ii+1;jj<Ndata;jj++) if (ii!=jj)   {// for every other datum
			j=(index==NULL) ? jj : index[jj];
			jind=j*dim;
			d=distInf2(Dim, &(XData[iind]),&(XData[jind]), &direction);
			tempneigbor->AddDistance(d,direction,jj);
		}  
		m_theneighbors[ii]->Pack(tempneigbor);
		tempneigbor->Reset();
	}

	delete tempneigbor;

	// all distances computed
	// sort every row
	for(i=0;i<Ndata;i++)  m_theneighbors[i]->SortAll();

	double temp1, temp2;
	k=1; j=-1;
	int dir;


#ifdef LPSOLVE
  MyLP = make_lp( 2*Ndata, 0);
  MyLP->do_presolve=FALSE;   
  set_verbose(MyLP,0);
#else
  MyLP = lpx_create_prob();

  lpx_add_rows(MyLP, 2*Ndata);
  lpx_set_int_parm(MyLP,LPX_K_MSGLEV,0);

  int Columns=0, cnt=0;

  lpx_set_obj_dir(MyLP, LPX_MIN);
#endif

	double row[5], row1[5];
	int    rowno[5];
	row[1]=1;
	row[2]=-1;
	row[3]=-1;
	row[4]=1;
	rowno[0]=0;
	row1[1]=-1;
	row1[2]=1;
	row1[3]=1;
	row1[4]=-1;

	for(ii=0;ii<Ndata;ii++) {
		m_theneighbors[ii]->ResetCounter();

		jj=m_theneighbors[ii]->GetNextJ(temp1,dir); // distance, direction
		i=(index==NULL) ? ii : index[ii];
		
		
		while(jj>=0) {
				m_theneighbors[ii]->RemoveReferences(m_theneighbors[jj],dir);
				j=(index==NULL) ? jj : index[jj];
				//temp1=dist(i,j);

				temp2=(YData[i]-YData[j]);

			if(LCf==1) 
				m_Lip=min__(LocalLC[i],LocalLC[j]);

// because the distance is symmetric
				// in columns
				row[0]=m_Lip*temp1 - temp2; // was RHS

				rowno[1]=ii+1;
				rowno[2]=jj+1;
				rowno[3]=ii+1+Ndata;
				rowno[4]=jj+1+Ndata;

				row1[0]=m_Lip*temp1 + temp2; // was RHS

#ifdef LPSOLVE
			add_columnex(MyLP, 5, row, rowno);
			add_columnex(MyLP, 5, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}
				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif

				k+=2;
			jj=m_theneighbors[ii]->GetNextJ(temp1,dir); // dist, direction
		} 
		m_theneighbors[ii]->ComputeSets(); // restore the sets for subseuent computations
		// I will never need the distances... remove them?
	}


// clear memory, everything is already in the hash table
	for (i=0;i<Ndata;i++) delete(m_theneighbors[i]);
	free(m_theneighbors);

	k--; // total columns
   m_number_constraints=k;


   if(UseOtherBounds) {
		for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
			i=(index==NULL) ? ii : index[ii];
			iind=i*dim;
			row[0]=    ExtraUpperBound(dim, &(XData[iind]), &m_Lip) - YData[i]; // Upper bound
			row1[0]= - ExtraLowerBound(dim, &(XData[iind]), &m_Lip) + YData[i]; // Lower bound
			rowno[1]=ii+1;
			rowno[2]=ii+1+Ndata;
#ifdef LPSOLVE
			add_columnex(MyLP, 3, row, rowno);
			add_columnex(MyLP, 3, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}
				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif
		}
   }




#ifdef LPSOLVE
#else
   lpx_add_cols(MyLP, Columns); // how many
   lpx_load_matrix(MyLP, cnt, &Row[0], &Col[0], &Vals[0]);//begin(), Col.begin(), Vals.begin());
   for(it=1;it<=Columns;it++) {
	   lpx_set_obj_coef(MyLP, it, Obj[it]);
	   lpx_set_col_bnds(MyLP, it, LPX_LO, 0.0,0.0);
   }
#endif

   Row.clear();
   Col.clear();
   Vals.clear();
   Obj.clear();


// rhs
#ifdef LPSOLVE
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for u
		set_constr_type(MyLP,ii+1,GE);
	}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for v
		set_constr_type(MyLP,ii+1,GE);}
   } else
	for(i=0;i<Ndata*2;i++) {	
		set_rh(MyLP,i+1, -1 ); // all weights the same
		set_constr_type(MyLP,i+1,GE);}

	
	set_minim(MyLP); // well, we always do that
#else
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
			}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
		}
   } else
	for(i=0;i<Ndata*2;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  -1.0, 0.0);
		 // all weights the same
	}
#endif

	m_lasterror=0;
	int		res;
	double	minval;
	double  u,v;

#ifdef LPSOLVE
	int RR=get_Nrows(MyLP);
	int CC=get_Ncolumns(MyLP);
	double *sol=(double*)malloc(sizeof(double)*(1 + RR + CC));

	res=solve(MyLP);

#else
#ifdef IPT
	res=lpx_interior(MyLP);
#else
	res=lpx_simplex(MyLP);
#endif
#endif

#ifdef LPSOLVE
	if(res==OPTIMAL) {
#ifdef DUAL
		get_dual_solution(MyLP, sol);
#else
		get_primal_solution(MyLP, sol);
#endif
		m_minvalue = minval = -get_objective(MyLP) ;  // minimum

		for(i=1;i<=Ndata;i++)
		{
#ifdef DUAL
			u= sol[i] ;
			v= sol[i+Ndata] ;
#else
			u= sol[RR+i] ;
			v= sol[RR+i+Ndata] ;
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		
			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 
	free (sol);

#else

	if(res==LPX_E_OK) {
#ifdef IPT
		m_minvalue = minval = -lpx_ipt_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#else
		m_minvalue = minval = -lpx_get_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#endif

		for(i=1;i<=Ndata;i++)
		{
#ifdef IPT
			u=lpx_ipt_row_dual(MyLP,i);
			v=lpx_ipt_row_dual(MyLP,i+Ndata);
#else
			u=lpx_get_row_dual(MyLP,i);
			v=lpx_get_row_dual(MyLP,i+Ndata);
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		

			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 

	lpx_delete_prob(MyLP);

#endif
}




void	SLipIntBasic::SmoothLipschitz2internal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf, int Cf,
											double* LocalLC, double* W, int* Cons, int region, double* Region, int* index)

{
	int k, iind, jind ,it,ii,jj;
	int Ndata=npts;
	int Dim=dim;
	double m_Lip = *LocalLC;


	double d,d1;

	double  temp2;

	k=1; j=-1;

	vector<int> Row,Col;
	vector<double> Vals, Obj;

	Vals.reserve(4*Ndata*(Ndata-1)+1);
	Row.reserve(4*Ndata*(Ndata-1)+1);
	Col.reserve(4*Ndata*(Ndata-1)+1);
	Obj.reserve(Ndata*Ndata+1);

	Vals.push_back(0);
	Row.push_back(0);
	Col.push_back(0);
	Obj.push_back(0);

#ifdef LPSOLVE
	if(KeepCVProblem==3)	{ delete_lp(MyLP);KeepCVProblem=0;}

  MyLP = make_lp( 2*Ndata, 0);
  MyLP->do_presolve=FALSE;   
  set_verbose(MyLP,0); 
#else
  if(KeepCVProblem==3)	{ lpx_delete_prob(MyLP);KeepCVProblem=0;} 

  MyLP = lpx_create_prob();

  lpx_add_rows(MyLP, 2*Ndata);
  lpx_set_int_parm(MyLP,LPX_K_MSGLEV,0);

  int Columns=0, cnt=0;

  lpx_set_obj_dir(MyLP, LPX_MIN);
#endif


	double row[5], row1[5];
	int    rowno[5];
// all pairs

	row[1]=1;
	row[2]=-1;
	row[3]=-1;
	row[4]=1;
	rowno[0]=0;
	row1[1]=-1;
	row1[2]=1;
	row1[3]=1;
	row1[4]=-1;

	for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
		i=(index==NULL) ? ii : index[ii];
		iind=i*dim;

		for(jj=ii+1;jj<Ndata;jj++) //if (ii!=jj)  // because const. distance is asymetric
		{ // for every other datum
			j=(index==NULL) ? jj : index[jj];
			jind=j*dim;

			temp2=(YData[i]-YData[j]);

			if(Cf==1) {
				if(region==0)
				{
					d= dist(Dim,  &(XData[iind]), &(XData[jind]),Cons);
					d1=dist(Dim,  &(XData[jind]), &(XData[iind]),Cons);
				} else if (region==1) // left
				{
					d=distLeftRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1=distLeftRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons, Region);
				} else // right
				{
					d = distRightRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1 = distRightRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons,Region);
				}
			}
			else
				d=d1=dist(Dim,  &(XData[iind]),&(XData[jind]));

			if(LCf==1) 
				m_Lip=min__(LocalLC[i],LocalLC[j]);


			row[0]=m_Lip*d - temp2; // was RHS
			row1[0]=m_Lip*d1 + temp2; // was RHS

			rowno[1]=ii+1;
			rowno[2]=jj+1;
			rowno[3]=ii+1+Ndata;
			rowno[4]=jj+1+Ndata;

#ifdef LPSOLVE
			add_columnex(MyLP, 5, row, rowno);
			add_columnex(MyLP, 5, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}
				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif
			k+=2;
		}  
	}

	k--; // total columns
   m_number_constraints=k;

   if(UseOtherBounds) {
		for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
			i=(index==NULL) ? ii : index[ii];
			iind=i*dim;
			row[0]=    ExtraUpperBound(dim, &(XData[iind]), &m_Lip) - YData[i]; // Upper bound
			row1[0]= - ExtraLowerBound(dim, &(XData[iind]), &m_Lip) + YData[i]; // Lower bound
			rowno[1]=ii+1;
			rowno[2]=ii+1+Ndata;
#ifdef LPSOLVE
			add_columnex(MyLP, 3, row, rowno);
			add_columnex(MyLP, 3, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}
				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif
		}
   }



#ifdef LPSOLVE
#else
   lpx_add_cols(MyLP, Columns); // how many

   lpx_load_matrix(MyLP, cnt, &Row[0], &Col[0], &Vals[0]);
   for(it=1;it<=Columns;it++) {
	   lpx_set_obj_coef(MyLP, it, Obj[it]);
	   lpx_set_col_bnds(MyLP, it, LPX_LO, 0.0,0.0);
   }
#endif

   Row.clear();
   Col.clear();
   Vals.clear();
   Obj.clear();

// rhs
#ifdef LPSOLVE
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for u
		set_constr_type(MyLP,ii+1,GE);
	}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for v
		set_constr_type(MyLP,ii+1,GE);}
   } else
	for(i=0;i<Ndata*2;i++) {	
		set_rh(MyLP,i+1, -1 ); // all weights the same
		set_constr_type(MyLP,i+1,GE);}

	
	set_minim(MyLP); // well, we always do that
#else
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
			}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
		}
   } else
	for(i=0;i<Ndata*2;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  -1.0, 0.0);
		 // all weights the same
	}
#endif

	m_lasterror=0;

	int		res;
	double	minval;
	double  u,v;

#ifdef LPSOLVE
	int RR=get_Nrows(MyLP);
	int CC=get_Ncolumns(MyLP);
	double *sol=(double*)malloc(sizeof(double)*(1 + RR + CC));
	res=solve(MyLP);

#else

#ifdef IPT
	res=lpx_interior(MyLP);
#else
	res=lpx_simplex(MyLP);
#endif
#endif


#ifdef LPSOLVE
	if(res==OPTIMAL) {
#ifdef DUAL
		get_dual_solution(MyLP, sol);
#else
		get_primal_solution(MyLP, sol);
#endif
		m_minvalue = minval = -get_objective(MyLP) ;  // minimum

		for(i=1;i<=Ndata;i++)
		{
#ifdef DUAL
			u= sol[i] ;
			v= sol[i+Ndata] ;
#else
			u= sol[RR+i] ;
			v= sol[RR+i+Ndata] ;
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		
			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 
	free (sol);
	if(KeepCVProblem==0) delete_lp(MyLP);

#else

	if(res==LPX_E_OK) {
#ifdef IPT
		m_minvalue = minval = -lpx_ipt_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#else
		m_minvalue = minval = -lpx_get_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#endif

		for(i=1;i<=Ndata;i++)
		{
#ifdef IPT
			u=lpx_ipt_row_dual(MyLP,i);
			v=lpx_ipt_row_dual(MyLP,i+Ndata);
#else
			u=lpx_get_row_dual(MyLP,i);
			v=lpx_get_row_dual(MyLP,i+Ndata);
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		

			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 

	if(KeepCVProblem==0) lpx_delete_prob(MyLP);
#endif

}

// TODO for lp-solve
void	SLipIntBasic::SmoothLipschitz2internalUpdate(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf, int Cf,
											double* LocalLC, double* W, int* Cons, int region, double* Region, int* index)

{ // we only need to change RHS and Obj function
	int i,j,k, iind, jind ,it,ii,jj;
	int Ndata=npts;
	int Dim=dim;
	double m_Lip = *LocalLC;

	double d,d1;

	double  temp2;

	k=1; j=-1;

	vector<double>  Obj;

	Obj.reserve(Ndata*Ndata+2);

	Obj.push_back(0);
	int Columns=0;

// all pairs

	for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
		i=(index==NULL) ? ii : index[ii];
		iind=i*dim;

		for(jj=ii+1;jj<Ndata;jj++) //if (i!=j) 
		{ // for every other datum
			j=(index==NULL) ? jj : index[jj];
			jind=j*dim;

			temp2=(YData[i]-YData[j]);

			if(Cf==1) {
				if(region==0)
				{
					d= dist(Dim,  &(XData[iind]), &(XData[jind]),Cons);
					d1=dist(Dim,  &(XData[jind]), &(XData[iind]),Cons);
				} else if (region==1) // left
				{
					d=distLeftRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1=distLeftRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons, Region);
				} else // right
				{
					d = distRightRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1 = distRightRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons,Region);
				}
			}
			else
				d=d1=dist(Dim, &(XData[iind]),&(XData[jind]));

			if(LCf==1) 
				m_Lip=min__(LocalLC[i],LocalLC[j]);

			Columns++;
			Obj.push_back(m_Lip*d - temp2);

			Columns++;
			Obj.push_back(m_Lip*d1 + temp2);
		}  
	}


   if(UseOtherBounds) 
		for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
			i=(index==NULL) ? ii : index[ii];
			iind=i*dim;
			Columns++;
			Obj.push_back(ExtraUpperBound(dim, &(XData[iind]), &m_Lip) - YData[i]);
			Columns++;
			Obj.push_back(- ExtraLowerBound(dim, &(XData[iind]), &m_Lip) + YData[i]);
   }



#ifdef LPSOLVE
	set_obj_fn(MyLP,&Obj[0]);
#else
   for(it=1;it<=Columns;it++) {
	   lpx_set_obj_coef(MyLP, it, Obj[it]);
	   lpx_set_col_bnds(MyLP, it, LPX_LO, 0.0,0.0);
   }
#endif

   Obj.clear();

// rhs - keep this as weights might change
#ifdef LPSOLVE
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for u
		set_constr_type(MyLP,ii+1,GE);
	}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for v
		set_constr_type(MyLP,ii+1,GE);}
   } else
	for(i=0;i<Ndata*2;i++) {	
		set_rh(MyLP,i+1, -1 ); // all weights the same
		set_constr_type(MyLP,i+1,GE);}
	
	set_minim(MyLP); // well, we always do that
#else
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
			}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
		}
   } else
	for(i=0;i<Ndata*2;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  -1.0, 0.0);
		 // all weights the same
	}
#endif


	m_lasterror=0;

	int		res;
	double	minval;
	double  u,v;

#ifdef LPSOLVE

	int RR=get_Nrows(MyLP);
	int CC=get_Ncolumns(MyLP);
	double *sol=(double*)malloc(sizeof(double)*(1 + RR + CC));

	res=solve(MyLP);

#else  // glpk

	lpx_warm_up(MyLP);
#ifdef IPT
	res=lpx_interior(MyLP);
#else
	res=lpx_simplex(MyLP);
#endif
#endif


#ifdef LPSOLVE
	if(res==OPTIMAL) {
#ifdef DUAL
		get_dual_solution(MyLP, sol);
#else
		get_primal_solution(MyLP, sol);
#endif
		m_minvalue = minval = -get_objective(MyLP) ;  // minimum

		for(i=1;i<=Ndata;i++)
		{
#ifdef DUAL
			u= sol[i] ;
			v= sol[i+Ndata] ;
#else
			u= sol[RR+i] ;
			v= sol[RR+i+Ndata] ;
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		
			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 
	free (sol);

#else

	if(res==LPX_E_OK) {
#ifdef IPT
		m_minvalue = minval = -lpx_ipt_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#else
		m_minvalue = minval = -lpx_get_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#endif

		for(i=1;i<=Ndata;i++)
		{
#ifdef IPT
			u=lpx_ipt_row_dual(MyLP,i);
			v=lpx_ipt_row_dual(MyLP,i+Ndata);
#else
			u=lpx_get_row_dual(MyLP,i);
			v=lpx_get_row_dual(MyLP,i+Ndata);
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		

			TData[ii] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 

#endif
}





void	SLipIntInf::SmoothLipschitzSimpinternal(int dim, int npts,  double* XData, double* YData, double* TData,  int Wf,
											double LC, double* W, int* index)
{
	int i,j,k, iind, jind, direction, it,ii,jj ;
	int Ndata=npts;
	int Dim=dim;
	double m_Lip=LC;
	
	OneRow**	m_theneighbors;

	vector<int> Row,Col;
	vector<double> Vals, Obj;

	Vals.push_back(0);
	Row.push_back(0);
	Col.push_back(0);
	Obj.push_back(0);

	double d;
	m_theneighbors=(OneRow**) malloc(sizeof(OneRow*)*Ndata);
	for(i=0;i<Ndata;i++) m_theneighbors[i]=new OneRow(Dim+1,Ndata);  // memory needed is Ndata*Ndata*(1)

	// tempneighbor is used as a temp storage for all the distances, which are then
	// transferred to m_theneighbors in a compact form
	OneRow* tempneigbor=new OneRow(Dim+1,Ndata*(Dim+1));
	tempneigbor->Reset(); // reset all the indices

	for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
		i=(index==NULL) ? ii : index[ii];
		iind=i*dim;

		for(jj=ii+1;jj<Ndata;jj++) if (ii!=jj) { // for every other datum
			j=(index==NULL) ? jj : index[jj];
			jind=j*dim;
			d=distSimp(Dim,  &(XData[iind]),&(XData[jind]), &direction);
//			m_theneighbors[i]->AddDistance(d,direction,j);
			tempneigbor->AddDistance(d,direction,jj);
		}  
		m_theneighbors[ii]->Pack(tempneigbor);
		tempneigbor->Reset();
	}

	delete tempneigbor;

	// all distances computed
	// sort every row
	for(i=0;i<Ndata;i++)  m_theneighbors[i]->SortAll();

	// use setrowmode !!! hash s memory hungry

	double temp1, temp2;
	k=1; j=-1;
	int dir;

#ifdef LPSOLVE
  MyLP = make_lp( 2*Ndata, 0);
  MyLP->do_presolve=FALSE;   
  set_verbose(MyLP,0);
#else
  MyLP = lpx_create_prob();

  lpx_add_rows(MyLP, 2*Ndata);
  lpx_set_int_parm(MyLP,LPX_K_MSGLEV,0);

  int Columns=0, cnt=0;

  lpx_set_obj_dir(MyLP, LPX_MIN);
#endif

	double row[5], row1[5];
	int    rowno[5];
// all pairs

	row[1]=1; 
	row[2]=-1;
	row[3]=-1;
	row[4]=1;
	rowno[0]=0;
	row1[1]=-1;
	row1[2]=1;
	row1[3]=1;
	row1[4]=-1;

	for(ii=0;ii<Ndata;ii++) {
		m_theneighbors[ii]->ResetCounter();

		jj=m_theneighbors[ii]->GetNextJ(temp1,dir); // distance, direction
		i=(index==NULL) ? ii : index[ii];			
		
		while(jj>=0) {
				m_theneighbors[ii]->RemoveReferences(m_theneighbors[jj],dir);
				j=(index==NULL) ? jj : index[jj];
				temp2=(YData[i]-YData[j]);

				// in columns
				row[0]=m_Lip*temp1 - temp2; // was RHS
				row1[0]=m_Lip*temp1 + temp2; // was RHS

				rowno[1]=ii+1;
				rowno[2]=jj+1;
				rowno[3]=ii+1+Ndata;
				rowno[4]=jj+1+Ndata;

#ifdef LPSOLVE
			add_columnex(MyLP, 5, row, rowno);
			add_columnex(MyLP, 5, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}

				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=4;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif
			k+=2;

			jj=m_theneighbors[ii]->GetNextJ(temp1,dir); // dist, direction
		} 
		m_theneighbors[ii]->ComputeSets(); // restore the sets for subseuent computations
		// I will never need the distances... remove them?
	}

// clear memory, everything is already in the hash table
	for (i=0;i<Ndata;i++) delete(m_theneighbors[i]);
	free(m_theneighbors);

	k--; // total columns
   m_number_constraints=k;


   if(UseOtherBounds) {
		for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
			i=(index==NULL) ? ii : index[ii];
			iind=i*dim;
			row[0]=    ExtraUpperBound(dim, &(XData[iind]), &m_Lip) - YData[i]; // Upper bound
			row1[0]= - ExtraLowerBound(dim, &(XData[iind]), &m_Lip) + YData[i]; // Lower bound
			rowno[1]=ii+1;
			rowno[2]=ii+1+Ndata;
#ifdef LPSOLVE
			add_columnex(MyLP, 3, row, rowno);
			add_columnex(MyLP, 3, row1, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
					cnt++;
				}
				Columns++;
				Obj.push_back(row1[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row1[it]);
					cnt++;
				}
#endif
		}
   }



#ifdef LPSOLVE
#else
   lpx_add_cols(MyLP, Columns); // how many
   lpx_load_matrix(MyLP, cnt, &Row[0], &Col[0], &Vals[0]);
   for(it=1;it<=Columns;it++) {
	   lpx_set_obj_coef(MyLP, it, Obj[it]);
	   lpx_set_col_bnds(MyLP, it, LPX_LO, 0.0,0.0);
   }
#endif

   Row.clear();
   Col.clear();
   Vals.clear();
   Obj.clear();

// rhs
#ifdef LPSOLVE
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for u
		set_constr_type(MyLP,ii+1,GE);
	}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for v
		set_constr_type(MyLP,ii+1,GE);}
   } else
	for(i=0;i<Ndata*2;i++) {	
		set_rh(MyLP,i+1, -1 ); // all weights the same
		set_constr_type(MyLP,i+1,GE);}

	
	set_minim(MyLP); // well, we always do that
#else
   if(Wf==1) {
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];		
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
			}
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii - Ndata: index[ii - Ndata];
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], 0.0);
		}
   } else
	for(i=0;i<Ndata*2;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  -1.0, 0.0);
		 // all weights the same
	}
	lpx_set_obj_dir(MyLP, LPX_MIN);
#endif

	m_lasterror=0;

	int		res;
	double	minval;
	double  u,v;

#ifdef LPSOLVE
	double *sol=(double*)malloc(sizeof(double)*(1 + get_Nrows(MyLP) + get_Ncolumns(MyLP)));

	int RR=get_Nrows(MyLP);
	int CC=get_Ncolumns(MyLP);
	res=solve(MyLP);

#else
#ifdef IPT
	res=lpx_interior(MyLP);
#else
	res=lpx_simplex(MyLP);
#endif
#endif

#ifdef LPSOLVE
	if(res==OPTIMAL) {
#ifdef DUAL
		get_dual_solution(MyLP, sol);
#else
		get_primal_solution(MyLP, sol);
#endif
		m_minvalue = minval = -get_objective(MyLP) ;  // minimum

		for(i=1;i<=Ndata;i++)
		{
#ifdef DUAL
			u= sol[i] ;
			v= sol[i+Ndata] ;
#else
			u= sol[RR+i] ;
			v= sol[RR+i+Ndata] ;
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		
			TData[i-1] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 
	free (sol);

#else

	if(res==LPX_E_OK) {
#ifdef IPT
		m_minvalue = minval = -lpx_ipt_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#else
		m_minvalue = minval = -lpx_get_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#endif

		for(i=1;i<=Ndata;i++)
		{
#ifdef IPT
			u=lpx_ipt_row_dual(MyLP,i);
			v=lpx_ipt_row_dual(MyLP,i+Ndata);
#else
			u=lpx_get_row_dual(MyLP,i);
			v=lpx_get_row_dual(MyLP,i+Ndata);
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		

			TData[i-1] = YData[ii] + (u-v);
		}
	} else m_lasterror=res; // LP not solved 

	lpx_delete_prob(MyLP);

#endif

}
















/* ***************************  Lipschitz Classifier *********************************
class SLipClass: public SLipInt {

// Returns the value of the Classifier and performs smoothingwith a special objective fucntion
};
*/


int		SLipClass::ValueConsClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   double LipConst)
{ 	if(ValueCons( dim,  npts,   Cons,  x,  XData,  YData,    LipConst) >= 0) return +1; else return -1; }

int		SLipClass::ValueClass(int dim, int npts, double* x, double* XData, double* YData,   double LipConst)
{	if(Value( dim,  npts,   x,  XData,  YData,    LipConst) >= 0) return +1; else return -1; }

int		SLipClass::ValueConsLeftRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* LeftRegion)
{	if(ValueConsLeftRegion( dim,  npts,   Cons,  x,  XData,  YData,    LipConst, LeftRegion) >= 0) return +1; else return -1; }

int		SLipClass::ValueConsRightRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* RightRegion)
{	if(ValueConsRightRegion( dim,  npts,   Cons,  x,  XData,  YData,    LipConst, RightRegion) >= 0) return +1; else return -1;}

int		SLipClass::ValueLocalClass(int dim, int npts, double* x, double* XData, double* YData)
{	if(ValueLocal( dim,  npts,  x,  XData,  YData) >= 0) return +1; else return -1;}

int		SLipClass::ValueLocalConsClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData)
{	if(ValueLocalCons( dim,  npts,   Cons,  x,  XData,  YData) >= 0) return +1; else return -1;}

int		SLipClass::ValueLocalConsLeftRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
												   double* LeftRegion)
{	if(ValueLocalConsLeftRegion( dim,  npts,   Cons,  x,  XData,  YData,   LeftRegion) >= 0) return +1; else return -1;}

int		SLipClass::ValueLocalConsRightRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
												  double* RightRegion)
{	if(ValueLocalConsRightRegion( dim,  npts,   Cons,  x,  XData,  YData,   RightRegion) >= 0) return +1; else return -1;}


void	SLipClass::SmoothLipschitzClass(int dim, int npts,  double* XData, double* YData, double* TData, double *LC)
{	SmoothLipschitz2Classinternal(dim,npts, XData,  YData,  TData, 0,0,0, LC,  LC, NULL);}

void	SLipClass::SmoothLipschitzWClass(int dim, int npts,  double* XData, double* YData, double* TData, double *LC, double* W)
{	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,1,0, LC,  W,  NULL);}

void	SLipClass::SmoothLipschitzConsClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC)
{ 	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,0, 1, LC,  LC, Cons); }

void	SLipClass::SmoothLipschitzWConsClass(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double *LC, double* W)
{	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,1,1, LC,  W, Cons);}

void	SLipClass::SmoothLipschitzConsLeftRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* LeftRegion)
{ 	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,0, 1, LC,  LC, Cons, 1, LeftRegion); }

void	SLipClass::SmoothLipschitzConsRightRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* RightRegion)
{ 	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,0, 1, LC,  LC, Cons, 2, RightRegion); }

void	SLipClass::SmoothLipschitzWConsLeftRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* W, double* LeftRegion)
{ 	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,1, 1, LC,  W, Cons, 1, LeftRegion); }

void	SLipClass::SmoothLipschitzWConsRightRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* W, double* RightRegion)
{ 	SmoothLipschitz2Classinternal(dim,npts,XData,  YData,  TData, 0,1, 1, LC,  W, Cons, 2, RightRegion); }



int sign(double a) {return (a>=0)? 1:-1;}
void	SLipClass::SmoothLipschitz2Classinternal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,int Cf,
											double* LocalLC, double* W, int* Cons, int region, double* Region, int* index)
{// assumes YData are + or -1
	int i,j,k, iind, jind ,it, ii,jj;
	int Ndata=npts;
	int Dim=dim;
	double m_Lip = *LocalLC;

	if(SmoothingParam>0) m_Lip=0; // will find it myself

	double d,d1;
	double  temp2;

	k=1; j=-1;

#ifdef LPSOLVE
  if(SmoothingParam==0)
	  MyLP = make_lp( 2*Ndata, 0); 
  else 
	  MyLP = make_lp( 2*Ndata+1, 0); 
  MyLP->do_presolve=FALSE;   
  set_verbose(MyLP,0);

#else

  MyLP = lpx_create_prob();

  if(SmoothingParam==0)
	 lpx_add_rows(MyLP, 2*Ndata);
  else
	 lpx_add_rows(MyLP, 2*Ndata+1);

  lpx_set_int_parm(MyLP,LPX_K_MSGLEV,0);

  int Columns=0;

  lpx_set_obj_dir(MyLP, LPX_MIN);
#endif

	double row[4], row1[4],rowT[4], row1T[4];
	int    rowno[4];
	rowT[1]=row[1]=1;  // jst cache these values
	rowT[2]=row[2]=-1;
	row1T[1]=row1[1]=-1;
	row1T[2]=row1[2]=1;
	rowno[0]=0;

	{
	vector<int> Row,Col;
	vector<double> Vals, Obj, Lastrow;

	Vals.reserve(2*Ndata*(Ndata+1)+1);
	Row.reserve(2*Ndata*(Ndata+1)+1);
	Col.reserve(2*Ndata*(Ndata+1)+1);
	Obj.reserve(2*Ndata*(Ndata+1)+1);
	Lastrow.reserve(2*Ndata*Ndata+1);

	Vals.push_back(0);
	Row.push_back(0);
	Col.push_back(0);
	Obj.push_back(0);
	Lastrow.push_back(0);
// all pairs

	for(ii=0;ii<Ndata;ii++) { // for every datum (row) // parallelisation can be done here
		i=(index==NULL) ? ii : index[ii];
		iind=i*dim;

		for(jj=ii+1;jj<Ndata;jj++) //if (ii!=jj) 
		{ // for every other datum
			j= (index==NULL) ? jj : index[jj];
			jind=j*dim; 

			temp2=(YData[i]-YData[j]);

			if(Cf==1) {
				if(region==0)
				{
					d= dist(Dim,  &(XData[iind]), &(XData[jind]),Cons);
					d1=dist(Dim, &(XData[jind]), &(XData[iind]),Cons);
				} else if (region==1) // left
				{
					d=distLeftRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1=distLeftRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons, Region);
				} else // right
				{
					d = distRightRegion(Dim,  &(XData[iind]), &(XData[jind]), Cons, Region);
					d1 = distRightRegion(Dim,  &(XData[jind]), &(XData[iind]), Cons,Region);
				}
			}
			else
				d=d1=dist(Dim, &(XData[iind]),&(XData[jind]));

			if(LCf==1) 
				m_Lip=min__(LocalLC[i],LocalLC[j]);

			row[0]=m_Lip*d - temp2; // was RHS
			row1[0]=m_Lip*d1 + temp2; // was RHS
			rowno[1]=ii+1;
			rowno[2]=jj+1;
#ifdef LPSOLVE
			row[1] = rowT[1]* -sign(YData[i]) ;
			row[2] = rowT[2]* -sign(YData[j]) ;
			row1[1] = row1T[1]* -sign(YData[i]) ;
			row1[2] = row1T[2]* -sign(YData[j]) ;
			if(SmoothingParam==0){
				add_columnex(MyLP, 3, row, rowno);
				add_columnex(MyLP, 3, row1, rowno);
			}
			else { //extra row
				row[3]=-d; row1[3]=-d1;
				rowno[3]=2*Ndata+1;
				add_columnex(MyLP, 4, row, rowno);
				add_columnex(MyLP, 4, row1, rowno);
			}
#else
				Columns++;
				Obj.push_back(row[0]);

					Row.push_back(rowno[1]);
					Col.push_back(Columns);
					Vals.push_back(row[1]* -sign(YData[i]) );
					Row.push_back(rowno[2]);
					Col.push_back(Columns);
					Vals.push_back(row[2]* -sign(YData[j]) );

				Lastrow.push_back(d);

				Columns++;
				Obj.push_back(row1[0]);
					Row.push_back(rowno[1]);
					Col.push_back(Columns);
					Vals.push_back(row1[1]*-sign(YData[i]) );
					Row.push_back(rowno[2]);
					Col.push_back(Columns);
					Vals.push_back(row1[2]*-sign(YData[j]) );
				Lastrow.push_back(d1);
#endif
			k+=2;
		}  
	}


#ifdef LPSOLVE
#else
   int ColTem=Columns;  // to know how many nonzeroes
   if(SmoothingParam>0) {
	   k=2*Ndata+1; // last row
	   for(i=1;i<=Columns;i++) {
		   Row.push_back(k);
		   Col.push_back(i);
		   Vals.push_back(-Lastrow[i]);
	   }
   }
#endif
   
   row[0]=0;
   row[2]=-1;
   rowno[0]=0;
   row[1]=1;
   for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];
			//row[1]=1*(YData[i]>=0? -1 : +1);;
			rowno[1]=ii+1;
			rowno[2]=ii+1+Ndata;
				k++;
#ifdef LPSOLVE
				add_columnex(MyLP, 3, row, rowno);
#else
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
				}
#endif
	}
	row[0]=Penalty; // was RHS
	row[2]=-1;
	row[1]=(Penalty+1);
	rowno[0]=0;
	for(ii=0;ii<Ndata;ii++) {
		i=(index==NULL) ? ii : index[ii];
//			row[1]=(Penalty+1)*(YData[i]>=0? -1 : +1);
			rowno[1]=ii+1;
			rowno[2]=ii+1+Ndata;
				k++;
#ifdef LPSOLVE
				add_columnex(MyLP, 3, row, rowno);
#else				
				Columns++;
				Obj.push_back(row[0]);
				for(it=1;it<=2;it++) {
					Row.push_back(rowno[it]);
					Col.push_back(Columns);
					Vals.push_back(row[it]);
				}
#endif
	}

	k--; // total columns

   m_number_constraints=k;

#ifdef LPSOLVE
#else
   lpx_add_cols(MyLP, Columns); // how many

   if(SmoothingParam>0) 
	   lpx_load_matrix(MyLP, 2 * Columns+ColTem, &Row[0], &Col[0], &Vals[0]);
	else
	   lpx_load_matrix(MyLP, 2 * Columns, &Row[0], &Col[0], &Vals[0]);

   for(it=1;it<=Columns;it++) {
	   lpx_set_obj_coef(MyLP, it, Obj[it]);
	   lpx_set_col_bnds(MyLP, it, LPX_LO, 0.0,0.0);  //>=0
   }
#endif
   Row.clear();
   Col.clear();
   Vals.clear();
   Obj.clear();
   Lastrow.clear();
}

// rhs

#ifdef LPSOLVE
   if(Wf==1) {
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii-Ndata : index[ii-Ndata];	
		set_rh(MyLP,ii+1, -W[i] ); // we can put weights here for u
		set_constr_type(MyLP,ii+1,GE); //????????????
		}
   } else
	for(i=Ndata;i<Ndata*2;i++) {	
		set_rh(MyLP,i+1, -1.0 ); // we can put weights here for v
		set_constr_type(MyLP,i+1,GE); //????????????
		 // all weights the same
	}

	if(SmoothingParam>0) {
		set_rh(MyLP,Ndata*2+1, -Ndata*SmoothingParam ); // we can put weights here for v
		set_constr_type(MyLP,Ndata*2+1,GE); //????????????
	}
	
	set_minim(MyLP); // well, we always do that
#else
	for(i=0;i<Ndata;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  0.0, 0.0);
	}


   if(Wf==1) {
	for(ii=Ndata;ii<Ndata*2;ii++) {	
		i=(index==NULL) ? ii-Ndata : index[ii-Ndata];		
		lpx_set_row_bnds(MyLP, ii+1, LPX_LO,  -W[i], -W[i]);
		}
   } else
	for(i=Ndata;i<Ndata*2;i++) {	
		lpx_set_row_bnds(MyLP, i+1, LPX_LO,  -1.0, -1.0);
		 // all weights the same
	}

   if(SmoothingParam>0) 
		lpx_set_row_bnds(MyLP, Ndata*2+1, LPX_LO,  -Ndata*SmoothingParam, 0.0);
#endif

	m_lasterror=0;

	int		res;
	double	minval;
	double  u,v;


#ifdef LPSOLVE
	int RR=get_Nrows(MyLP);
	int CC=get_Ncolumns(MyLP);
	double *sol=(double*)malloc(sizeof(double)*(1 + RR + CC));
	res=solve(MyLP);

#else  // glpk
	lpx_warm_up(MyLP);
#ifdef IPT
	res=lpx_interior(MyLP);
#else
	res=lpx_simplex(MyLP);
#endif
#endif


#ifdef LPSOLVE
	if(res==OPTIMAL) {
#ifdef DUAL
		get_dual_solution(MyLP, sol);
#else
		get_primal_solution(MyLP, sol);
#endif
		m_minvalue = minval = -get_objective(MyLP) ;  // minimum

		for(i=1;i<=Ndata;i++)
		{
#ifdef DUAL
			u= sol[i] ;
#else
			u= sol[RR+i] ;
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		
			TData[ii] = YData[ii] + (u)*-sign(YData[ii]);

		if(SmoothingParam>0) { // Optimal Lipschitz constant?
			LocalLC[0] = sol[2*Ndata]; 
		}

		}
	} else m_lasterror=res; // LP not solved 
	free (sol);
// glpk
#else
	if(res==LPX_E_OK) {
#ifdef IPT
		m_minvalue = minval = -lpx_ipt_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#else
		m_minvalue = minval = -lpx_get_obj_val(MyLP);// -get_objective(MyLP) ;  // minimum
#endif

		for(i=1;i<=Ndata;i++)
		{
#ifdef IPT
			u=lpx_ipt_row_dual(MyLP,i);
#else
			u=lpx_get_row_dual(MyLP,i);
#endif
			ii=(index==NULL) ? i-1 : index[i-1];		

			TData[ii] = YData[ii] + (u)*-sign(YData[ii]);
		}

		if(SmoothingParam>0) { // Optimal Lipschitz constant?
#ifdef IPT
			LocalLC[0] = lpx_ipt_row_dual(MyLP,2*Ndata+1);
#else
			LocalLC[0] = lpx_get_row_dual(MyLP,2*Ndata+1);
#endif
		}

	} else m_lasterror=res; // LP not solved 

	lpx_delete_prob(MyLP);
#endif
}







/* This function implements the golden section search
   algorithm. 
   
     Input: A,B: the interval ends, B > A.
     Output: The argument that minimizes f()

   The makes successive calls to f() until it has found
   the minimum with the desired accuracy
*/
double SLipIntBasic::golden(double A, double B) 
{
  double alf1=0, alf2=0;
  double falf1=0, falf2=0;
  //double G=0.618034;  /* The golden section */
  double   G=(-1.+sqrt(5.))/2.; /* The golden section */
  double tol;


  int  N1=80; //Number of Fibonacci search steps
  int  N2=17; //Number of Golden section search steps
  
  double increment;
  int i,sign =0;
  double q0,q1,q2,alfA,alfB;
  double delta=0.05;  //factor for determining the initial interval of uncertainty
	

  /* Make sure B > A */
  if (B<=A) {
  //  fprintf(stderr,"golden: error: Illegal interval ends\n");
    return A;
  }
 /*Fibonacci search*/
  q0=A;
  q1= q0 + delta;
  increment=delta*(G+1); 
  q2=q1+increment;
 
  falf1=Fun(q1);
  if (falf1>(Fun(q0))) {
    //printf("Error: please choose a new starting point A");
    //return 0;
	  A=q0; B=q1;
	  goto Gold;
  }
  
  for (i=0;i<N1;i++){
	  falf2=Fun(q2);
    if (falf2 < falf1){
      q0=q1;
      q1=q2;
		falf1=falf2;
      increment=increment*(G+1);
      q2=q2+increment;
      
    }
    else { //Find the interval bracketing a minimum [q0,q2]
      alfA=q0;
      alfB=q2;
      sign=1;   // successful sign
      break;
    } 
  }
  if (sign==0) { // unsuccessful search
  //  printf("Error: please increase the cycle number N! to find a valid initial interval");
  }
  
 
  A=alfA;  //New interval found through the fibonacci search
  B=alfB;   
 
Gold:
  /* Iinitialize variables */
  tol=(B-A)*0.0001;
  alf1 = A + (1-G)*(B-A);
  alf2 = B - (1-G)*(B-A);
  falf1=(Fun)(alf1);
  falf2=(Fun)(alf2);

  for (i =0;i<N2;i++){
    
    // Use the left hand interval, if the function value at the
    // right hand golden point is the larger
    if (falf2>falf1) {

      
      // Shift re-usable results left
      B = alf2;
      alf2 = alf1;
      falf2 = falf1;
      
      // Compute new alf1 and function value
      alf1 = A + (1-G)*(B-A);
      falf1 = (Fun)(alf1);
    }
    
    // otherwise, use the right hand interval
    else {
      // Shift re-usable results left
      A = alf1;
      alf1 = alf2;
      falf1 = falf2;
      
      // Compute new Alpha2 and function value
      alf2 = B - (1-G)*(B-A);
      falf2 = (Fun)(alf2);
    }
  }  // Golden section loop

  /* Return the midpoint of the interval when it is small enough */
  return (alf1+alf2)/2;

}  /* golden */




double SLipIntBasic::Fun(double x)
{
	switch(	TypeLipEstimate) {
	case 1: return MinFuncCV(x);
	case 2: return MinFuncLocalSplit(x);
	case 0:
	default: return MinFuncSplit(x);
	}
}

double SLipIntBasic::MinFuncSplit(double x)
{
	M=x; // set up LipConst
	ComputeSmoothenedSplit();

	return ComputeFitIndex();
}

double SLipIntBasic::MinFuncCV(double x)
{
	M=x; // set up LipConst
	return ComputeFitIndexCV();
}

double SLipIntBasic::MinFuncLocalSplit(double x)
{
	M=x; // set up LipConst
	return 0 ;
	//	ComputeSmoothenedLocalSplit();
}



/************** The methods below implement sample splitting and cross-validation ***********/


int SLipIntBasic::ComputeSmoothenedSplit()
{
// called by Problem.fv()
// interpretation of the flags
	int Wf= (LocalW==0?0:1);
	int Tf=(type>0?1:0);
	int region=0;
	if(type ==2) region=1; 
	else if(type==3) region=2;

	SmoothLipschitz2internal(Dim, Indexsize,  LocalXData, LocalYData, LocalTData, 0,  Wf, Tf,
							&M, LocalW, LocalCons, region, LocalRegion, Index);

	return 0;
}


int SLipIntBasic::ComputeFitLipschitzCV(int excluded)
{
// called by golden section

	int Wf= (LocalW==0?0:1);
	int Tf=(type>0?1:0);
	int region=0;
	if(type ==2) region=1; 
	else if(type==3) region=2;

	int iL,iI;

// set up index arrays
	iI=0;
	for(iL=0;iL<NPTS;iL++) {
		if(iL<excluded) {Index[iI]=iL; iI++;}
		else if(iL==excluded) IndexComp[0]=iL;
		else  {Index[iI]=iL; iI++;} //(i>excluded)
	}

	if(KeepCVProblem==1) {
		SmoothLipschitz2internal(Dim, Indexsize,  LocalXData, LocalYData, LocalTData, 0,  Wf, Tf,
							&M, LocalW, LocalCons, region, LocalRegion, Index);
		KeepCVProblem=2;
	} else if (KeepCVProblem==2)
		SmoothLipschitz2internalUpdate(Dim, Indexsize,  LocalXData, LocalYData, LocalTData, 0,  Wf, Tf,
							&M, LocalW, LocalCons, region, LocalRegion, Index);
	return 0;
}



int		SLipIntBasic::ComputeLipschitzFinal()
// called once the optimal Lipschitz constant has been found. Smoothen the whole data set.
{
	int Wf= (LocalW==0?0:1);
	int Tf=(type>0?1:0);
	int region=0;
	if(type ==2)
		region=1; 
	else if(type==3) region=2;

	SmoothLipschitz2internal(Dim, NPTS,  LocalXData, LocalYData, LocalTData, 0,  Wf, Tf,
							&M, LocalW, LocalCons, region, LocalRegion,NULL);
	return 0;
}



int SLipIntInf::ComputeSmoothenedSplit()
{
	if(Dim>=5)  return SLipIntBasic::ComputeSmoothenedSplit();

	int Wf= (LocalW==0?0:1);

	SmoothLipschitzInfinternal(Dim, Indexsize,  LocalXData, LocalYData, LocalTData, 0,  Wf,
							&M, LocalW, Index);
	return 0;
}


int SLipIntInf::ComputeFitLipschitzCV(int excluded)
{
	if(Dim>=5)  return SLipIntBasic::ComputeFitLipschitzCV(excluded);

	int Wf= (LocalW==0?0:1);

	int iL,iI;

// set up index arrays
	iI=0;
	for(iL=0;iL<NPTS;iL++) {
		if(iL<excluded) {Index[iI]=iL; iI++;}
		else if(iL==excluded) IndexComp[0]=iL;
		else  {Index[iI]=iL; iI++;} //(i>excluded)
	}

	SmoothLipschitzInfinternal(Dim, Indexsize,  LocalXData, LocalYData, LocalTData, 0,  Wf, 
							&M, LocalW,   Index);
	return 0;
}



int		SLipIntInf::ComputeLipschitzFinal()
// called once the optimal Lipschitz constant has been found. Smoothen the whole data set.
{
	int Wf= (LocalW==0?0:1);

	if(Dim < 5)
		SmoothLipschitzInfinternal(Dim, NPTS,  LocalXData, LocalYData, LocalTData, 0,  Wf, 
							&M, LocalW, NULL);
	else
		return SLipIntBasic::ComputeLipschitzFinal();

	return 0;
}



/********* Generic methods for CV and sample splitting ******/

double SLipIntBasic::ComputeFitIndex()
{//	compute goodness of fit
	int i,idx;
	double r1,r=0;
	for(i=0;i< IndexsizeComp; i++) {
		idx=IndexComp[i]*Dim;
		r1=value(Dim, Indexsize,   &(LocalXData[idx]), LocalXData, LocalTData,   M, Index, type, LocalCons, LocalRegion);

		r+= sqr__(LocalYData[IndexComp[i]] - r1);
	}
	return r;
}

double SLipIntBasic::ComputeFitIndexCV()
{//	compute goodness of fit using CV, called from golden section
	int i,idx, excl;
	double r1,r=0;

	for(excl=0;excl<NPTS;excl++) {

		ComputeFitLipschitzCV(excl); // prepares Index, sets up LP, solves it

		for(i=0;i< IndexsizeComp; i++) {
			idx=IndexComp[i]*Dim;
			r1=value(Dim, Indexsize,   &(LocalXData[idx]), LocalXData, LocalTData,   M, Index, type, LocalCons, LocalRegion);

			r+= sqr__(LocalYData[IndexComp[i]] - r1);
		}
	}

	return r;
}



// entry points to CV and split: these are generic methods

void	SLipIntBasic::ComputeLipschitzSplit(int dim, int npts, double* XData, double* YData, double* TData, double ratio,
			int type, int* Cons, double* Region, double *W)
{
	// entry point to Split
	// set up split ratio
	Dim=dim;
	NPTS=npts;
	LocalXData=XData;
	LocalYData=YData;
	LocalTData=TData;
	LocalW=W;
	LocalRegion=Region;
	LocalCons=Cons;
	SLipIntBasic::type=type;

	PrepareLipschitzSplit(ratio);
	// declare problem and quasi

	// compute max Lip.Const
	ComputeLipschitz(dim,npts,XData, YData);
	TypeLipEstimate=0; //0 sample splitting, 1 CV

	M = MaxLipConst = golden(0,MaxLipConst); // golden section + Fibonnacci search

	ComputeLipschitzFinal();

	free(Index);
	free(IndexComp);
}

void	SLipIntBasic::ComputeLipschitzCV(int dim, int npts, double* XData, double* YData, double* TData,
			int type, int* Cons, double* Region, double *W)
{
// entry point to CV

	Dim=dim;
	NPTS=npts;
	LocalXData=XData;
	LocalYData=YData;
	LocalTData=TData;
	LocalW=W;
	LocalRegion=Region;
	LocalCons=Cons;
	SLipIntBasic::type=type;

	PrepareLipschitzCV();

	KeepCVProblem=1;

	ComputeLipschitz(dim,npts,XData, YData);
	TypeLipEstimate=1; //0 sample splitting, 1 CV

	M = MaxLipConst = golden(0,MaxLipConst);

	KeepCVProblem=3; //desroys LP first
	ComputeLipschitzFinal();

	KeepCVProblem=0; //back to default value

	free(Index);
	free(IndexComp);
}


// returns 1 with probability p and 0 with 1-p
int RandomBin(double p)
{
	int i=rand();
	if(i<=RAND_MAX*p) return 1; else return 0;
}


void SLipIntBasic::PrepareLipschitzSplit(double SplitP)
{
// prepares the arrays and splits randomly the sample
	int i,j,k;

	Indexsize= (int)(SplitP * NPTS); // or something
	ND=Indexsize;
	IndexsizeComp=NPTS-Indexsize;

	Index=(int*) malloc(Indexsize*sizeof(int));
	IndexComp=(int*) malloc(IndexsizeComp*sizeof(int));
	
	for(i=0;i<Indexsize;i++) Index[i]=0;
	for(i=0;i<IndexsizeComp;i++) IndexComp[i]=0;

	// set  to 1 2 5 7 8 - ie the elements that are in
	j=k=0;
	for(i=0;i<NPTS;i++)
	{
		if((j<Indexsize) && (k>=IndexsizeComp || RandomBin(SplitP)) ) { Index[j]=i; j++;}
		else { IndexComp[k]=i; k++; }
	}
}

void SLipIntBasic::PrepareLipschitzCV()
{
// prepares the arrays for CV
	Indexsize= NPTS - 1; 
	ND=Indexsize;
	IndexsizeComp=NPTS-Indexsize;

	Index=(int*) malloc(Indexsize*sizeof(int));
	IndexComp=(int*) malloc(IndexsizeComp*sizeof(int));
	
	for(i=0;i<Indexsize;i++) Index[i]=0;
	for(i=0;i<IndexsizeComp;i++) IndexComp[i]=0;
}



