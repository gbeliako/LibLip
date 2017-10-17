/**************************************************************************

 ***************************************************************************/

#include "interpol.h"

// Gobal variables

MemoryBlock<SVSetNode>	MBSV;
Interpolant*	Parent;


// merges 2 doubles into one, for lexicographic ordering
// this helps to reduce the complexity of testing for repeated vector components
double merge(double a, double b)
{
	double c;
	ULINT  d=(d2ulint(a) & ULLMASK) | ((d2ulint(b)>>32) & ~ULLMASK);  // fixed in version 1.2
	c=ulint2d(d);
	return c;
}


// aux global function, to avoid duplicating data points
int TestPresent(support_vector & sv,SVDeque& seq, int Dim)
{
//	return 0;
	int i,j;
	int k;
	support_vector* svp;
	for(i=Dim;i<seq.size(); i++) {   // skip the first Dim reserved vectors
		k=1;
		svp=&(seq[i]);
		for(j=sv.vec.size()-2;j>0;j--)  // in the opposite direction
			if(sv.vec[j] != svp->vec[j]) 
					{ k=0; break;}
		if(k) return 1;
	}
	return 0;
}
// some preprocessing resulting in a faster version: check first whether there is a vector
// with first component the same as that of the tested vector, using binary search in a sorted array
// and then test the second component, using binary search, and only then test all other components
// using linear search
int TestPresent(support_vector & sv, vector<double> &first, vector<double> &second, SVDeque& seq, int dim)
{
//	double* location = lower_bound(first.begin(), first.end(), sv.vec[0]) ;
// g++ does not like it
	double* location = lower_bound(&first[0], &first[first.size()], sv.vec[0]) ;

	if(*location == *(location+1)) {
		if(dim>1) {
//			location = lower_bound(second.begin(), second.end(), merge(sv.vec[0],sv.vec[1])) ;
			location = lower_bound(&second[0], &second[second.size()], merge(sv.vec[0],sv.vec[1])) ;
			if(*location == *(location+1)) 
				return TestPresent(sv,seq,dim+1);
			else return 0;
		} else
			return 1; // only one component, already see that it is duplicated
	}
		
// not duplicated
	return 0;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Interpolant::Interpolant()
{
	Parent=this;
	C_minus=1;
	m_Constants.newsize(1);
	m_Constants=0;
}; // constructor


Interpolant::~Interpolant()
{	
	HeapPossibleMin.EraseAll();
}

void Interpolant::FreeMem()
{	
	HeapPossibleMin.EraseAll();
	SVectors.clear();
	LastLabel=0;
}

void	Interpolant::Init(int dim)
{
	Dim=dim;
	m_lastindex.newsize(Dim);
	m_index.newsize(Dim);
	m_initvec.newsize(Dim);
   m_TempSiVec.newsize(Dim);

	Parent=this;
}


void	Interpolant::QueryDyn(real* x)
{
	Parent=this;

	m_sv.SVForm(x,0);
	HeapPossibleMin.ProcessAllDyn(&m_sv);
	m_indexset=&(HeapPossibleMin.m_indexset);
}

/*--------------------------------------------------
	Calculates the value of the lower interpolant
	using fast query method. If query unsuccessful,
	uses slower and safer explicit evaluation.

----------------------------------------------------*/
real	Interpolant::FValueL(real* x)
{
	int i;
	real g1,g1m,t1;
	Match=0;
	// search for the triangles
	QueryDyn(x);
	g1m=-Infinity; 
	for(i=0;i<Dim;i++)	m_indexset->erase(i); // erase basis vectors

	//  here explicit search
	if(m_indexset->empty()) { 
	//	cout<< "explicit"<<x<<endl;
		return FValueExplicit(x);
	}

	support_vector* SVP;
	m_indexiter=m_indexset->begin();
//	cout << m_indexset->size()<<endl;

	while(m_indexiter!=m_indexset->end())
	{
		g1=Infinity;
		SVP = &(SVectors[*m_indexiter]);
			for(i=0;i<Dim;i++) 
			{ // every component
				t1 = SVP->vec[i];

				t1=m_Constants[i] * (t1 + x[i]); 
			
				if(t1<g1) g1=t1;	 // min
				if(g1<=g1m) break; // no point to continue, already smaller than g1m (max is needed)
			}
			if(g1m<g1) { g1m=g1; } // max
		m_indexiter++; 	
	}
	return g1m ;
}





real	Interpolant::FValueExplicit(real* x)
{
	// exhaustive evaluation
	int i,j;
	real g1,g1m,t1;

	Match=0;

	g1m=-Infinity; 

	j=LastLabel;
	while(j>=Dim)  // >=Dim?
	{		j--;
			g1=Infinity;
			for(i=0;i<Dim;i++) 
			{ // every component
				t1 = SVectors[j].vec[i];

				t1=m_Constants[i] * (t1 + x[i]); 
			
				if(t1<g1) g1=t1;	 // min
				if(g1<=g1m) break; // no point to continue, already smaller than g1m (max is needed)
			}
			if(g1m<g1) { g1m=g1; } // max
	}
	return g1m ;
}


void Interpolant::Construct()
{
	Match=0;
	Init(Dim);

	ComputeCminus();
	HeapPossibleMin.Init();

	InitPopulateBoundaryPoints();

// build the tree			
	LoadAdditionalPoints();
}

void Interpolant::ConstructInv()
{
	for(int i=0;i < SVectors.size();i++) 	SVectors[i].funvalue = -SVectors[i].funvalue; // just invert the data

	Construct();
}

void Interpolant::ConstructExplicit()
{
	Match=0;
	Init(Dim);

	ComputeCminus();
//	HeapPossibleMin.Init();
	for(int i=Dim;i<SVectors.size();i++)  
		SVectors[i].SVForm(LastLabel);

	LastLabel=SVectors.size();
}

void Interpolant::ConstructInvExplicit()
{
	for(int i=0;i < SVectors.size();i++) 	SVectors[i].funvalue = -SVectors[i].funvalue; // just invert the data

	ConstructExplicit();
}



void	Interpolant::InitPopulateBoundaryPoints() 
{
	LastLabel=0;
	temp_x.newsize(Dim);

	temp_x=0.0;
	int i,j;
	support_vector	sv;

	for(i=0;i<Dim;i++)
	{
		SVectors[i].vec[i]= SInfinity;
		for(j=0;j<Dim;j++)
			if(i != j)
				SVectors[i].vec[j] = (1.-SInfinity)/(Dim-1.);
		SVectors[i].SVForm(SVectors[i].vec,SVectors[i].funvalue);
		SVectors[i].label=LastLabel;
		LastLabel++;
	}

	// create the root node of the tree

	// now generate new SVset, with every support vector we just constructed
	SVSetNodePtr n=MBSV.GetNextFree();
	SVSetNode*	S=MBSV.GetAt(n);
	S->Init();
		//new SVSetNode;

	for(i=0;i<Dim;i++)	HeapPossibleMin.temp_index[i]=i;

	HeapPossibleMin.AddRootNode(n);

	LastLabel=Dim;
}


void	Interpolant::LoadAdditionalPoints()
{
	// load points and function values, 
	unsigned int i;
	int sz1, sz2;

	for(i=Dim;i<SVectors.size();i++)  {
		SVectors[i].SVForm(LastLabel);

		LastLabel++;

		HeapPossibleMin.ProcessAll(&(SVectors[i]));

/*		would do log file if required
		sz1=1;
		sz2=HeapPossibleMin.GetSize() ;
		if(logfile!=NULL)
		{fprintf(logfile," %d %d %d %d %f %f %f %f\n",SVectors.size(),sz1,sz2,0,
				ElapsedTime(),0,0,0); //#endif
			fflush(logfile);
		}
*/

	}	

}

void	Interpolant::ComputeCminus()
{
	Parent=this; // just in case for the interpolant
	int i,j;
// just to ensure right consts
	if(m_Constants.size()!=Dim) {
		dVec temp=m_Constants;
		m_Constants.newsize(Dim);
		j=min_(temp.size(),Dim);
		for( i=0;i<j;i++)
		if(temp[i]>0)
			m_Constants[i]=temp[i];
		else m_Constants[i]=1;
		for(i=j;i<Dim;i++) m_Constants[i]=m_Constants[j-1];
	}

	C_minus=0;
	for( i=0;i<Dim;i++) C_minus += 1.0/m_Constants[i];
}

void	Interpolant::SetConstants(dVec& newconst) {m_Constants=newconst;}
void	Interpolant::SetConstants(real newconst, int dim) {m_Constants.newsize(dim); m_Constants=newconst; m_Constants[dim-1] *=sqrt(dim-1.0);}



/* --------------------------Implementation of STCInterpolant class -------------------------------- */
STCInterpolant::STCInterpolant()
{	m_lower=new Interpolant;
    m_upper=new Interpolant;
	Dim=0;
	Lip1=Lip2=aux=0;
}

STCInterpolant::~STCInterpolant()
{	
	delete m_lower;
	delete m_upper;
	Dim=0;  
	if(!Lip1) free(Lip1);
	if(!Lip2) free(Lip2);
	if(!aux) free(aux);
}

void	STCInterpolant::SetConstants(real newconst, int dim)
{
	m_lower->SetConstants(newconst,dim);
	m_upper->SetConstants(newconst,dim);
}
void	STCInterpolant::SetConstants()
{	SetConstants(LipschitzConst,Dim); }

void	STCInterpolant::SetConstants(real newconst)
{
	if(Dim>0)
		SetConstants(newconst,Dim);
	else SetConstants(newconst, 1);
}

void	STCInterpolant::SetData(int dim, int K, real* x, real* y, int test)
{
	m_lasterr=0;
	Dim=dim+1;

	m_lower->Init(Dim);
	m_upper->Init(Dim);

	int i,j,k;
	real v,u;
	support_vector	SV;
	SV.vec.newsize(Dim);
	SV.vec=SV.funvalue=0;

	for(i=0;i<Dim;i++)   // reserve for the basis vectors
	{
		m_lower->SVectors.push_back(SV); 
		m_upper->SVectors.push_back(SV);
	}

/* this part is to accelerate searching for repeated data, using binary search */
	vector<double> firstcomponent;
	vector<double> secondcomponent;

  if(test) {
	firstcomponent.reserve(K);
	if(dim>1) { // use binary search wrt the first two components, if dim>1
		secondcomponent.reserve(K);
		for(i=0;i<K;i++) {
			firstcomponent.push_back(x[i*dim]);
//			secondcomponent.push_back(x[i*dim+1]);
			// some trickery to merge 2 vector components into one double
			secondcomponent.push_back(merge(x[i*dim],x[i*dim+1]));
		}
		sort(firstcomponent.begin(),firstcomponent.end());
		sort(secondcomponent.begin(),secondcomponent.end());
	} else { // only the first component, if dim=1
		for(i=0;i<K;i++) {
			firstcomponent.push_back(x[i*dim]);
		}
		sort(firstcomponent.begin(),firstcomponent.end());
	}
  } // if test


// copy the data into SVectors, excluding repeated values
	k=0;
	for(i=0;i<K;i++) {
		u=0;
		for(j=0;j< Dim-1;j++) {
			v=SV.vec[j]=x[k]; k++;
			u+=v;
		}
		SV.vec[Dim-1]=1.0-u;

		SV.funvalue=y[i];

// to avoid repeated data points
		if(!test || 
		  !TestPresent(SV, firstcomponent, secondcomponent, m_lower->SVectors, dim)) {
			m_lower->SVectors.push_back(SV);
			m_upper->SVectors.push_back(SV);
		}
	}
    
// ensure Lipschitz constants were set correctly for dim=Dim

	if(m_lower->m_Constants.size()<Dim ) {  // if it was set but with incorrect dimension
		if(m_lower->m_Constants[0]>0)       // if it was actually set (0 means not set before)
			SetConstants(m_lower->m_Constants[0]);  // now we have right Dim
	}
}


void	STCInterpolant::SetDataColumn(int dim, int K, real* x, real* y, int test)
{
	m_lasterr=0;
	Dim=dim+1;

	m_lower->Init(Dim);
	m_upper->Init(Dim);

	int i,j,k;
	real v,u;
	support_vector	SV;
	SV.vec.newsize(Dim);
	SV.vec=SV.funvalue=0;

	for(i=0;i<Dim;i++)   // reserve for the basis vectors
	{
		m_lower->SVectors.push_back(SV); 
		m_upper->SVectors.push_back(SV);
	}

/* this part is to accelerate searching for repeated data, using binary search */
	vector<double> firstcomponent;
	vector<double> secondcomponent;

  if(test) {
	firstcomponent.reserve(K);
	if(dim>1) { // use binary search wrt the first two components, if dim>1
		secondcomponent.reserve(K);
		for(i=0;i<K;i++) {
			firstcomponent.push_back(x[i]);  // forst components of all data points live in the first column
			// some trickery to merge 2 vector components into one double
			secondcomponent.push_back(merge(x[i],x[i+K])); // first and second columns
		}
		sort(firstcomponent.begin(),firstcomponent.end());
		sort(secondcomponent.begin(),secondcomponent.end());
	} else { // only the first component, if dim=1
		for(i=0;i<K;i++) {
			firstcomponent.push_back(x[i]);  // first components of all data
		}
		sort(firstcomponent.begin(),firstcomponent.end());
	}
  } // if test


// copy the data into SVectors, excluding repeated values
	k=0;
	for(i=0;i<K;i++) {
		u=0;
		for(j=0;j< Dim-1;j++) {
			v=SV.vec[j]=x[j*K+i]; //k++;  // extract rows from columnwise representation
			u+=v;
		}
		SV.vec[Dim-1]=1.0-u;

		SV.funvalue=y[i];

// to avoid repeated data points
		if(!test || 
		  !TestPresent(SV, firstcomponent, secondcomponent, m_lower->SVectors, dim)) {
			m_lower->SVectors.push_back(SV);
			m_upper->SVectors.push_back(SV);
		}
	}
    
// ensure Lipschitz constants were set correctly for dim=Dim

	if(m_lower->m_Constants.size()<Dim ) {  // if it was set but with incorrect dimension
		if(m_lower->m_Constants[0]>0)       // if it was actually set (0 means not set before)
			SetConstants(m_lower->m_Constants[0]);  // now we have right Dim
	}
}



void	STCInterpolant::Construct()
{
	// assumed that getData was called and there is data in both m_lower and m_upper
	aux=(double*) malloc(sizeof(double)*Dim) ; // just to be sure it was initiated

	m_lower->Construct();
	m_upper->ConstructInv();

	if(m_lower->Match==1 || m_upper->Match==1)
		// LipConst is too small
		m_lasterr=ERR_LIP_LOW;
}

void	STCInterpolant::ConstructExplicit()
{
	// assumed that getData was called and there is data in both m_lower and m_upper
	aux=(double*) malloc(sizeof(double)*Dim) ; // just to be sure it was initiated

	m_lower->ConstructExplicit();
	m_upper->ConstructInvExplicit();
}


void	STCInterpolant::ComputeSlack(real* x)
{
	real u=0;
	for(int i=0;i<Dim-1;i++) { aux[i]=x[i]; u+=x[i];}
	aux[Dim-1]=1.-u;
}


real	STCInterpolant::Value(int dim, real* x)
{
	m_lasterr=0;
	if(dim>=Dim) // slack variable  computed
	{
		Lo=m_lower->FValueL(x);      // change from exhaustive to norma;
		Up= - m_upper->FValueL(x);
		return (Lo+Up)*0.5;
	} else { //compute slack
		ComputeSlack(x);
		Lo=m_lower->FValueL(aux);      // change from exhaustive to norma;
		Up= - m_upper->FValueL(aux);
		return (Lo+Up)*0.5;
	}
}
	
real	STCInterpolant::ValueExplicit(int dim, real* x)
{
	m_lasterr=0;
	if(dim>=Dim) // slack variable  computed
	{
		Lo=m_lower->FValueExplicit(x);      // change from exhaustive to norma;
		Up= - m_upper->FValueExplicit(x);
		return (Lo+Up)*0.5;
	} else { //compute slack
		ComputeSlack(x);
		Lo=m_lower->FValueExplicit(aux);      // change from exhaustive to norma;
		Up= - m_upper->FValueExplicit(aux);
		return (Lo+Up)*0.5;
	}
}


real	STCInterpolant::DetermineLipschitz()
{
// computes an estimate of the Lipschitz constant in simplicial distance from the data
// by computing distances and differences of function values for all pairs of data.
		Lip1=(double*) malloc(sizeof(double)*Dim) ; // just to be sure it was initiated
		Lip2=(double*) malloc(sizeof(double)*Dim) ; // just to be sure it was initiated

		support_vector *SU, *SV;

		int i,j,k,k1;
		real u,v;		
		
		for(i=0;i<Dim;i++) {
			Lip1[i]=0;
			Lip2[i]=0;
		}


		for(i=Dim;i<m_lower->SVectors.size();i++) {
			SU=&(m_lower->SVectors[i]);
			for(j=Dim;j<m_lower->SVectors.size();j++) 
				if(i != j) 
				{
					SV=&(m_lower->SVectors[j]);
					v=0;
					for(k=0;k<Dim;k++) {
						 u=SU->vec[k]-SV->vec[k]; // v is the distance in polyhedral norm
						 if(v<u) {v=u; k1=k;}  // and k1 is the direction
					}
					u=SU->funvalue - SV->funvalue; 
				  if(v > 0) // otherwise we are in trouble
				  { u= u/v;
					if(u>=0) {
						if(Lip1[k1] < u) 
							Lip1[k1]=u;
					} else { 
						if(Lip2[k1] < -u) 
							Lip2[k1]=-u;
					}
				  } // if v > 0
				} // for j
		} // for i

		u=0;
		for(i=0;i<Dim;i++) {
			Lip2[i]=Lip1[i]=max_(Lip1[i],Lip2[i])*1.000001;
			u=max_(u,Lip1[i]);
		}

		LipschitzConst=max_(u,1e-5);
		return LipschitzConst;
}


void	STCInterpolant::FreeMemory()
{
	m_upper->FreeMem();
	m_lower->FreeMem();
}


