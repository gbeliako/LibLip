/**************************************************************************

    begin                : April 30 2004
	version				 : 1.0 
    copyright            : (C) 2004 by Gleb Beliakov
    email                : gleb@deakin.edu.au


 *   This file contains declarations of two classes: Interpolant and       *
 *    STCInterpolant. STCInterpolant implements the method of multivariate *
 *   interpolation of Lipschitz functions using scattered data             *
 *                                                                         *
	Interpolant is a worker class which implements on-sided Lipschitz 
	interpolation of a function (from below), from scattered data.

	STCInterpolant is the class that implements Lipschitz interpolation
	which uses upper and lower interpolation. 

	The API interface is provided through STCInterpolant class 

    
STCInterpolant performs several functions:
	Receives the data for interpolation		SetData(dim,K,x,y)
	Constructs the interpolant				Construct() or ConstructExplicit()
	Computes the value of the interpolant   Value(x) or ValueExplicit(x)
	Computes the Lipschitz constant of the data set DetermineLipschitz()
	Sets the Lipschitz constant				SetConstants(LipConst)


STCInterpolant works as follows:
	After the data is received, it computes the slack variables for all data,
	and constructs the upper and lower interpolants. When Value(x) is needed,
	it evaluates the upper and lower interpolants, and takes the average.
	This value is the best approximation to the function f, which it interpolates,
	in the worst case scenario.

    There are 2 modes of evaluation: explicit (exhaustive comparison of
	K support function to compute their maximum or minimum), and "fast" method,
	which involves building a tree of local minima of the lower saw-tooth cover
	interpolant (maxima of the upper interpolant). This method takes the
	logarithmic time of the number of data points, but requires preprocessing,
	which is exponential in dim. Can be used for small number of variables <6,
	because otherwise the number of local minimizers is just too big, and
	explicit method becomes more efficient. The explicit method takes
	linear time of the number of data points.

Example of usage:

	STCInterpolant MyInt;
	MyInt.SetData(dim,K,x,y);
	where dim is dimension, K is the number of data points,
	x is a matrix containing data abscissae (in rows)
	y is a vector of function values

    MyInt.SetConstants(LipConst,dim);
    MyInt.Construct();    // or ConstructExplicit()

    double x[dim+1];
	x[0]=1; x[1]=3; ... x[dim]=1- sum x[i] (slack variable)

    r=MyInt.Value(dim+1,x);   // computes the value (or ValueExplicit(x))

  if necessary, Lipschitz constant can be computed from the data
    MyInt.DetermineLipschitz();
    
     

	See documentation about the methods used for further information

 *                                                                         *
 *  © Gleb Beliakov, 2004												   *
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



#if !defined(INTERPOL)
#define INTERPOL

#include "forest.h"

#define ERR_BOTH_FAIL 3
#define ERR_LO_FAIL 1
#define ERR_UP_FAIL 2
#define ERR_LIP_LOW 10
#define ERR_WRONG_LIP_LEN 11


/* Interpolant is a worker class, and is not part of API */
class Interpolant
{
public:
	SVDeque				SVectors;	// here support vectors live
	int					LastLabel;  // How many support vectors are there
	int					Dim;		// problem size +1 (slack variable)
	int					Match;		// flag indicating success of a query
	dVec				m_Constants;// here Lipschitz constants live

// SVSets live here
	Forest				HeapPossibleMin;  // needs to be public, accessed directly by SVSetnode

private:
// aux staff

	SVDeque::iterator	iter;
	SVDeque::iterator	t_iter; // just iterators

	real				C_minus; 

	dVec				temp_x;
	siVec				m_lastindex, m_TempSiVec;
	int					iteration, Iters;
	siVec				m_index, m_initvec; // the index of SVs 
	support_vector		m_sv;

	indexset			*m_indexset;  // the result of a query is returned here: the subset of s.functions
									  // to compute fValue
	indexsetiter		m_indexiter;  // iterator

public:
	Interpolant();
	virtual ~Interpolant();
	void	FreeMem();
	void	Init(int dim);


	// sets Lipschitz constant of the function
	void	SetConstants(dVec& newconst);
	// sets Lipschitz constant of the function (allows for different values for different coordinates
	void	SetConstants(real newconst, int dim); // create a vector, the last one *= sqrt(n-1)

// computes the value of the interpolant 
	virtual real	FValueL(real* x);		  // "fast" method for a small number of variables <6
	virtual real	FValueExplicit(real* x);  // just exhaustive search

//  construct the interpolant from the data already stored in SVectors 
	virtual	void	Construct();
//  inverts the sign of function values (for the upper interpolation)
	virtual	void	ConstructInv();

// as above, but does not construct the tree of local minima, will use explicit function evaluation
	virtual	void	ConstructExplicit();
//  inverts the sign of function values (for the upper interpolation)
	virtual	void	ConstructInvExplicit();

private:
// used by the fast method
	void	QueryDyn(real* x);
	void	ComputeCminus();

//  processes the data points and created the tree of local minima
	void	InitPopulateBoundaryPoints(); // points at infinity to kick start the algorithm
	void	LoadAdditionalPoints();

};


/*******************************************************************
	Class STCInterpolant

	Computes the value of the piecewise linear interpolant to the 
	multivariate scattered data using 2 methods:
	1) fast method requiring preprocessing
	2) slower direct method (no preprocessing)


********************************************************************/
#if !defined(STCINTERPOLANT)
#define STCINTERPOLANT

class STCInterpolant {

public:
	int				Dim;				// dimension +1 (slack variable)

private:
	Interpolant		*m_lower, *m_upper;	// the upper and lower interpolants

	real			LipschitzConst;		// LipschitzConstant of the function
	real			Lo,Up;				// lower and upper interpolant values

	int				m_lasterr;
	double			*aux, *Lip1, *Lip2;		// to compute Lipschitz constant of the data set
	double			*m_Constants;		// here Lipschitz constants live

public:
	STCInterpolant();
	~STCInterpolant(); // currently does nothing

// to set the Lipschitz constant of the function
	void	SetConstants(real newconst);			
	void	SetConstants();							// use computed lipschitz constants
	void	SetConstants(real newconst, int dim); // internal routine

// received the data set of dimension dim, of K data points
// test indicates the necessity to test whether all data are different (may be slow)
	void	SetData(int dim, int K, real* x, real* y, int test=0);

// the same as above, but uses fortran conventions for storing matrices (in columns)
	void	SetDataColumn(int dim, int K, real* x, real* y, int test=0);

// compute from the data set (KxKxdim operations)
	real	DetermineLipschitz();

// construct the interpolant (for small dimension <6)
	void	Construct();
// does not create the ree of local minima, just preproceses support vectors for
// subsequent explicit evaluation of the value
	void	ConstructExplicit();


// the member functions below perform the same computation, but use slightly different syntaxis

// computes the value of the interpolant, assuming that x already contains the slack variable
	real	Value(dVec& x);			// fast evaluation
	real	ValueExplicit(dVec& x);	// explicit evaluation, does not require preprocessing

// as above, but automatically computes the slack variable
	real	ValueSlack(dVec& x);			// fast evaluation
	real	ValueSlackExplicit(dVec& x);
	
// computes the value without using TNT library
	real	Value(int dim, real* x);
	real	ValueExplicit(int dim, real* x);

// computes the slack variable and stores it in Lip1
	void	ComputeSlack(dVec& x);
	void	ComputeSlack(real* x);


	int		LastError() {return m_lasterr;};

// this method is called when the interpolant is no longer needed, but is not automatically destroyed
// within the scope of its definition. Since the memory occupied can be fairly large, the user may wish
// to free the memory before the destructor does it. No other members can be called subsequently.
	void	FreeMemory();
};

#endif


#endif 
