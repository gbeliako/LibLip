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
#if !defined(SLIPINTERPOL)
#define SLIPINTERPOL

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;

#define DUAL

//#define LPSOLVE
#ifdef LPSOLVE
#include "lp_solve/lp_lib.h"
#else
extern "C" {
//#include <glpk.h>
#include "glpk/glpk.h"
}
#endif



// add procedural interface

#ifndef SLIPINTERPOL1


class SLipIntBasic {
public:
	double*   LipConst, MaxLipConst;  // an array of Lipschitz constants and the largest Lipschitz constant
// computed automatically in ComputeLipschitz

	int		*neighbors, *pneighbors;

	float	*GridR, *GridVal;
	int	*GridLim;

	double* Scaling;
	int	UseOtherBounds;

// diagnistics
	int		m_lasterror;
	int		m_number_constraints; 
	double  m_minvalue; 

	double OptimalPenalty;


	SLipIntBasic() {
		LipConst=0;
		MaxLipConst=0;
		Dim=0;
		NPTS=0;
		Scaling=0;
		m_lasterror=0;
		KeepCVProblem=0;
		UseOtherBounds=0;

		Indexsize=IndexsizeComp=ND=0;
		Index=NULL; IndexComp=NULL;
		LocalXData=LocalYData=LocalTData=YH=NULL; // pointers
		type=0;
		LocalCons=NULL;
		LocalRegion=NULL; LocalW=NULL;
		pneighbors=neighbors=NULL;

		GridR=GridVal=NULL;
		GridLim=NULL;
	}

	~SLipIntBasic() { 
		free(LipConst); 
		if(Scaling!=NULL) free(Scaling);
		if(pneighbors!=NULL) free(pneighbors);
		if(neighbors!=NULL) free(neighbors);
		if(GridR!=NULL) free(GridR);
		if(GridVal!=NULL) free(GridVal);
		if(GridLim!=NULL) free(GridLim);
	}

//************* MUST be implemented in the derived class ***********
// Computes the smallest Lipschitz constant, compatible with the data
	virtual	void	ComputeLipschitz(int dim, int npts, double* XData, double* YData)=0;

	virtual double dist(int dim,  double* x, double* xk, double* param=NULL) {return 0;};
	virtual double dist(int dim,  double* x, double* xk, int* Cons, double* param=NULL){return 0;}; // constrained
	virtual double distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL){return 0;};
	virtual double distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param=NULL){return 0;};
	virtual double distAll(int dim,  int type, double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL) {return 0;};
//************* MUST be implemented in the derived class ***********


// entry points
// type: 0 - usual interpolation, 1-constrained, 2- constrained Left , 3 - constrained right region
	void	ComputeLipschitzSplit(int dim, int npts, double* XData, double* YData, double* TData,  double ratio=0.5,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

	void	ComputeLipschitzCV(int dim, int npts, double* XData, double* YData, double* TData, 
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

// Computes the local Lipschitz constants in any norm, compatible with the data
	void	ComputeLocalLipschitz(int dim, int npts, double* XData, double* YData);
	void	ComputeLocalLipschitzCons(int dim, int npts, int _type, int* Cons, double* XData, double* YData, double* Region=NULL);



// Returns the value of the interpolant 
	double	Value(int dim, int npts,   double* x, double* XData, double* YData,   double LipConst, int* index=NULL);

// Returns the value of the interpolant, with the Lipschitz constant
// computed from the data. Can be used after ComputeLipschitz
	double	Value(int dim, int npts,   double* x, double* XData, double* YData, int* index=NULL);

// Returns the value of the interpolant, with the local Lipschitz constants
// computed from the data. Can be used after ComputeLocalLipschitz
	double	ValueLocal(int dim, int npts,   double* x, double* XData, double* YData);

	int	FindVoronoi(int dim, int npts,   double* x, double* XData, double &d);

	int	ComputeScaling(int dim, int npts, double* XData, double* YData);

// *** methods below refer to Monotone interpolation **

// returns 1 if x >> y wrt Cons
	int		Dominates(int dim, double* x, double * y, int* Cons);

// Returns the value of the monotone interpolant 
	double	ValueCons(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   double LipConst, int* index=NULL);

// Returns the value of the monotone interpolant, with the Lipschitz constant
// computed from the data. Can be used after ComputeLipschitz
	double	ValueCons(int dim, int npts, int* Cons,  double* x, double* XData, double* YData);

// Returns the value of the interpolant , assuming it is monotone for x<< LeftRegion
	double	ValueConsLeftRegion(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* LeftRegion, int* index=NULL);
// Returns the value of the interpolant in l_2 norm,  assuming it is monotone for x>> RightRegion
	double	ValueConsRightRegion(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* RightRegion, int* index=NULL);


// Returns the value of the monotone interpolant, with the local Lipschitz constants
// computed from the data. Can be used after ComputeLocalLipschitz
	double	ValueLocalCons(int dim, int npts, int* Cons, double* x, double* XData, double* YData);
// same but monotonicity is for x<< LeftRegion or x>>RightRegion
	double	ValueLocalConsLeftRegion(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, double* Region);
	double	ValueLocalConsRightRegion(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, double* Region);

// Verifies the data is monotone wrt specified variables
	int		VerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC=10e20, double eps=1e-7);

// Verifies the data is monotone wrt specified variables in the region On x<< LeftBoundary 
	int		VerifyMonotonicityLeftRegion (int dim, int npts, int* Cons,  double* XData, double* YData, 
															 double* LeftRegion, double LC=10e20,double eps=1e-7);
// Verifies the data is monotone wrt specified variables in the region On x >> Rightboundary 
	int		VerifyMonotonicityRightRegion (int dim, int npts, int* Cons,  double* XData, double* YData, 
															  double* RightRegion, double LC=10e20,double eps=1e-7);

// the working horses
	double	ValueLocal2Consinternal(int dim, int npts, int* Cons,  double* x, double* XData, double* YData, int reg, double* Region);

	void	SmoothLipschitz2internal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,int Cf,
											double* LC, double* W, int* Cons, int region=0, double* Region=NULL, int* index=NULL);
	void	SmoothLipschitz2internalUpdate(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,int Cf,
											double* LC, double* W, int* Cons, int region=0, double* Region=NULL, int* index=NULL);



// Assumes the data in XData are stored columnwise, as in Fortran.
// Used for compatibility of this library with other packages, e.g.,Matlab, R, which
// may use column format. 
	void	ConvertXData(int dim, int npts,  double* XData);
	void	ConvertXData(int dim, int npts,  double* XData, double* auxStorage);


// implement the methods for sample splitting and CV. used internally
	virtual int		ComputeSmoothenedSplit();
	virtual int		ComputeLipschitzFinal();
	virtual int		ComputeFitLipschitzCV(int excluded);

	virtual double  ExtraUpperBound(int dim, double* x, double * param) {return 10e20;};  // derived classes may overwrite these
	virtual double  ExtraLowerBound(int dim, double* x, double * param) {return -10e20;};

	virtual double	Fun(double x); // for golden section algorithm
	double MinFuncSplit(double x);
	double MinFuncCV(double x);
	double MinFuncLocalSplit(double x);

	virtual double value(int dim, int npts,   double* x, double* XData, double* YData,   double LipConst, int* index=NULL,
		 int type=0, int* Cons=NULL, double* Region=NULL	);
// various parameters
// interpretation depends on the derived class. In this class:
// type: 0 - usual interpolation, 1-constrained, 2- constrained Left , 3 - constrained right region
	double valuelocal(int dim, int npts,   double* x, double* XData, double* YData, int type, int* Cons, double* Region);


// called internally
	double	ComputeFitIndexCV();
	double	ComputeFitIndex();
	void	PrepareLipschitzSplit(double SplitP);
	void	PrepareLipschitzCV();
	double	golden(double A, double B); 

	int BinSearch(double r, float* Arr, int le, int ri);

//these are private, but need to be inherited, so declared as public
	double	M; // temp value of the LipConst
	double g1,g2,d1,d2,d3;
	int i,j,i1;
	int Dim,NPTS; 
	int	TotalNeighbors;

// these vars are for the cross-validation/ sample splitting
	int		Indexsize,IndexsizeComp,ND;
	int		*Index, *IndexComp;
	double	*LocalXData, *LocalYData, *LocalTData; // pointers
	double  *YH;
// parameters to be passed to value() method and other CV and splitting routines
	int		type;
	int		*LocalCons;
	double	*LocalRegion, *LocalW;
	double  *AuxXData; 

	int		TypeLipEstimate; //0 sample splitting, 1 CV
	int		KeepCVProblem;

//	double Gamma;

#ifdef LPSOLVE
	lprec		*MyLP; 
#else
	LPX			*MyLP; 
#endif

};


class SLipInt:public SLipIntBasic {
public:

// In the methods below, XData contain the abscissae of data points x^k (arranged
// in rows (C-convention)) and YData contain y^k. x is the point at which g(x) is needed.

// Computes the smallest Lipschitz constant in l_2 norm, compatible with the data
	virtual void	ComputeLipschitz(int dim, int npts, double* XData, double* YData);



// Methods below refer to Lipschitz smoothing **************************

// Smooth the data subject to given Lipschitz constant in Euclidean norm
	void	SmoothLipschitz(int dim, int npts,  double* XData, double* YData, double* TData, double LC);
// Smooth the data subject to given Lipschitz constant in Euclidean norm, subject to weightings
	void	SmoothLipschitzW(int dim, int npts,  double* XData, double* YData, double* TData, double LC, double* W);


// same, subject to monotonicity constraints
	void	SmoothLipschitzCons(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC);
	void	SmoothLipschitzWCons(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double LC, double* W);
	void	SmoothLipschitzConsLeftRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* LeftRegion);
	void	SmoothLipschitzConsRightRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* RightRegion);
	void	SmoothLipschitzWConsLeftRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* W, double* LeftRegion);
	void	SmoothLipschitzWConsRightRegion(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC, double* W, double* RightRegion);



	virtual double dist(int dim,  double* x, double* xk, double* param=NULL);
	virtual double dist(int dim,  double* x, double* xk, int* Cons, double* param=NULL); // constrained
	virtual double distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);
	virtual double distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param=NULL);
	virtual double distInf1(int dim,  double* x, double* xk, int* dir);
	virtual double distAll(int dim,  int type, double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);

private:


// Computes the smallest Lipschitz constant in l_inf norm, compatible with the data, used internally for scaling
	void	ComputeLipschitzInf(int dim, int npts, double* XData, double* YData);
};



// Same as SLipInt, but uses l_infty rather than Euclidean norm


class SLipIntInf:public SLipIntBasic {
public:

// In the methods below, XData contain the abscissae of data points x^k (arranged
// in rows (C-convention)) and YData contain y^k. x is the point at which g(x) is needed.

// Computes the smallest Lipschitz constant in l_inf norm, compatible with the data
	virtual void	ComputeLipschitz(int dim, int npts, double* XData, double* YData);


// Same, but uses an array of Lipschitz constants wrt each variable
	double	ValueDir(int dim, int npts, double* x, double* XData, double* YData,   double* LipConst, int* index=NULL);

// Returns the value of the interpolant, with the Lipschitz constant
// computed from the data. Can be used after ComputeLipschitz
	double	ValueDir(int dim, int npts, double* x, double* XData, double* YData);


// Returns the value of the interpolant , with the Lipschitz constant
// computed from the data. Can be used after ComputeLipschitz
	double	ValueConsDir(int dim, int npts, int* Cons, double* x, double* XData, double* YData);



// Methods below refer to Lipschitz smoothing

// Smooth the data subject to given Lipschitz constant in l-infy norm
	void	SmoothLipschitz(int dim, int npts,  double* XData, double* YData, double* TData, double LC);
// Smooth the data subject to given Lipschitz constant in l-infy norm, subject to weightings
	void	SmoothLipschitzW(int dim, int npts,  double* XData, double* YData, double* TData, double LC, double* W);

// same, subject to monotonicity constraints
	void	SmoothLipschitzCons(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double LC);
	void	SmoothLipschitzWCons(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double LC, double* W);


// Same in simplicial distance
	void	SmoothLipschitzSimp(int dim, int npts,  double* XData, double* YData, double* TData, double LC);
// Same in simplicial distance, subject to weightings
	void	SmoothLipschitzSimpW(int dim, int npts,  double* XData, double* YData, double* TData, double LC, double* W);

	virtual double dist(int dim,  double* x, double* xk, double* param=NULL);
	virtual double dist(int dim,  double* x, double* xk, int* Cons, double* param=NULL); // constrained
	virtual double distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);
	virtual double distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param=NULL);

	virtual double distDir(int dim, double* x, double* xk, int* dir); // if we need the direction (just the coord)
	virtual double distInfDir(int dim,  double* x, double* xk, int* dir); // if we need the direction + left/right

	virtual double distSimp(int dim,  double* x, double* xk, int* dir);
	virtual double distAll(int dim,  int type, double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);


	virtual int		ComputeSmoothenedSplit();
	virtual int		ComputeLipschitzFinal();
	virtual int		ComputeFitLipschitzCV(int excluded);

private:
// the working horses
	void	SmoothLipschitzInfinternal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,
											double* LC, double* W, int* index=NULL);
	void	SmoothLipschitzSimpinternal(int dim, int npts,  double* XData, double* YData, double* TData,  int Wf,
											 double LC, double* W, int* index=NULL);
};



class SLipIntLp:public SLipInt {
public:
	double m_P, m_P1;

	SLipIntLp(): SLipInt() { m_P=m_P1=1;}

	void      SetP(double p) {m_P=p; if(m_P<= 1.001) m_P=1.0; m_P1=1.0/m_P; };
	double    GetP() {return m_P;};

	virtual double dist(int dim,  double* x, double* xk, double* param=NULL);
	virtual double dist(int dim,  double* x, double* xk, int* Cons, double* param=NULL); // constrained
	virtual double distLeftRegion(int dim,  double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);
	virtual double distRightRegion(int dim,  double* x, double* xk, int* Cons, double* RightRegion, double* param=NULL);
	virtual double distAll(int dim,  int type, double* x, double* xk, int* Cons, double* LeftRegion, double* param=NULL);

};


// not documented
// a classifier based on Lipschitz interpolation
class SLipClass: public SLipInt {

public:
	double	Penalty, SmoothingParam;

	SLipClass(): SLipInt() {
		SmoothingParam=0;
		Penalty=1;
	}
	
	// Returns the value of the Classifier 
	int		ValueClass(int dim, int npts, double* x, double* XData, double* YData,   double LipConst);
	int		ValueConsClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   double LipConst);
	int		ValueConsLeftRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* LeftRegion);
	int		ValueConsRightRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double LipConst, double* RightRegion);

	int		ValueLocalClass(int dim, int npts, double* x, double* XData, double* YData);
	int		ValueLocalConsClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData);
	int		ValueLocalConsLeftRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													   double* LeftRegion);
	int		ValueLocalConsRightRegionClass(int dim, int npts,  int* Cons, double* x, double* XData, double* YData,   
													  double* RightRegion);

// performs large margin classifier smoothing of the data
	void	SmoothLipschitzClass(int dim, int npts,  double* XData, double* YData, double* TData, double *LC);
	void	SmoothLipschitzWClass(int dim, int npts,  double* XData, double* YData, double* TData, double *LC, double* W);
	void	SmoothLipschitzConsClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC);
	void	SmoothLipschitzWConsClass(int dim, int npts,  int* Cons, double* XData, double* YData, double* TData, double *LC, double* W);
	void	SmoothLipschitzConsLeftRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* LeftRegion);
	void	SmoothLipschitzConsRightRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* RightRegion);
	void	SmoothLipschitzWConsLeftRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* W, double* LeftRegion);
	void	SmoothLipschitzWConsRightRegionClass(int dim, int npts, int* Cons, double* XData, double* YData, double* TData,  double *LC, double* W, double* RightRegion);

// the working horse
	void	SmoothLipschitz2Classinternal(int dim, int npts,  double* XData, double* YData, double* TData, int LCf, int Wf,int Cf,
											double* LC, double* W, int* Cons, int region=0, double* Region=NULL, int* index=NULL);


};






#endif


#endif 


/*
#if defined(__cplusplus)
extern "C"
{
#endif
#include <qhull.h>
#include <mem.h>
#include <qset.h>
#include <geom.h>
#include <merge.h>
#include <poly.h>
#include <io.h>
#include <stat.h>
#if defined(__cplusplus)
}
#endif
*/