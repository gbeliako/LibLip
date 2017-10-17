/**************************************************************************

 ***************************************************************************/

#include "liblip.h"

#include "liblipc.h"

// global variables

// an instance of the interpolant class
STCInterpolant		gl;

// simple Lipschitz interpolant
SLipInt					sli;
// simple Lipschitz interpolant
SLipIntInf				slii;

// Lipschitz constant (not yet set)
double	GlobalLip=0;



void	STCSetLipschitz(real x) {GlobalLip=x;}

int	STCBuildLipInterpolant(int Dim, int Ndata, double* x, double* y)
{
	gl.SetData(Dim,Ndata,x,y);

	// Lipschitz constants live here
	if(GlobalLip<=0) {
		gl.DetermineLipschitz();
		gl.SetConstants();			// automatic
	} else
		gl.SetConstants(GlobalLip,Dim+1);  // if it was specified

	gl.Construct();

	return gl.LastError();
//	if(gl.LastError()==ERR_LIP_LOW) cout << "Lipschitz const low or data coincide" << endl;
}

int	STCBuildLipInterpolantExplicit(int Dim, int Ndata,  double* x, double* y)
{
	gl.SetData(Dim,Ndata,x,y);

	// Lipschitz constants live here
	if(GlobalLip<=0) {
		gl.DetermineLipschitz();
		gl.SetConstants();			// automatic, but slow
	} else
		gl.SetConstants(GlobalLip,Dim+1);

	gl.ConstructExplicit();

	return gl.LastError();

//	if(gl.LastError()==ERR_LIP_LOW) cout << "Lipschitz const low or data coincide" << endl;
}

// the methods below are identical to the above, but use columnwise storage of matrices
int	STCBuildLipInterpolantColumn(int Dim, int Ndata, double* x, double* y)
{
	gl.SetDataColumn(Dim,Ndata,x,y);

	// Lipschitz constants live here
	if(GlobalLip<=0) {
		gl.DetermineLipschitz();
		gl.SetConstants();			// automatic
	} else
		gl.SetConstants(GlobalLip,Dim+1);  // if it was specified

	gl.Construct();

	return gl.LastError();
//	if(gl.LastError()==ERR_LIP_LOW) cout << "Lipschitz const low or data coincide" << endl;
}

int	STCBuildLipInterpolantExplicitColumn(int Dim, int Ndata,  double* x, double* y)
{
	gl.SetDataColumn(Dim,Ndata,x,y);

	// Lipschitz constants live here
	if(GlobalLip<=0) {
		gl.DetermineLipschitz();
		gl.SetConstants();			// automatic, but slow
	} else
		gl.SetConstants(GlobalLip,Dim+1);

	gl.ConstructExplicit();

	return gl.LastError();

//	if(gl.LastError()==ERR_LIP_LOW) cout << "Lipschitz const low or data coincide" << endl;
}


 double	STCValue( double* x )
{
	return gl.Value(gl.Dim-1,x); // need to compute the slack variable
}

 double	STCValueExplicit( double* x )
{
	return gl.ValueExplicit(gl.Dim-1,x);
}


 void	STCFreeMemory() {gl.FreeMemory();}




/*--------------------------------------------------------*/
/* interface to the members of SLipInt class */
double	LipIntValue(int Dim, int Ndata, double* x, double* Xd,double* y,  double Lipconst, int* Index)
{ return sli.Value(Dim, Ndata, x, Xd, y, Lipconst, Index); }

double	LipIntValueAuto(int Dim, int Ndata, double* x,double* Xd, double* y, int* Index)
{ return sli.Value(Dim, Ndata, x, Xd,y, Index); }

double	LipIntValueCons(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, int* Index)
{ return sli.ValueCons(Dim, Ndata, Cons, x, Xd, y, Lipconst, Index); }

double	LipIntValueConsLeftRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, double* Region, int* Index)
{ return sli.ValueConsLeftRegion(Dim, Ndata, Cons, x, Xd, y, Lipconst, Region, Index); }

double	LipIntValueConsRightRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, double* Region, int* Index)
{ return sli.ValueConsRightRegion(Dim, Ndata, Cons, x, Xd, y, Lipconst, Region, Index); }


double	LipIntValueLocal(int Dim, int Ndata, double* x, double* Xd,double* y)
{ return sli.ValueLocal(Dim, Ndata, x, Xd,y); }

double	LipIntValueLocalCons(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y)
{ return sli.ValueLocalCons(Dim, Ndata, Cons, x, Xd,y); }

double	LipIntValueLocalConsLeftRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region)
{ return sli.ValueLocalConsLeftRegion(Dim, Ndata, Cons, x, Xd,y,Region); }

double	LipIntValueLocalConsRightRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region)
{ return sli.ValueLocalConsRightRegion(Dim, Ndata, Cons, x, Xd,y,Region); }


void	LipIntComputeLipschitz(int Dim, int Ndata, double* x, double* y)
{  sli.ComputeLipschitz(Dim, Ndata, x, y); }

void	LipIntComputeLocalLipschitz(int Dim, int Ndata, double* x, double* y)
{ sli.ComputeLocalLipschitz(Dim, Ndata, x, y);}

void	LipIntComputeLipschitzCV(int Dim, int Ndata, double* Xd, double* y, double* T,
			int type, int* Cons, double* Region, double *W)
{	sli.ComputeLipschitzCV(Dim,  Ndata, Xd,  y,  T, type,  Cons,  Region,  W); }

void	LipIntComputeLipschitzSplit(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio,
			int type, int* Cons, double* Region, double *W)
{	sli.ComputeLipschitzSplit(Dim,  Ndata, Xd,  y,  T, ratio, type,  Cons,  Region,  W); }

void	LipIntComputeLocalLipschitzSmooth(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio,
			int type, int* Cons, double* Region, double *W)
{	sli.ComputeLocalLipschitzSmooth(Dim,  Ndata, Xd,  y,  T, ratio, type,  Cons,  Region,  W); }

void	LipIntComputeLocalLipschitzPenalty(int Dim, int Ndata, double* Xd, double* y, double* T, double penalty,
			int type, int* Cons, double* Region, double *W)
{	sli.ComputeLocalLipschitzPenalty(Dim,  Ndata, Xd,  y,  T, penalty, type,  Cons,  Region,  W); }

void	LipIntSmoothLipschitz(int Dim, int Ndata,  double* Xd, double* y, double* T,  double LC, 
							  int fW, int fC, int fR, double* W, int* Cons, double* Region)
{ // fR is 0, 1-left, 2-right
	sli.SmoothLipschitz2internal(Dim,Ndata,Xd,  y,  T, 0,fW, fC, &LC,  W, Cons, fR, Region);
}

void	LipIntSetGamma(double g)
{	sli.Gamma=g;  }

double	LipIntGetLipConst() 
{ return sli.MaxLipConst; }

void	LipIntGetScaling(double *S) 
{	int i;
	for(i=0;i<sli.NPTS;i++) 
	S[i]=sli.Scaling[i]; 
}


int		LipIntComputeScaling(int Dim, int Ndata, double* XData, double* YData)
{	return sli.ComputeScaling(Dim, Ndata, XData,YData); }



void	ConvertXData(int dim, int npts,  double* XData)
{    sli.ConvertXData(dim, npts, XData); }

void	ConvertXData(int dim, int npts,  double* XData, double *auxdata)
{    sli.ConvertXData(dim, npts, XData,auxdata); }


int		LipIntVerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC, double eps)
{	return sli.VerifyMonotonicity(dim,npts,Cons,XData,YData,LC,eps); }

int		LipIntVerifyMonotonicityLeftRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, double LC, double eps)
{	return sli.VerifyMonotonicityLeftRegion(dim,npts,Cons,XData,YData,Region,LC,eps); }

int		LipIntVerifyMonotonicityRightRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, double LC, double eps)
{	return sli.VerifyMonotonicityRightRegion(dim,npts,Cons,XData,YData,Region,LC,eps); }




/* interface to the members of SLipIntInf class ====================================== */
double	LipIntInfValue(int Dim, int Ndata, double* x, double* Xd,double* y,  double Lipconst, int* Index)
{ return slii.Value(Dim, Ndata, x, Xd, y, Lipconst, Index); }

double	LipIntInfValueAuto(int Dim, int Ndata, double* x,double* Xd, double* y, int* Index)
{ return slii.Value(Dim, Ndata, x, Xd,y, Index); }

double	LipIntInfValueCons(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, int* Index)
{ return slii.ValueCons(Dim, Ndata, Cons, x, Xd, y, Lipconst, Index); }

double	LipIntInfValueConsLeftRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, double* Region, int* Index)
{ return sli.ValueConsLeftRegion(Dim, Ndata, Cons, x, Xd, y, Lipconst, Region, Index); }

double	LipIntInfValueConsRightRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, double* Region, int* Index)
{ return slii.ValueConsRightRegion(Dim, Ndata, Cons, x, Xd, y, Lipconst, Region, Index); }


double	LipIntInfValueLocal(int Dim, int Ndata, double* x, double* Xd,double* y)
{ return slii.ValueLocal(Dim, Ndata, x, Xd,y); }

double	LipIntInfValueLocalCons(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y)
{ return slii.ValueLocalCons(Dim, Ndata, Cons, x, Xd,y); }

double	LipIntInfValueLocalConsLeftRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region)
{ return slii.ValueLocalConsLeftRegion(Dim, Ndata, Cons, x, Xd,y,Region); }

double	LipIntInfValueLocalConsRightRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region)
{ return slii.ValueLocalConsRightRegion(Dim, Ndata, Cons, x, Xd,y,Region); }


void	LipIntInfComputeLipschitz(int Dim, int Ndata, double* x, double* y)
{  slii.ComputeLipschitz(Dim, Ndata, x, y); }

void	LipIntInfComputeLocalLipschitz(int Dim, int Ndata, double* x, double* y)
{ slii.ComputeLocalLipschitz(Dim, Ndata, x, y);}

void	LipIntInfComputeLipschitzCV(int Dim, int Ndata, double* Xd, double* y, double* T,
			int type, int* Cons, double* Region, double *W)
{	slii.ComputeLipschitzCV(Dim,  Ndata, Xd,  y,  T, type,  Cons,  Region,  W); }

void	LipIntInfComputeLipschitzSplit(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio,
			int type, int* Cons, double* Region, double *W)
{	slii.ComputeLipschitzSplit(Dim,  Ndata, Xd,  y,  T, ratio, type,  Cons,  Region,  W); }

void	LipIntInfComputeLocalLipschitzSmooth(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio,
			int type, int* Cons, double* Region, double *W)
{	slii.ComputeLocalLipschitzSmooth(Dim,  Ndata, Xd,  y,  T, ratio, type,  Cons,  Region,  W); }

void	LipIntInfComputeLocalLipschitzPenalty(int Dim, int Ndata, double* Xd, double* y, double* T, double penalty,
			int type, int* Cons, double* Region, double *W)
{	slii.ComputeLocalLipschitzPenalty(Dim,  Ndata, Xd,  y,  T, penalty, type,  Cons,  Region,  W); }

void	LipIntInfSmoothLipschitz(int Dim, int Ndata,  double* Xd, double* y, double* T,  double LC, 
							  int fW, int fC, int fR, double* W, int* Cons, double* Region)
{ // fR is 0, 1-left, 2-right
	slii.SmoothLipschitz2internal(Dim,Ndata,Xd,  y,  T, 0,fW, fC, &LC,  W, Cons, fR, Region);
}

void	LipIntInfSetGamma(double g)
{	slii.Gamma=g;  }

double	LipIntInfGetLipConst() 
{ return slii.MaxLipConst; }

void	LipIntInfGetScaling(double *S) 
{	int i;
	for(i=0;i<sli.NPTS;i++) 
	S[i]=slii.Scaling[i]; 
}


int		LipIntInfComputeScaling(int Dim, int Ndata, double* XData, double* YData)
{	return slii.ComputeScaling(Dim, Ndata, XData,YData); }




int		LipIntInfVerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC, double eps)
{	return slii.VerifyMonotonicity(dim,npts,Cons,XData,YData,LC,eps); }

int		LipIntInfVerifyMonotonicityLeftRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, double LC, double eps)
{	return slii.VerifyMonotonicityLeftRegion(dim,npts,Cons,XData,YData,Region,LC,eps); }

int		LipIntInfVerifyMonotonicityRightRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, double LC, double eps)
{	return slii.VerifyMonotonicityRightRegion(dim,npts,Cons,XData,YData,Region,LC,eps); }







void	LipIntInfSmoothLipschitzSimp(int dim, int npts,  double* XData, double* YData, double* TData,  double LC)
{	slii.SmoothLipschitzSimp(dim,npts,XData,YData,TData,LC);}

void	LipIntInfSmoothLipschitzSimpW(int dim, int npts,  double* XData, double* YData, double* TData,  double LC, double* W)
{	slii.SmoothLipschitzSimpW(dim,npts,XData,YData,TData,LC,W);}




