/**************************************************************************

 ***************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif


#define NULL 0

/* interface to the members of SLipInt class ===================== */
double	LipIntValue(int Dim, int Ndata, double* x, double* Xd,double* y,  double Lipconst, int* Index=NULL);

double	LipIntValueAuto(int Dim, int Ndata, double* x,double* Xd, double* y, int* Index=NULL);

double	LipIntValueCons(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, int* Index=NULL);

double	LipIntValueConsLeftRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, 
								  double* Region, int* Index=NULL);

double	LipIntValueConsRightRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, 
								   double* Region, int* Index=NULL);


double	LipIntValueLocal(int Dim, int Ndata, double* x, double* Xd,double* y);

double	LipIntValueLocalCons(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y);

double	LipIntValueLocalConsLeftRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);

double	LipIntValueLocalConsRightRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);


void	LipIntComputeLipschitz(int Dim, int Ndata, double* x, double* y);

void	LipIntComputeLocalLipschitz(int Dim, int Ndata, double* x, double* y);

void	LipIntComputeLipschitzCV(int Dim, int Ndata, double* Xd, double* y, double* T,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntComputeLipschitzSplit(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio=0.5,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntComputeLocalLipschitzSmooth(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio=0.5,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntComputeLocalLipschitzPenalty(int Dim, int Ndata, double* Xd, double* y, double* T, double penalty,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntSmoothLipschitz(int Dim, int Ndata,  double* Xd, double* y, double* T,  double LC, 
							  int fW=0, int fC=0, int fR=0, double* W=NULL, int* Cons=NULL, double* Region=NULL);
 // fR is 0, 1-left, 2-right

void	LipIntSetGamma(double g);
double	LipIntGetLipConst() ;
void	LipIntGetScaling(double *S) ;
int		LipIntComputeScaling(int Dim, int Ndata, double* XData, double* YData);
void	ConvertXData(int dim, int npts,  double* XData);
void	ConvertXDataAux(int dim, int npts,  double* XData, double *auxdata);

int		LipIntVerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC=10e20, double eps=1e-7);
int		LipIntVerifyMonotonicityLeftRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, 
										   double LC=10e20, double eps=1e-7);
int		LipIntVerifyMonotonicityRightRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, 
											double LC=10e20, double eps=1e-7);


/* interface to the members of SLipIntInf class ====================================== */

double	LipIntInfValue(int Dim, int Ndata, double* x, double* Xd,double* y,  double Lipconst, int* Index=NULL);

double	LipIntInfValueAuto(int Dim, int Ndata, double* x,double* Xd, double* y, int* Index=NULL);

double	LipIntInfValueCons(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, int* Index=NULL);

double	LipIntInfValueConsLeftRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, 
								  double* Region, int* Index=NULL);

double	LipIntInfValueConsRightRegion(int Dim, int Ndata, int* Cons, double* x, double* Xd,double* y,  double Lipconst, 
								   double* Region, int* Index=NULL);


double	LipIntInfValueLocal(int Dim, int Ndata, double* x, double* Xd,double* y);

double	LipIntInfValueLocalCons(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y);

double	LipIntInfValueLocalConsLeftRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);

double	LipIntInfValueLocalConsRightRegion(int Dim, int Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);


void	LipIntInfComputeLipschitz(int Dim, int Ndata, double* x, double* y);

void	LipIntInfComputeLocalLipschitz(int Dim, int Ndata, double* x, double* y);

void	LipIntInfComputeLipschitzCV(int Dim, int Ndata, double* Xd, double* y, double* T,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntInfComputeLipschitzSplit(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio=0.5,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntInfComputeLocalLipschitzSmooth(int Dim, int Ndata, double* Xd, double* y, double* T, double ratio=0.5,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntInfComputeLocalLipschitzPenalty(int Dim, int Ndata, double* Xd, double* y, double* T, double penalty,
			int type=0, int* Cons=NULL, double* Region=NULL, double *W=NULL);

void	LipIntInfSmoothLipschitz(int Dim, int Ndata,  double* Xd, double* y, double* T,  double LC, 
							  int fW=0, int fC=0, int fR=0, double* W=NULL, int* Cons=NULL, double* Region=NULL);
 // fR is 0, 1-left, 2-right

void	LipIntInfSetGamma(double g);
double	LipIntInfGetLipConst() ;
void	LipIntInfGetScaling(double *S) ;
int		LipIntInfComputeScaling(int Dim, int Ndata, double* XData, double* YData);

int		LipIntInfVerifyMonotonicity(int dim, int npts, int* Cons,  double* XData, double* YData, double LC=10e20, double eps=1e-7);
int		LipIntInfVerifyMonotonicityLeftRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, 
										   double LC=10e20, double eps=1e-7);
int		LipIntInfVerifyMonotonicityRightRegion(int dim, int npts, int* Cons,  double* XData, double* YData, double* Region, 
											double LC=10e20, double eps=1e-7);


void	LipIntInfSmoothLipschitzSimp(int dim, int npts,  double* XData, double* YData, double* TData,  double LC);
void	LipIntInfSmoothLipschitzSimpW(int dim, int npts,  double* XData, double* YData, double* TData,  double LC, double* W);




/* interface to the members of STCInterpolant class ====================================== */

// supplies the data to Interpolant and constructs the interpolant
// assuming a given Lipschitz constant, supplied by SetLipschitz
// if LipConstant was not supplied, tries to find it from the data
// assumes that all data are different. 

 int	STCBuildLipInterpolant(int Dim, int Ndata, double* x, double* y);

// as above, but for explicit evaluation, needs no preprocessing, but may be slower
 int	STCBuildLipInterpolantExplicit(int Dim, int Ndata, double* x, double* y);

// in the methods above, the coordinates of the data points in x are stored in rows

// the following methods store data in columns (like in fortran or Matlab)
// they use the transposed of the matrix x 
 int	STCBuildLipInterpolantColumn(int Dim, int Ndata, double* x, double* y);

// as above, but for explicit evaluation, needs no preprocessing, but may be slower
 int	STCBuildLipInterpolantExplicitColumn(int Dim, int Ndata, double* x, double* y);


// specify the Lipschitz constant for your function
 void	STCSetLipschitz(double x);

// computes the value of the interpolant at any given point x
 double	STCValue( double* x );

// same but using explicit evaluation with no preprocessing
 double	STCValueExplicit( double* x );

 void	STCFreeMemory();






#ifdef __cplusplus
}
#endif
