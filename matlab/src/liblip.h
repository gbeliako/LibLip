/**************************************************************************

  Procedural intervace to the methods of Lipschitz interpolant classes
 ***************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif

typedef double ( *USER_FUNCTION)(int *, double *, double *);
void LipIntSetBounds(USER_FUNCTION low, USER_FUNCTION up);
void LipIntClearBounds();
//#define NULL 0

/* interface to the members of SLipInt class ===================== */
double	LipIntValue(int* Dim, int* Ndata, double* x, double* Xd,double* y,  double* Lipconst, int* Index);

double	LipIntValueAuto(int* Dim, int* Ndata, double* x,double* Xd, double* y, int* Index);

double	LipIntValueCons(int* Dim, int* Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, int* Index);

double	LipIntValueConsLeftRegion(int* Dim, int* Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, 
								  double* Region, int* Index);

double	LipIntValueConsRightRegion(int* Dim, int* Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, 
								   double* Region, int* Index);


double	LipIntValueLocal(int *Dim, int *Ndata, double* x, double* Xd,double* y);

double	LipIntValueLocalCons(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y);

double	LipIntValueLocalConsLeftRegion(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);

double	LipIntValueLocalConsRightRegion(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);


void	LipIntComputeLipschitz(int *Dim, int *Ndata, double* x, double* y);

void	LipIntComputeLocalLipschitz(int *Dim, int *Ndata, double* x, double* y);

void	LipIntComputeLipschitzCV(int *Dim, int *Ndata, double* Xd, double* y, double* T,
			int* type, int* Cons, double* Region, double *W);

void	LipIntComputeLipschitzSplit(int *Dim, int *Ndata, double* Xd, double* y, double* T, double* ratio,
			int* type, int* Cons, double* Region, double *W);


void	LipIntSmoothLipschitz(int *Dim, int *Ndata,  double* Xd, double* y, double* T,  double* LC, 
							  int* fW, int* fC, int* fR, double* W, int* Cons, double* Region);
 // fR is 0, 1-left, 2-right

void	LipIntSetGamma(double* g);
double	LipIntGetLipConst() ;
void	LipIntGetScaling(double *S) ;
int		LipIntComputeScaling(int *Dim, int *Ndata, double* XData, double* YData);
void	ConvertXData(int *Dim, int* npts,  double* XData);
void	ConvertXDataAux(int *Dim, int* npts,  double* XData, double *auxdata);

int		LipIntVerifyMonotonicity(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* LC, double* eps);
int		LipIntVerifyMonotonicityLeftRegion(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* Region, 
										   double* LC, double* eps);
int		LipIntVerifyMonotonicityRightRegion(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* Region, 
											double* LC, double* eps);


/* interface to the members of SLipIntInf class ====================================== */

double	LipIntInfValue(int *Dim, int *Ndata, double* x, double* Xd,double* y,  double* Lipconst, int* Index);

double	LipIntInfValueAuto(int *Dim, int *Ndata, double* x,double* Xd, double* y, int* Index);

double	LipIntInfValueCons(int *Dim, int *Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, int* Index);

double	LipIntInfValueConsLeftRegion(int *Dim, int *Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, 
								  double* Region, int* Index);

double	LipIntInfValueConsRightRegion(int *Dim, int *Ndata, int* Cons, double* x, double* Xd,double* y,  double* Lipconst, 
								   double* Region, int* Index);


double	LipIntInfValueLocal(int *Dim, int *Ndata, double* x, double* Xd,double* y);

double	LipIntInfValueLocalCons(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y);

double	LipIntInfValueLocalConsLeftRegion(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);

double	LipIntInfValueLocalConsRightRegion(int *Dim, int *Ndata,int* Cons, double* x, double* Xd,double* y, double* Region);


void	LipIntInfComputeLipschitz(int *Dim, int *Ndata, double* x, double* y);

void	LipIntInfComputeLocalLipschitz(int *Dim, int *Ndata, double* x, double* y);

void	LipIntInfComputeLipschitzCV(int *Dim, int *Ndata, double* Xd, double* y, double* T,
			int* type, int* Cons, double* Region, double *W);

void	LipIntInfComputeLipschitzSplit(int *Dim, int *Ndata, double* Xd, double* y, double* T, double* ratio,
			int* type, int* Cons, double* Region, double *W);


void	LipIntInfSmoothLipschitz(int *Dim, int *Ndata,  double* Xd, double* y, double* T,  double* LC, 
							  int* fW, int* fC, int* fR, double* W, int* Cons, double* Region);
 // fR is 0, 1-left, 2-right

double	LipIntInfGetLipConst() ;
void	LipIntInfGetScaling(double *S) ;
int		LipIntInfComputeScaling(int *Dim, int *Ndata, double* XData, double* YData);

int		LipIntInfVerifyMonotonicity(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* LC, double* eps);
int		LipIntInfVerifyMonotonicityLeftRegion(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* Region, 
										   double* LC, double* eps);
int		LipIntInfVerifyMonotonicityRightRegion(int *Dim, int* npts, int* Cons,  double* XData, double* YData, double* Region, 
										   double* LC, double* eps);


void	LipIntInfSmoothLipschitzSimp(int *Dim, int* npts,  double* XData, double* YData, double* TData,  double* LC);
void	LipIntInfSmoothLipschitzSimpW(int *Dim, int* npts,  double* XData, double* YData, double* TData,  double* LC, double* W);




/* interface to the members of STCInterpolant class ====================================== */

// supplies the data to Interpolant and constructs the interpolant
// assuming a given Lipschitz constant, supplied by SetLipschitz
// if LipConstant was not supplied, tries to find it from the data
// assumes that all data are different. 

 int	STCBuildLipInterpolant(int *Dim, int *Ndata, double* x, double* y);

// as above, but for explicit evaluation, needs no preprocessing, but may be slower
 int	STCBuildLipInterpolantExplicit(int *Dim, int *Ndata, double* x, double* y);

// in the methods above, the coordinates of the data points in x are stored in rows

// the following methods store data in columns (like in fortran or Matlab)
// they use the transposed of the matrix x 
 int	STCBuildLipInterpolantColumn(int *Dim, int *Ndata, double* x, double* y);

// as above, but for explicit evaluation, needs no preprocessing, but may be slower
 int	STCBuildLipInterpolantExplicitColumn(int *Dim, int *Ndata, double* x, double* y);


// specify the Lipschitz constant for your function
 void	STCSetLipschitz(double* x);

// computes the value of the interpolant at any given point x
 double	STCValue( double* x );

// same but using explicit evaluation with no preprocessing
 double	STCValueExplicit( double* x );

 void	STCFreeMemory();






#ifdef __cplusplus
}
#endif
