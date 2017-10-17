#include "mex.h"
#include "liblip.h"
#include <string.h>

typedef struct _sParams
{
    char   m_pName[100];
    int    m_Dim;
    int    m_N;
    int    m_N1, m_M1;
    int    *m_pCons;
    double *m_pX;
    double *m_pXd;
    double *m_pYd;
    double m_Ratio;
    int    m_Type;
    double m_LC;
    double *m_pRegion;
    double *m_pW;
    int    m_fW;
    int    m_fC;
    int    m_fR;
    double m_eps;
    int    *m_pIndex;
//boundaries callbacks:  pass lowerBoundaryCallback and upperBoundaryCallback callback function to liblip interface
    //upper boundary
    mxArray *m_prhs_cb_UB[4];
    mxArray *m_plhs_cb_UB[1];
    double  *m_pDistParamArr_UB;

	double  *pParUB,*pParLB;
	double  *pDimUB,*pDimLB;

    //lower boundary
    mxArray *m_prhs_cb_LB[4];
    mxArray *m_plhs_cb_LB[1];
    double  *m_pDistParamArr_LB;

}sParams;

static sParams gParams;

static int boundaries_set=0;

void my_strlwr_s(char* s)
{
	int i=0;
	while(s[i]!=0) {
		if (s[i]>='A' && s[i]<='Z') s[i]=s[i]+('a'-'A');
	}
}

void debugParams(sParams *pParams)
{
    int i,j;
    if (pParams->m_pName)
        mexPrintf("Name=%s \n", pParams->m_pName);
    mexPrintf("Dim=%d \n", pParams->m_Dim);
    mexPrintf("N=%d \n", pParams->m_N);
    if (pParams->m_pXd)
    {
        mexPrintf("Xd= \n");
        for (i = 0; i < pParams->m_N; i++)
        {
            for (j = 0; j < pParams->m_Dim; j++)
                mexPrintf(" %f ", pParams->m_pXd[i*pParams->m_Dim +j]);
            printf("\n");
        }
    }
}

bool parseName(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    int buflen = 0;
    bool bRes = false;
    if (nrhs)
    {
        bRes = mxIsChar(prhs[paramNo]);
        if (bRes)
        {
            buflen = (mxGetM(prhs[paramNo]) * mxGetN(prhs[paramNo])) + 1;
            memset(pParams->m_pName, 0, buflen);
            if (mxGetString(prhs[paramNo], pParams->m_pName, buflen))
            {
                mexPrintf("error extracting function name\n", paramNo);
                bRes = false;
            } else
			     my_strlwr_s(pParams->m_pName);
        }
        else /*"Error: paramter #%d */
            mexPrintf("%d input must be function name\n", paramNo);
    }
    else 
        mexPrintf("Error: wrong number of parameters (%d). Must be 1 or more.\n", nrhs);
    return bRes;
}

bool parseDim(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) !=1 )
    {
        mexPrintf("%d input must be an integer (dimension)\n", paramNo);
        bRes = false;
    }
    else
        pParams->m_Dim = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseN(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d input must be an integer (the number of data)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_N = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseCons(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsInt32(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_Dim &&  mxGetM(prhs[paramNo])==1) || !(mxGetM(prhs[paramNo])== pParams->m_Dim &&  mxGetN(prhs[paramNo])==1))
    {
        mexPrintf("%d input must be a vector of integers size 1 x %d (monotonicity constraints)\n", paramNo, pParams->m_Dim);
        bRes = false;
    }
    else  //x check that below code works
    {
        pParams->m_pCons = (int *)mxGetPr(prhs[paramNo]);
    }
    return bRes;
}


bool parseX(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_Dim &&  mxGetM(prhs[paramNo])==1) || !(mxGetM(prhs[paramNo])== pParams->m_Dim &&  mxGetN(prhs[paramNo])==1))
    {
        mexPrintf("%d input must be a vector of doubles dimension size 1 x %d (data to interploate)\n", paramNo, pParams->m_Dim);
        bRes = false;
    }
    else
        pParams->m_pX = mxGetPr(prhs[paramNo]);
    return bRes;
}

bool parseXd(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    double d;
    if (!mxIsDouble(prhs[paramNo]) || mxGetM(prhs[paramNo])!= pParams->m_Dim ||  mxGetN(prhs[paramNo]) != pParams->m_N)
    {
        mexPrintf("%d input must be an array of doubles size %d x %d (data 'Xd')\n", paramNo, pParams->m_Dim, pParams->m_N);
        bRes = false;
    }
    else
    {
        pParams->m_pXd = mxGetPr(prhs[paramNo]);
        d = mxGetScalar(prhs[paramNo]);
    }
    return bRes;
}

bool parseY(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_N && mxGetM(prhs[paramNo])==1) ||  !(mxGetM(prhs[paramNo])== pParams->m_N && mxGetN(prhs[paramNo])==1) )
    {
        mexPrintf("%d input must be must be a vector of doubles dimension size 1 x %d (output of the function on data 'X')\n", paramNo, pParams->m_N);
        bRes = false;
    }
	else {
        pParams->m_pYd = mxGetPr(prhs[paramNo]);
		pParams->m_M1=mxGetM(prhs[paramNo]);
		pParams->m_N1=mxGetN(prhs[paramNo]);
	}
    return bRes;
}

bool parseRatio(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d dst input must be an integer (ratio)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_Ratio = mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseType(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsInt32(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d dst input must be an integer (type)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_Type = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseLC(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]))
    {
        mexPrintf("%d input must be a double (Lipschitz constant)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_LC = mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseRegion(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_Dim &&  mxGetM(prhs[paramNo])==1) || !(mxGetM(prhs[paramNo])== pParams->m_Dim &&  mxGetN(prhs[paramNo])==1))

  //  if (!mxIsDouble(prhs[paramNo]) || mxGetN(prhs[paramNo])!= pParams->m_Dim ||  mxGetM(prhs[paramNo])!=1)
    {
        mexPrintf("%d input must be a vector of doubles dimension size 1 x %d (Region)\n", paramNo, pParams->m_Dim);
        bRes = false;
    }
    else
        pParams->m_pRegion = mxGetPr(prhs[paramNo]);
    return bRes;
}

bool parseW(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_N &&  mxGetM(prhs[paramNo])==1) || !(mxGetM(prhs[paramNo])== pParams->m_N&&  mxGetN(prhs[paramNo])==1))
    {
        mexPrintf("%d input must be a vector of doubles dimension size 1 x %d (Weights)\n", paramNo, pParams->m_N);
        bRes = false;
    }
    else
        pParams->m_pW = mxGetPr(prhs[paramNo]);
    return bRes;
}

bool parsefW(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsInt32(prhs[paramNo]))
    {
        mexPrintf("%d input must be an int(fW)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_fW = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parsefC(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsInt32(prhs[paramNo]))
    {
        mexPrintf("%d input must be an int(fC)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_fC = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parsefR(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsInt32(prhs[paramNo]))
    {
        mexPrintf("%d input must be an int(fR)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_fR = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseEps(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]))
    {
        mexPrintf("%d input must be a double (eps)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_eps = mxGetScalar(prhs[paramNo]);
    return bRes;
}


bool parseIndex(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pIndex = NULL;
    if (!mxIsInt32(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_N &&  mxGetM(prhs[paramNo])==1) || !(mxGetM(prhs[paramNo])== pParams->m_N&&  mxGetN(prhs[paramNo])==1))

 //   if (!mxIsInt32(prhs[paramNo]) || mxGetN(prhs[paramNo])!= pParams->m_N ||  mxGetM(prhs[paramNo])!=1)
    {
        mexPrintf("%d input must be a vector of integers size 1 x %d (index)\n", paramNo, pParams->m_N);
        bRes = false;
    }
    else  //x check that below code works
    {
        pParams->m_pIndex = (int *)mxGetPr(prhs[paramNo]);
    }
    return bRes;
}

void freeBoundaryMem(void)
{
//upper boundary
    if (gParams.m_prhs_cb_UB[0])
        mxDestroyArray(gParams.m_prhs_cb_UB[0]);

    if (gParams.m_prhs_cb_UB[1])
        mxDestroyArray(gParams.m_prhs_cb_UB[1]);

    if (gParams.m_prhs_cb_UB[2])
        mxDestroyArray(gParams.m_prhs_cb_UB[2]);

    if (gParams.m_prhs_cb_UB[3])
        mxDestroyArray(gParams.m_prhs_cb_UB[3]);

    gParams.m_prhs_cb_UB[0] = gParams.m_prhs_cb_UB[1] = gParams.m_prhs_cb_UB[2] = gParams.m_prhs_cb_UB[3] =NULL;
//lower boudary
    if (gParams.m_prhs_cb_LB[0])
        mxDestroyArray(gParams.m_prhs_cb_LB[0]);

    if (gParams.m_prhs_cb_LB[1])
        mxDestroyArray(gParams.m_prhs_cb_LB[1]);

    if (gParams.m_prhs_cb_LB[2])
        mxDestroyArray(gParams.m_prhs_cb_LB[2]);

    if (gParams.m_prhs_cb_LB[3])
        mxDestroyArray(gParams.m_prhs_cb_LB[3]);

    gParams.m_prhs_cb_LB[0] = gParams.m_prhs_cb_LB[1] = gParams.m_prhs_cb_LB[2] = gParams.m_prhs_cb_LB[3] =NULL;

}  

bool parseBoundaryFuncRef(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    int buflen = 0;
    bool bRes = true;
    freeBoundaryMem();

    gParams.m_prhs_cb_UB[0] = mxDuplicateArray(prhs[paramNo]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_UB[0]);

    gParams.m_prhs_cb_UB[1] = mxCreateDoubleMatrix(1, 100 , mxREAL);
    gParams.m_pDistParamArr_UB = mxGetPr(gParams.m_prhs_cb_UB[1]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_UB[1]);
    
    gParams.m_prhs_cb_UB[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    gParams.pDimUB = mxGetPr(gParams.m_prhs_cb_UB[2]);
    gParams.pDimUB[0] = gParams.m_Dim;
    mexMakeArrayPersistent(gParams.m_prhs_cb_UB[2]);

    gParams.m_prhs_cb_UB[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    gParams.pParUB = mxGetPr(gParams.m_prhs_cb_UB[3]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_UB[3]);

//lower boudary
    paramNo++;

    gParams.m_prhs_cb_LB[0] = mxDuplicateArray(prhs[paramNo]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_LB[0]);

    gParams.m_prhs_cb_LB[1] = mxCreateDoubleMatrix(1, 100, mxREAL);
    gParams.m_pDistParamArr_LB = mxGetPr(gParams.m_prhs_cb_LB[1]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_LB[1]);
    
    gParams.m_prhs_cb_LB[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    gParams.pDimLB = mxGetPr(gParams.m_prhs_cb_LB[2]);
    gParams.pDimLB[0] = gParams.m_Dim;
    mexMakeArrayPersistent(gParams.m_prhs_cb_LB[2]);

	gParams.m_prhs_cb_LB[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    gParams.pParLB = mxGetPr(gParams.m_prhs_cb_LB[3]);
    mexMakeArrayPersistent(gParams.m_prhs_cb_LB[3]);

    return bRes;
}



/*
double upperBoundaryCallback(double* p, int dim)
{ // example: multivariate normal distribution
    
    double res;
    memcpy(gParams.m_pDistParamArr_UB, p, dim * sizeof(double));
    if (mexCallMATLAB(1 , gParams.m_plhs_cb_UB, 3, gParams.m_prhs_cb_LB, "feval")==0)
    {
        res = mxGetScalar(gParams.m_plhs_cb_UB[0]);
    }
    else 
        mexPrintf("callback error \n");
    mxDestroyArray(gParams.m_plhs_cb_UB[0]);
    return res;
}


double lowerBoundaryCallback(double* p, int dim, double* par)
{ // example: multivariate normal distribution
    
    double res;
    memcpy(gParams.m_pDistParamArr_LB, p, dim * sizeof(double));

    if (mexCallMATLAB(1 , gParams.m_plhs_cb_LB, 3, gParams.m_prhs_cb_LB, "feval")==0)
    {
        res = mxGetScalar(gParams.m_plhs_cb_LB[0]);
    }
    else 
        mexPrintf("callback error \n");
    mxDestroyArray(gParams.m_plhs_cb_LB[0]);
    return res;
}
*/

double upwrap(int* dim, double* x, double* p) {
    double res;

	if (!boundaries_set) return 10e20;
    memcpy(gParams.m_pDistParamArr_UB, p, *dim * sizeof(double));
	gParams.pParUB[0]= *p; 
	gParams.pDimUB[0]= *dim;

    if (mexCallMATLAB(1 , gParams.m_plhs_cb_UB, 4, gParams.m_prhs_cb_UB, "feval")==0)
    {
        res = mxGetScalar(gParams.m_plhs_cb_UB[0]);
    }
    else 
        mexPrintf("callback error \n");
    mxDestroyArray(gParams.m_plhs_cb_UB[0]);
    return res;
}

double lowwrap(int* dim, double* x, double* p) {
    double res;
	if (!boundaries_set) return -10e20;
    memcpy(gParams.m_pDistParamArr_LB, p, *dim * sizeof(double));
	gParams.pParLB[0]= *p; 
	gParams.pDimLB[0]= *dim;

    if (mexCallMATLAB(1 , gParams.m_plhs_cb_LB, 4, gParams.m_prhs_cb_LB, "feval")==0)
    {
        res = mxGetScalar(gParams.m_plhs_cb_LB[0]);
    }
    else 
        mexPrintf("callback error \n");
    mxDestroyArray(gParams.m_plhs_cb_LB[0]);
    return res;
}

bool parseParams(int nrhs, const mxArray *prhs[], sParams *pParams)
{
    bool bRes = true;
    int paramNo = 0;
    bRes = parseName(nrhs, prhs, pParams, paramNo);
    paramNo++;
    if (bRes)
    {
        if (strcmp(pParams->m_pName, "value") == 0 || strcmp(pParams->m_pName, "infvalue") == 0)
        {       
            if (nrhs == 7 || nrhs == 8)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                pParams->m_pIndex = NULL;
                if (nrhs == 8)
                    bRes = parseIndex(nrhs, prhs, pParams, paramNo) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 7 or 8.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "valueauto") == 0 || strcmp(pParams->m_pName, "infvalueauto") == 0)
        {       
            if (nrhs == 6 || nrhs == 7)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                pParams->m_pIndex = NULL;
                if (nrhs == 7)
                    bRes = parseIndex(nrhs, prhs, pParams, paramNo) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 6 or 7.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "valuecons") == 0 || strcmp(pParams->m_pName, "infvaluecons") == 0)
        {       
            if (nrhs == 8 || nrhs == 9)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                pParams->m_pIndex = NULL;
                if (nrhs == 9)
                    bRes = parseIndex(nrhs, prhs, pParams, paramNo) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 8 or 9.\n", nrhs);
                bRes = false;
            }
        }
        else  if ((strcmp(pParams->m_pName, "valueconsleftregion") == 0) || (strcmp(pParams->m_pName, "valueconsrightregion") == 0) || 
            (strcmp(pParams->m_pName, "infvalueconsleftregion") == 0) || (strcmp(pParams->m_pName, "valueinfconsrightregion") == 0))
        {       
            if (nrhs == 9 || nrhs == 10)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRegion(nrhs, prhs, pParams, paramNo++) && bRes;
                pParams->m_pIndex = NULL;
                if (nrhs == 10)
                    bRes = parseIndex(nrhs, prhs, pParams, paramNo) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 9 or 10.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "valuelocal") == 0 || strcmp(pParams->m_pName, "infvaluelocal") == 0)
        {       
            if (nrhs == 6)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 7.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "valuelocalcons") == 0 || strcmp(pParams->m_pName, "infvaluelocalcons") == 0)
        {       
            if (nrhs == 7)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 7.\n", nrhs);
                bRes = false;
            }
        }
        else  if ((strcmp(pParams->m_pName, "valuelocalconsleftregion") == 0) || (strcmp(pParams->m_pName, "valuelocalconsrightregion") == 0) ||
            (strcmp(pParams->m_pName, "infvaluelocalconsleftregion") == 0) || (strcmp(pParams->m_pName, "infvaluelocalconsrightregion") == 0))
        {       
            if (nrhs == 8)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRegion(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 8.\n", nrhs);
                bRes = false;
            }
        }
        else  if ((strcmp(pParams->m_pName, "computelipschitz") == 0) || (strcmp(pParams->m_pName, "computelocallipschitz") == 0) ||
            (strcmp(pParams->m_pName, "infcomputelipschitz") == 0) || (strcmp(pParams->m_pName, "infcomputelocallipschitz") == 0))
        {       
            if (nrhs == 5)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 5.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "computelipschitzcv") == 0 || strcmp(pParams->m_pName, "infcomputelipschitzcv") == 0)
        {       
            if (nrhs == 7 || nrhs == 8 || nrhs == 9)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseType(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs >= 6)
                    bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs >= 7)
                    bRes = parseRegion(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs == 8)
                    bRes = parseW(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 7 or 8 or 9.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "computelipschitzsplit") == 0 || strcmp(pParams->m_pName, "infcomputelipschitzsplit") == 0)
        {       
            if ( nrhs == 8 || nrhs == 9 || nrhs == 10 )
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRatio(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseType(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs >= 7)
                    bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs >= 8)
                    bRes = parseRegion(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs == 9)
                    bRes = parseW(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 8 or 9 or 10.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "smoothlipschitz") == 0 || strcmp(pParams->m_pName, "infsmoothlipschitz") == 0)
        {       
            if ( nrhs == 10 || nrhs == 11 || nrhs == 12 )
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parsefW(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parsefC(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parsefR(nrhs, prhs, pParams, paramNo++) && bRes;

                if (nrhs >=  paramNo  && pParams->m_fW)
                    bRes = parseW(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs >=  paramNo  && pParams->m_fC)
                    bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                if (nrhs ==  paramNo  && pParams->m_fR)
                    bRes = parsefR(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 10 or 11 or 12.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "getlipconst") == 0 || strcmp(pParams->m_pName, "infgetlipconst") == 0)
        {       
            if ( nrhs == 1 )
                bRes = true;
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 1.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "getscaling") == 0 || strcmp(pParams->m_pName, "infgetscaling") == 0)
        {       
            if ( nrhs == 1 )
                bRes = true;
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 1.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "computescaling") == 0 || strcmp(pParams->m_pName, "infcomputescaling") == 0)
        {       
            if ( nrhs == 5 )
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 5.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "verifymonotonicity") == 0 || strcmp(pParams->m_pName, "infverifymonotonicity") == 0 )
        {       
            if ( nrhs == 8)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseEps(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 8.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "verifymonotonicityleftregion") == 0 || strcmp(pParams->m_pName, "verifymonotonicityrightregion") == 0 ||
            strcmp(pParams->m_pName, "infverifymonotonicityleftregion") == 0 || strcmp(pParams->m_pName, "infverifymonotonicityrightregion") == 0)
        {       
            if ( nrhs == 9)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseCons(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRegion(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseEps(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 9.\n", nrhs);
                bRes = false;
            }
        } 
        else  if (strcmp(pParams->m_pName, "stcbuildlipinterpolant") == 0 || strcmp(pParams->m_pName, "stcbuildlipinterpolantexplicit") == 0 )
        {       
            if ( nrhs == 5 || nrhs == 6)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseN(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseXd(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseY(nrhs, prhs, pParams, paramNo++) && bRes;
				if( nrhs == 6 )
				  bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes; else
				  pParams->m_LC=0;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 5 or 6.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "stcsetlipschitz") == 0 )
        {       
            if ( nrhs == 2)
            {
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 2.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "stcvalue") == 0 || strcmp(pParams->m_pName, "stcvalueexplicit") == 0)
        {       
            if ( nrhs == 2)
            {
                bRes = parseX(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 2.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "stcfreememory") == 0)
        {       
            if ( nrhs == 1)
            {
                bRes = true;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 1.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "setbounds") == 0)
        {       
            if ( nrhs == 3)
            {
                bRes = parseBoundaryFuncRef(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = true;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 3.\n", nrhs);
                bRes = false;
            }
        }
        else  if (strcmp(pParams->m_pName, "clearbounds") == 0)
        {       
            if ( nrhs == 1)
            {
                bRes = true;
            }
            else 
            {
                mexPrintf("Error: wrong number of parameters (%d). Must be 0.\n", nrhs);
                bRes = false;
            }
        }
        else
        {
            mexPrintf("Error: unknown function name '%s' \n", pParams->m_pName);
            bRes = false;
        }
    }
    return bRes;
}

void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    bool bRes = true;
    double *pRes = NULL;
    mexAtExit(freeBoundaryMem);
// initialization of parameters
/*	if(!boundaries_set) {
      gParams.m_prhs_cb_UB[0] = gParams.m_prhs_cb_UB[1] = gParams.m_prhs_cb_UB[2] = NULL;
      gParams.m_prhs_cb_LB[0] = gParams.m_prhs_cb_LB[1] = gParams.m_prhs_cb_LB[2] = NULL;
	}
*/
	gParams.m_pCons=NULL;
	gParams.m_pIndex=NULL;
	gParams.m_pRegion=NULL;
	gParams.m_pW=NULL;
	gParams.m_Type=0;
	gParams.m_pIndex=NULL;

    bRes = parseParams(nrhs, prhs, &gParams);
    if (bRes)
    {
        if (strcmp(gParams.m_pName, "value")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValue(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
        }

		else if (strcmp(gParams.m_pName, "stcvalue")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = STCValue(gParams.m_pX);
        }
        else if (strcmp(gParams.m_pName, "stcvalueexplicit")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = STCValueExplicit(gParams.m_pX);
        }

        if (strcmp(gParams.m_pName, "infvalue")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValue(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "valueauto")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueAuto(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "infvalueauto")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueAuto(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pIndex);
        }

        else if (strcmp(gParams.m_pName, "valuecons")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "infvaluecons")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "valueconsleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "infvalueconsleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "valueconsrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "infvalueconsrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
        }
        else if (strcmp(gParams.m_pName, "valuelocal")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueLocal(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "infvaluelocal")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueLocal(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "valuelocalcons")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueLocalCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "infvaluelocalcons")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueLocalCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "valuelocalconsleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueLocalConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "infvaluelocalconsleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueLocalConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "valuelocalconsrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntValueLocalConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "infvaluelocalconsrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfValueLocalConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "computelipschitz")==0)
        {
            LipIntComputeLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "infcomputelipschitz")==0)
        {
            LipIntInfComputeLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "computelocallipschitz")==0)
        {
            LipIntComputeLocalLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "infcomputelocallipschitz")==0)
        {
            LipIntInfComputeLocalLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "computelipschitzcv")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntComputeLipschitzCV(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
        }
        else if (strcmp(gParams.m_pName, "infcomputelipschitzcv")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntInfComputeLipschitzCV(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
        }
        else if (strcmp(gParams.m_pName, "computelipschitzsplit")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntComputeLipschitzSplit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_Ratio, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
        }
        else if (strcmp(gParams.m_pName, "infcomputelipschitzsplit")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntInfComputeLipschitzSplit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_Ratio, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
        }
        else if (strcmp(gParams.m_pName, "smoothlipschitz")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntSmoothLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_LC, &gParams.m_fW, &gParams.m_fC, &gParams.m_fR, gParams.m_pW, gParams.m_pCons, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "infsmoothlipschitz")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_M1, gParams.m_N1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntInfSmoothLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, pRes, &gParams.m_LC, &gParams.m_fW, &gParams.m_fC, &gParams.m_fR, gParams.m_pW, gParams.m_pCons, gParams.m_pRegion);
        }
        else if (strcmp(gParams.m_pName, "getlipconst")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntGetLipConst();
        }
        else if (strcmp(gParams.m_pName, "infgetlipconst")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfGetLipConst();
        }
        else if (strcmp(gParams.m_pName, "getscaling")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, gParams.m_Dim, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntGetScaling(pRes);
        }
        else if (strcmp(gParams.m_pName, "infgetscaling")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, gParams.m_Dim, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntInfGetScaling(pRes);
        }
        else if (strcmp(gParams.m_pName, "computescaling")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, gParams.m_Dim, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntComputeScaling(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
			LipIntGetScaling(pRes);
        }
        else if (strcmp(gParams.m_pName, "infcomputescaling")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, gParams.m_Dim, mxREAL);
            pRes = mxGetPr(plhs[0]);
            LipIntComputeScaling(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
		    LipIntInfGetScaling(pRes);
        }
        else if (strcmp(gParams.m_pName, "verifymonotonicity")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntVerifyMonotonicity(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "infverifymonotonicity")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfVerifyMonotonicity(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "verifymonotonicityleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntVerifyMonotonicityLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "infverifymonotonicityleftregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfVerifyMonotonicityLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "verifymonotonicityrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntVerifyMonotonicityRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "infverifymonotonicityrightregion")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipIntInfVerifyMonotonicityRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
        }
        else if (strcmp(gParams.m_pName, "stcbuildlipinterpolant")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            STCSetLipschitz(&gParams.m_LC);
            *pRes = STCBuildLipInterpolant(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "stcbuildlipinterpolantexplicit")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
			STCSetLipschitz(&gParams.m_LC);
            *pRes = STCBuildLipInterpolantExplicit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
        }
        else if (strcmp(gParams.m_pName, "stcsetlipschitz")==0)
        {
            STCSetLipschitz(&gParams.m_LC);
        }
        else if (strcmp(gParams.m_pName, "stcfreememory")==0)
        {
            freeBoundaryMem();
            STCFreeMemory();
        }
        else if (strcmp(gParams.m_pName, "setbounds")==0)
        {
			boundaries_set=1;
			// now set members of sli and slii
			LipIntSetBounds(&lowwrap,&upwrap);
         }
        else if (strcmp(gParams.m_pName, "clearbounds")==0)
        {
			LipIntClearBounds();
			boundaries_set=0;
         }
    }

    if (!pRes && nlhs)
    {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        pRes = mxGetPr(plhs[0]);
        *pRes = 0;
    }

}
