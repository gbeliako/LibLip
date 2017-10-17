#include "maplec.h"
#include "liblip.h"

#ifndef bool
#define bool int
#endif
#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef NULL
#define NULL 0
#endif

#ifdef _WIN32
#define DllExport   __declspec( dllexport )
#else
#define DllExport
#endif

typedef struct _sParams
{
    int    m_Dim;
    int    m_N;
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
    double *m_pOutDoubleArr;
//boundaries callbacks:  pass lowerBoundaryCallback and upperBoundaryCallback callback function to liblip interface
//use SetBoundaries to callbacks
    MKernelVector m_MapleKernelVec;
    ALGEB   m_CallBackFunc_UB;
    ALGEB   m_CallBackFunc_LB;
}sParams;

static sParams gParams;
static int boundaries_set=0;
static double temp[100]; // just some working vector


bool parseDim(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Dim = MapleToInteger32(kv,(ALGEB) args[paramNo]);

// initialization of parameters
	if(!boundaries_set) {
		pParams->m_CallBackFunc_UB=NULL;
		pParams->m_CallBackFunc_LB=NULL;
	}
	pParams->m_pCons=NULL;
	pParams->m_pIndex=NULL;
	pParams->m_pRegion=NULL;
	pParams->m_pW=NULL;
	pParams->m_Type=0;
	pParams->m_pOutDoubleArr=NULL;
	pParams->m_pIndex=NULL;

    return bRes;
}

bool parseN(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_N = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}

bool parseCons(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pCons = (INTEGER32 *) RTableDataBlock( kv, (ALGEB) args[paramNo]);
    return bRes;
}


bool parseX(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pX = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseXd(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pXd = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseY(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pYd = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseRatio(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Ratio = MapleToFloat64( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseType(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Type = MapleToInteger32( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseLC(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_LC = MapleToFloat64( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseRegion(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pRegion = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseW(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pW = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parsefW(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_fW = MapleToInteger32( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parsefC(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_fC = MapleToInteger32( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parsefR(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_fR = MapleToInteger32( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseEps(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_eps = MapleToFloat64( kv, (ALGEB) args[paramNo]);
    return bRes;
}


bool parseIndex(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pIndex = (INTEGER32 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseOutDoubleArr(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pOutDoubleArr = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseBoundaryCallbacks(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{ 
    bool bRes = true;
    pParams->m_MapleKernelVec = kv;
	if( IsMapleProcedure(kv,(ALGEB)args[paramNo]) )
        pParams->m_CallBackFunc_UB = (ALGEB)args[paramNo];
    else 
    {
        MapleRaiseError(kv,"procedure expected for upper boundary");
        bRes = false;
    }
    paramNo++;

    if( IsMapleProcedure(kv,(ALGEB)args[paramNo]) )
        pParams->m_CallBackFunc_LB = (ALGEB)args[paramNo];
    else 
    {
        MapleRaiseError(kv,"procedure expected for lower boundary");
        bRes = false;
    }
    return bRes;
}




DllExport ALGEB MWRAP_LipIntValue(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    int paramNo = 1;
    bool bRes=true;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 6 || nargs == 7)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 7)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntValue(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 6 or 7.\n");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValue(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 6 || nargs == 7)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 7)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntInfValue(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 6 or 7.\n");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}



DllExport ALGEB MWRAP_LipIntValueAuto(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 5 || nargs == 6)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 6)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntValueAuto(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 5 or 6");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueAuto(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 5 || nargs == 6)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 6)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntInfValueAuto(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 5 or 6");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}


DllExport ALGEB MWRAP_LipIntValueCons(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 7 || nargs == 8)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 8)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntValueCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7 or 8");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueCons(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 7 || nargs == 8)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 8)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
        if (bRes)
            res = LipIntInfValueCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pIndex);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7 or 8");
        bRes = false;
    }
    return ToMapleFloat(kv, res);
}

bool parseRegionParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 8 || nargs == 9)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
		bRes = parseRegion(kv, args, pParams, paramNo++) && bRes;
        pParams->m_pIndex = NULL;
        if (nargs == 9)
            bRes = parseIndex(kv, args, pParams, paramNo) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7 or 8");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntValueConsLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseRegionParams(kv, args))
        res = LipIntValueConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntValueConsRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseRegionParams(kv, args))
        res = LipIntValueConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueConsLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseRegionParams(kv, args))
        res = LipIntInfValueConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueConsRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseRegionParams(kv, args))
        res = LipIntInfValueConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, gParams.m_pRegion, gParams.m_pIndex);
    return ToMapleFloat(kv, res);
}

bool parseLipIntValueLocalParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 5)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 5");
        bRes = false;
    }
    return bRes;
}


DllExport ALGEB MWRAP_LipIntValueLocal(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalParams(kv, args))
        res = LipIntValueLocal(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
    return ToMapleFloat(kv, res);
}


DllExport ALGEB MWRAP_LipIntInfValueLocal(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalParams(kv, args))
        res = LipIntInfValueLocal(&gParams.m_Dim, &gParams.m_N, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
    return ToMapleFloat(kv, res);
}

bool parseLipIntValueLocalConsParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 6)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 6");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntValueLocalCons(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsParams(kv, args))
        res = LipIntValueLocalCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueLocalCons(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsParams(kv, args))
        res = LipIntInfValueLocalCons(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd);
    return ToMapleFloat(kv, res);
}


bool parseLipIntValueLocalConsRegionParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 7)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRegion(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntValueLocalConsLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsRegionParams(kv, args))
        res = LipIntValueLocalConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntValueLocalConsRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsRegionParams(kv, args))
        res = LipIntValueLocalConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueLocalConsLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsRegionParams(kv, args))
        res = LipIntInfValueLocalConsLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfValueLocalConsRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntValueLocalConsRegionParams(kv, args))
        res = LipIntInfValueLocalConsRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pX, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion);
    return ToMapleFloat(kv, res);
}


bool parseLipIntComputeLipschitzParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 4)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 4");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntComputeLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzParams(kv, args))
        LipIntComputeLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntComputeLocalLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzParams(kv, args))
        LipIntComputeLocalLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntInfComputeLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzParams(kv, args))
        LipIntInfComputeLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntInfComputeLocalLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzParams(kv, args))
        LipIntInfComputeLocalLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
    return ToMapleNULL(kv);
}

bool parseLipIntComputeLipschitzCVParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 6 || nargs == 7 || nargs == 8 || nargs == 9) 
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        bRes = parseType(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= 7)
            bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= 8)
            bRes = parseRegion(kv, args, pParams, paramNo++) && bRes;
        if (nargs == 9)
            bRes = parseW(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 6 or 7 or 8 or 9");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntComputeLipschitzCV(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzCVParams(kv, args))
    {
            LipIntComputeLipschitzCV(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntInfComputeLipschitzCV(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzCVParams(kv, args))
    {
            LipIntInfComputeLipschitzCV(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
    }
    return ToMapleNULL(kv);
}

bool parseLipIntComputeLipschitzSplitParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 7 || nargs == 8 || nargs == 9 || nargs == 10 )//fixed
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRatio(kv, args, pParams, paramNo++) && bRes;
        bRes = parseType(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= 8)
            bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= 9)
            bRes = parseRegion(kv, args, pParams, paramNo++) && bRes;
        if (nargs == 10)
            bRes = parseW(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7 or 8 or 9 or 10");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntComputeLipschitzSplit(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzSplitParams(kv, args))
    {
        LipIntComputeLipschitzSplit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_Ratio, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
    }
    return ToMapleNULL(kv);
}



DllExport ALGEB MWRAP_LipIntInfComputeLipschitzSplit(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeLipschitzSplitParams(kv, args))
    {
        LipIntInfComputeLipschitzSplit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_Ratio, &gParams.m_Type, gParams.m_pCons, gParams.m_pRegion, gParams.m_pW);
    }
    return ToMapleNULL(kv);
}

bool parseLipIntSmoothLipschitzParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 9 || nargs == 10 || nargs == 11 || nargs == 12 ) //fixed
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        bRes = parsefW(kv, args, pParams, paramNo++) && bRes;
        bRes = parsefC(kv, args, pParams, paramNo++) && bRes;
        bRes = parsefR(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= paramNo && pParams->m_fW)
            bRes = parseW(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= paramNo && pParams->m_fC)
            bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        if (nargs >= paramNo && pParams->m_fR)
            bRes = parsefR(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 9 or 10 or 11 or 12");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntSmoothLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntSmoothLipschitzParams(kv, args))
            LipIntSmoothLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_LC, &gParams.m_fW, &gParams.m_fC, &gParams.m_fR, gParams.m_pW, gParams.m_pCons, gParams.m_pRegion);
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntInfSmoothLipschitz(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntSmoothLipschitzParams(kv, args))
        LipIntSmoothLipschitz(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd, gParams.m_pOutDoubleArr, &gParams.m_LC, &gParams.m_fW, &gParams.m_fC, &gParams.m_fR, gParams.m_pW, gParams.m_pCons, gParams.m_pRegion);
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_LipIntGetLipConst(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = LipIntGetLipConst();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0.\n");
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfGetLipConst(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = LipIntInfGetLipConst();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0.\n");
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntGetScaling(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
    {
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        LipIntGetScaling(gParams.m_pOutDoubleArr);
    }
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0");
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_LipIntInfGetScaling(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
    {
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        LipIntInfGetScaling(gParams.m_pOutDoubleArr);
    }
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0");
    return ToMapleNULL(kv);
}



bool parseLipIntComputeScalingParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 5 )
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
		bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 5");
        bRes = false;
    }
    return bRes;
}


DllExport ALGEB MWRAP_LipIntComputeScaling(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeScalingParams(kv, args))
    {
            LipIntComputeScaling(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
			LipIntGetScaling(gParams.m_pOutDoubleArr);
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LipIntInfComputeScaling(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntComputeScalingParams(kv, args))
    {
            LipIntInfComputeScaling(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
			LipIntInfGetScaling(gParams.m_pOutDoubleArr);
    }
    return ToMapleNULL(kv);
}


bool parseLipIntVerifyMonotonicityParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 7)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        bRes = parseEps(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 7");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_LipIntVerifyMonotonicity(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityParams(kv, args))
        res = LipIntVerifyMonotonicity(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfVerifyMonotonicity(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityParams(kv, args))
        res = LipIntInfVerifyMonotonicity(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, &gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

bool parseLipIntVerifyMonotonicityRegionParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 8)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseCons(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRegion(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        bRes = parseEps(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 8");
        bRes = false;
    }
    return bRes;
}


DllExport ALGEB MWRAP_LipIntVerifyMonotonicityLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityRegionParams(kv, args))
        res = LipIntVerifyMonotonicityLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntVerifyMonotonicityRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityRegionParams(kv, args))
        res = LipIntVerifyMonotonicityRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfVerifyMonotonicityLeftRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityRegionParams(kv, args))
        res = LipIntInfVerifyMonotonicityLeftRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_LipIntInfVerifyMonotonicityRightRegion(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseLipIntVerifyMonotonicityRegionParams(kv, args))
        res = LipIntInfVerifyMonotonicityRightRegion(&gParams.m_Dim, &gParams.m_N, gParams.m_pCons, gParams.m_pXd, gParams.m_pYd, gParams.m_pRegion,&gParams.m_LC, &gParams.m_eps);
    return ToMapleFloat(kv, res);
}

bool parseSTCBuildLipInterpolantParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 5 || nargs == 4)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseN(kv, args, pParams, paramNo++) && bRes;
        bRes = parseXd(kv, args, pParams, paramNo++) && bRes;
        bRes = parseY(kv, args, pParams, paramNo++) && bRes;
		if( nargs == 5 )
		  bRes = parseLC(kv, args, pParams, paramNo++) && bRes; else
		pParams->m_LC=0;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 4 or 5");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_STCBuildLipInterpolant(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
	if (parseSTCBuildLipInterpolantParams(kv, args)) {
		STCSetLipschitz(&gParams.m_LC);
        res = STCBuildLipInterpolant(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
	}
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_STCBuildLipInterpolantExplicit(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
	if (parseSTCBuildLipInterpolantParams(kv, args)){
		STCSetLipschitz(&gParams.m_LC);
        res = STCBuildLipInterpolantExplicit(&gParams.m_Dim, &gParams.m_N, gParams.m_pXd, gParams.m_pYd);
	}
    return ToMapleFloat(kv, res);
}


DllExport ALGEB MWRAP_STCSetLipschitz(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 1)
    {
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
            STCSetLipschitz(&gParams.m_LC);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 1");
        bRes = false;
    }
    return ToMapleNULL(kv);
}


bool parseSTCValueParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 1)
    {
        bRes = parseX(kv, args, pParams, paramNo++) && bRes;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 1");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_STCValue(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseSTCValueParams(kv, args))
        res = STCValue(gParams.m_pX);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_STCValueExplicit(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    sParams *pParams = &gParams;
    if (parseSTCValueParams(kv, args))
        res = STCValueExplicit(gParams.m_pX);
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_STCFreeMemory(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 0) //fixed
    {
        STCFreeMemory();
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0");
        bRes = false;
    }
    return ToMapleNULL(kv);
}



// fix LipConst
double upperBoundaryCallback(double* p, int dim)
{ 
    return EvalhfMapleProc(gParams.m_MapleKernelVec,(ALGEB)gParams.m_CallBackFunc_UB, dim, (p-1));
}

double lowerBoundaryCallback(double* p, int dim)
{ 
    return EvalhfMapleProc(gParams.m_MapleKernelVec,(ALGEB)gParams.m_CallBackFunc_LB, dim, (p-1));
}
double upwrap(int* dim, double* x, double* p) {
	int i;
	if (!boundaries_set) return 10e20;
	for(i=0;i<*dim;i++) temp[i]=x[i];
	temp[*dim]=*p;
	return upperBoundaryCallback(temp,(*dim)+1);
}
double lowwrap(int* dim, double* x, double* p) {
	int i;
	if (!boundaries_set) return -10e20;
	for(i=0;i<*dim;i++) temp[i]=x[i];
	temp[*dim]=*p;
	return lowerBoundaryCallback(temp,(*dim)+1);
}
DllExport ALGEB MWRAP_SetBoundaries(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 2)
    {
        bRes = parseBoundaryCallbacks(kv, args, pParams, paramNo++) && bRes;
		if(bRes) { 
			boundaries_set=1;
			// now set members of sli and slii
			LipIntSetBounds(&lowwrap,&upwrap);
		}
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 2");
        bRes = false;
    }
    return ToMapleNULL(kv);
}



DllExport ALGEB MWRAP_ClearBoundaries(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if ( nargs == 0)
    {
		boundaries_set=0;
		pParams->m_CallBackFunc_UB=NULL;
		pParams->m_CallBackFunc_LB=NULL;  
		LipIntClearBounds();
		// now clear members of sli and slii
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters. Must be 0");
        bRes = false;
    }
    return ToMapleNULL(kv);
}


