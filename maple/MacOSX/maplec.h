#ifndef _MAPLE_EXTERNAL_C_H_
#define _MAPLE_EXTERNAL_C_H_

/*
**      Windows DLL Magic
*/
/* These are used to provide special annotations
   of the .dll exports on Windows. */
#if !defined(EXT_DECL)
#  if defined _MSC_VER || defined __WATCOMC__
#    ifdef WMI_WINNT
#      define EXT_DECL __declspec(dllexport)
#    else
#      define EXT_DECL __declspec(dllimport)
#    endif /* WMI_WINNT */
#  else
#    define EXT_DECL
#  endif /* _MSC_VER */
#endif /* EXT_DECL */

#if !defined(M_DECL)
#  if (defined _MSC_VER || defined __WATCOMC__) && !defined(__LLVM_COMPILER)
#    define M_DECL __stdcall
#  else
#    define M_DECL
#  endif /* _MSC_VER */
#endif /* M_DECL */
#if !defined(M_CDECL)
#  if defined _MSC_VER || defined __WATCOMC__
#    define M_CDECL __cdecl
#  else
#    define M_CDECL
#  endif /* _MSC_VER */
#endif /* M_CDECL */

#include "mplshlib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------- */
/*                     Initialize a new Maple Session                     */
/* ---------------------------------------------------------------------- */

/* CallBack structure used by StartMaple for I/O */
typedef struct {
    void (M_DECL *textCallBack) ( void *data, int tag, const char *output );
    void (M_DECL *errorCallBack) ( void *data, M_INT offset, const char *msg );
    void (M_DECL *statusCallBack) ( void *data, long kilobytesUsed,
                                       long kilobytesAlloc, double cpuTime );
    char * (M_DECL *readLineCallBack) ( void *data, M_BOOL debug );
    M_BOOL (M_DECL *redirectCallBack) ( void *data, const char *name,
                                              const char *mode ) ;
    char * (M_DECL *streamCallBack) ( void *data, const char *stream,
                                         int nargs, char **args );
    M_BOOL (M_DECL *queryInterrupt) ( void *data );
    char * (M_DECL *callBackCallBack) ( void *data, char *output );
} MCallBackVectorDesc, *MCallBackVector;

/* Possible values for the tag parameter to the textCallBack function. */
#define MAPLE_TEXT_DIAG 1
#define MAPLE_TEXT_MISC 2
#define MAPLE_TEXT_OUTPUT 3
#define MAPLE_TEXT_QUIT 4
#define MAPLE_TEXT_WARNING 5
#define MAPLE_TEXT_ERROR 6
#define MAPLE_TEXT_STATUS 7
#define MAPLE_TEXT_PRETTY 8
#define MAPLE_TEXT_HELP 9
#define MAPLE_TEXT_DEBUG 10

/*
    This function is used to initialze the Maple kernel and return a 
    MKernelVector.  This function should not be called more than once.
    It should not be executed if Maple is already running (eg via an
    external function described by define_external() in Maple), since
    external functions provide their own kernel vector. 

    argc:  number of arguments in argv
    argv:  argument list (same arguments accepted by the command line
                           version of maple; see ?maple for complete list)
    cb:  I/O callbacks defined above.  All unused callbacks 
         must be set to 0.  If all callbacks are 0, stdout is used
         for output, and certain functions like readline() raise errors.
    user_data:  pointer to anything; gets passed through to the I/O callbacks
    void *info: version information -- always set to NULL
    errstr:  preallocated error buffer;  filled in if something goes wrong
             at startup.

*/
extern EXT_DECL MKernelVector M_DECL StartMaple( int argc, char *argv[],
                MCallBackVector cb, void *user_data, void *info, char *errstr );

extern EXT_DECL void M_DECL StopMaple( MKernelVector kv );

extern EXT_DECL M_BOOL M_DECL RestartMaple( MKernelVector kv, char *errstr );

extern EXT_DECL M_BOOL M_DECL RestartAndPauseMaple( MKernelVector kv, char *errstr );


/* ---------------------------------------------------------------------- */
/*                        Set/Query Library Path                          */
/* ---------------------------------------------------------------------- */

/* Returns an EXPSEQ of the current library path.  If the argument 'expseq'
   is non NULL, then the library path is set to the given value.  
   See ?libname for more detail.  

   To set libname to a single directory use:

       MapleLibName(kv,ToMapleString(kv,"/usr/local/lib"));

   To specify multiple paths, use:

       MapleLibName(kv,ToMapleExpressionSequence(kv,2,
	    ToMapleToMapleString(kv,"/usr/local/lib"),
	    ToMapleToMapleString(kv,"/home/mylib")));
*/
EXT_DECL ALGEB M_DECL MapleLibName( MKernelVector kv, ALGEB expseq );

/* ---------------------------------------------------------------------- */
/*                       Query/Set Kernel Options                         */
/* ---------------------------------------------------------------------- */

/* To query an option pass in val = NULL.  
   See ?kernelopts for a list of possible options and their possible values.
   Returns the value of the given setting prior to update (or the current
   value if val=NULL). 
*/
EXT_DECL ALGEB M_DECL MapleKernelOptions( MKernelVector kv, char *option, ALGEB val );


/* ---------------------------------------------------------------------- */
/*                       Assign to Maple Variables                        */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_DECL MapleAssign( MKernelVector kv, ALGEB lhs, ALGEB rhs );

EXT_DECL ALGEB M_DECL MapleAssignIndexed( MKernelVector kv, ALGEB lhs, 
			M_INT dim, M_INT *ind, ALGEB rhs );


/* ---------------------------------------------------------------------- */
/*             Error Handling and User Information                        */
/* ---------------------------------------------------------------------- */

/* Raising Errors */
EXT_DECL void M_DECL MapleRaiseError( MKernelVector kv, char *msg );
EXT_DECL void M_DECL MapleRaiseError1( MKernelVector kv, char *msg, ALGEB arg1 );
EXT_DECL void M_DECL MapleRaiseError2( MKernelVector kv, char *msg, ALGEB arg1, ALGEB arg2 );

/* Evaluate a C Function, Trapping Any Raised Errors */
EXT_DECL void* M_DECL MapleTrapError( MKernelVector kv, 
      void *(M_DECL *proc) ( void *data ), void *data, M_BOOL *errorflag );

/* Push an error handling routine on the error stack.  When an exception is raised, each
 * function on the error handling stack will be called.  The first argument is a
 * string containing the exception.  The second argument is the value that was passed
 * into MaplePushErrorProc when the function was pushed onto the stack.  This second
 * argument can be used to pass data into the error handling function.  
 *
 * As an example, imagine you have malloc'ed a piece of memory and want to make sure it
 * is freed, even if an execption occurs.  You could push an error handling
 * function that calls free on its second argument, and then pass the memory block as
 * the third argument to MaplePushErrorProc.  If an error does not occur, call
 * MaplePopErrorProc and then free the memory block directly */

/* push an error handing procedure and argument onto the error handling stack */
EXT_DECL void M_DECL MaplePushErrorProc( MKernelVector kv, void (M_DECL *errorproc)( const char *, void *), void *data );

/* Pop an error handling routine off the error stack */
EXT_DECL void M_DECL MaplePopErrorProc( MKernelVector kv );

/* Provide Run-time Information */
EXT_DECL void M_DECL MapleUserInfo( MKernelVector kv, int level, char *name, 
                         char *msg );   

/* output to Maple's output stream */
EXT_DECL int M_CDECL MaplePrintf( MKernelVector kv, const char *, ... );
EXT_DECL int M_CDECL MapleALGEB_Printf( MKernelVector kv, const char *, ... );
EXT_DECL ALGEB M_CDECL MapleALGEB_SPrintf( MKernelVector kv, const char *, ... );
EXT_DECL int M_DECL MapleALGEB_Printf0( MKernelVector kv, const char *format );
EXT_DECL int M_DECL MapleALGEB_Printf1( MKernelVector kv, const char *format,
     ALGEB a1 );
EXT_DECL int M_DECL MapleALGEB_Printf2( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2 );
EXT_DECL int M_DECL MapleALGEB_Printf3( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2, ALGEB a3 );
EXT_DECL int M_DECL MapleALGEB_Printf4( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2, ALGEB a3, ALGEB a4 );
EXT_DECL ALGEB M_DECL MapleALGEB_SPrintf0( MKernelVector kv, const char *format);
EXT_DECL ALGEB M_DECL MapleALGEB_SPrintf1( MKernelVector kv, const char *format,
     ALGEB a1 );
EXT_DECL ALGEB M_DECL MapleALGEB_SPrintf2( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2 );
EXT_DECL ALGEB M_DECL MapleALGEB_SPrintf3( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2, ALGEB a3 );
EXT_DECL ALGEB M_DECL MapleALGEB_SPrintf4( MKernelVector kv, const char *format,
     ALGEB a1, ALGEB a2, ALGEB a3, ALGEB a4 );

/* obsolete */
EXT_DECL void M_DECL InitMaplePrintf( MKernelVector kv );
EXT_DECL int M_DECL OldMaplePrintf( const char *, ... );

/* ---------------------------------------------------------------------- */
/*        Evaluate a Maple Procedure or DAG using hardware floats         */
/* ---------------------------------------------------------------------- */

EXT_DECL double M_DECL MapleEvalhf( MKernelVector kv, ALGEB s );

/* first argument here is args[1] */
EXT_DECL double M_DECL EvalhfMapleProc( MKernelVector kv, ALGEB fn, 
                            int nargs, double *args );
/* first argument here is args[1] */
EXT_DECL hfdata M_DECL EvalhfDataProc( MKernelVector kv, ALGEB fn, 
                            int nargs, hfdata *args );

/* ---------------------------------------------------------------------- */
/*                 Evaluate a Maple Procedure or statement                */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_CDECL EvalMapleProc( MKernelVector kv, ALGEB fn, int nargs, 
                /* ALGEB arg1, ALGEB arg2, */ ... );

EXT_DECL ALGEB M_DECL EvalMapleProcedure( MKernelVector kv, ALGEB fn, ALGEB args );

EXT_DECL ALGEB M_DECL EvalMapleStatement( MKernelVector kv, char *statement );

EXT_DECL ALGEB M_DECL MapleEval( MKernelVector kv, ALGEB s );

/* returns the value assigned to the given name 's' */
EXT_DECL ALGEB M_DECL MapleNameValue( MKernelVector kv, ALGEB s );

/* ---------------------------------------------------------------------- */
/*                          Arithmetic Functions                          */
/* ---------------------------------------------------------------------- */
EXT_DECL ALGEB M_DECL MapleNumericAdd( MKernelVector kv, ALGEB a, ALGEB b );
   
EXT_DECL ALGEB M_DECL MapleNumericSubtract( MKernelVector kv, ALGEB a, 
					    ALGEB b );

EXT_DECL ALGEB M_DECL MapleNumericMultiply( MKernelVector kv, ALGEB a, 
					    ALGEB b );

EXT_DECL ALGEB M_DECL MapleNumericPower( MKernelVector kv, ALGEB a, 
					 M_INT n );

EXT_DECL ALGEB M_DECL MapleIntegerDivide( MKernelVector kv, ALGEB a, 
					  ALGEB b, ALGEB *r );
    
EXT_DECL ALGEB M_DECL MapleFloatDivide( MKernelVector kv, ALGEB a, ALGEB b );

EXT_DECL M_INT M_DECL MapleIntegerAbsCompare( MKernelVector kv, ALGEB a, 
					      ALGEB b );

/* ---------------------------------------------------------------------- */
/*                          Data Queries                                  */
/* ---------------------------------------------------------------------- */

EXT_DECL M_BOOL M_DECL IsMapleAssignedName( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleComplexNumeric( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleComplex64( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleNumeric( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleFloat64( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleFunction( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleInteger( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleInteger8( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleInteger16( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleInteger32( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleInteger64( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleList( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleExpressionSequence( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleName( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleNULL( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMaplePointer( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMaplePointerNULL( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleProcedure( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleRTable( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleSet( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleStop( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleString( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleTable( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleTableBasedArray( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleUnassignedName( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL IsMapleUnnamedZero( MKernelVector kv, ALGEB s );

EXT_DECL MapleID M_DECL GetMapleID( MKernelVector kv, ALGEB dag );


/* ---------------------------------------------------------------------- */
/*                   Expression Sequence Manipulation                     */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_DECL NewMapleExpressionSequence( MKernelVector kv, int nargs);

EXT_DECL void M_DECL MapleExpseqAssign( MKernelVector kv, ALGEB expseq, 
                            M_INT i, ALGEB val );

EXT_DECL ALGEB M_DECL MapleExpseqSelect( MKernelVector kv, ALGEB expseq, 
                            M_INT i );

EXT_DECL ALGEB* M_DECL MapleGetElems( MKernelVector kv, ALGEB container, 
                M_INT *num_elems );

/* ---------------------------------------------------------------------- */
/*                             List Manipulation                          */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_DECL MapleListAlloc( MKernelVector kv, M_INT num_members );

EXT_DECL void M_DECL MapleListAssign( MKernelVector kv, ALGEB list, M_INT i, 
                            ALGEB val );

EXT_DECL ALGEB M_DECL MapleListSelect( MKernelVector kv, ALGEB list, M_INT i );


/* ---------------------------------------------------------------------- */
/*                     Conversion from Maple Objects                      */
/* ---------------------------------------------------------------------- */

EXT_DECL COMPLEXF32 M_DECL MapleToComplexFloat32( MKernelVector kv, ALGEB s ); 

EXT_DECL COMPLEXF64 M_DECL MapleToComplexFloat64( MKernelVector kv, ALGEB s ); 

EXT_DECL void M_DECL MapleToComplexSingle( MKernelVector kv, ALGEB s, COMPLEXF32 *c ); 

EXT_DECL void M_DECL MapleToComplexDouble( MKernelVector kv, ALGEB s, COMPLEXF64 *c );

EXT_DECL CXDAG M_DECL MapleToComplexFloatDAG( MKernelVector kv, ALGEB s );

EXT_DECL FLOAT32 M_DECL MapleToFloat32( MKernelVector kv, ALGEB s );

EXT_DECL FLOAT64 M_DECL MapleToFloat64( MKernelVector kv, ALGEB s );

EXT_DECL INTEGER8 M_DECL MapleToInteger8( MKernelVector kv, ALGEB s );

EXT_DECL INTEGER16 M_DECL MapleToInteger16( MKernelVector kv, ALGEB s );

EXT_DECL INTEGER32 M_DECL MapleToInteger32( MKernelVector kv, ALGEB s );

EXT_DECL INTEGER64 M_DECL MapleToInteger64( MKernelVector kv, ALGEB s );

EXT_DECL mpz_ptr M_DECL MapleToGMPInteger( MKernelVector kv, ALGEB s );

EXT_DECL M_BOOL M_DECL MapleToM_BOOL( MKernelVector kv, ALGEB s );

EXT_DECL M_INT M_DECL MapleToM_INT( MKernelVector kv, ALGEB s );

EXT_DECL void* M_DECL MapleToPointer( MKernelVector kv, ALGEB s );

EXT_DECL char* M_DECL MapleToString( MKernelVector kv, ALGEB s );

EXT_DECL void M_DECL MapleTohfData( MKernelVector kv, ALGEB s, hfdata *d );
EXT_DECL void M_DECL DoubleTohfData( MKernelVector kv, double re, double im, hfdata *d );

EXT_DECL double M_DECL RealhfData( MKernelVector kv, hfdata d );

/* Determine the number of arguments in a Maple object */
EXT_DECL M_INT M_DECL MapleNumArgs( MKernelVector kv, ALGEB expr );

/* ---------------------------------------------------------------------- */
/*       Rectangular Table (Vector, Matrix, Array) Manipulation           */
/* ---------------------------------------------------------------------- */

EXT_DECL void M_DECL RTableAppendAttribute( MKernelVector kv, RTableSettings *s, char *name );

EXT_DECL void M_DECL RTableAppendIndFn( MKernelVector kv, RTableSettings *s, 
                           ALGEB indfn );

EXT_DECL RTableData M_DECL RTableAssign( MKernelVector kv, ALGEB rt, 
                 M_INT *index, RTableData val );

EXT_DECL void M_DECL RTableAssignRef( MKernelVector kv, ALGEB rt, 
                 M_INT *index, RTableData *val );

EXT_DECL void M_DECL RTableBounds( MKernelVector kv, ALGEB rt, M_INT *bounds );

EXT_DECL M_INT M_DECL RTableLowerBound( MKernelVector kv, ALGEB rt, M_INT dim );

EXT_DECL M_INT M_DECL RTableUpperBound( MKernelVector kv, ALGEB rt, M_INT dim );

EXT_DECL ALGEB M_DECL RTableCopy( MKernelVector kv, RTableSettings *s, 
                            ALGEB rt );

EXT_DECL ALGEB M_DECL RTableCopyImPart( MKernelVector kv, RTableSettings *s, 
                       ALGEB rt );

EXT_DECL ALGEB M_DECL RTableCopyRealPart( MKernelVector kv, RTableSettings *s, 
                              ALGEB rt );

EXT_DECL ALGEB M_DECL RTableCreate( MKernelVector kv, RTableSettings *s, 
                         void *pdata, M_INT *bounds );

EXT_DECL ALGEB M_DECL RTableCreateFromDataBlock( MKernelVector kv, 
   RTableSettings *rts, M_INT *bounds, void *data_block, M_INT block_datatype );

EXT_DECL ALGEB M_DECL RTableCreateFromALGEB( MKernelVector kv, 
        RTableSettings *s, ALGEB initializer, M_INT *bounds );

EXT_DECL void* M_DECL RTableDataBlock( MKernelVector kv, ALGEB rt );

EXT_DECL void M_DECL RTableGetDefaults( MKernelVector kv, RTableSettings *s );

EXT_DECL void M_DECL RTableGetSettings( MKernelVector kv, RTableSettings *s, 
                           ALGEB rt );

EXT_DECL M_INT M_DECL RTableIndFn( MKernelVector kv, ALGEB rt, M_INT num );

EXT_DECL ALGEB M_DECL RTableIndFnArgs( MKernelVector kv, ALGEB rt, M_INT num );

EXT_DECL M_BOOL M_DECL RTableIsReal( MKernelVector kv, ALGEB rt );

EXT_DECL M_INT M_DECL RTableNumElements( MKernelVector kv, ALGEB rt );

EXT_DECL M_INT M_DECL RTableNumDimensions( MKernelVector kv, ALGEB rt );

EXT_DECL RTableData M_DECL RTableSelect( MKernelVector kv, ALGEB rt, 
                 M_INT *index );

EXT_DECL void M_DECL RTableSetAttribute( MKernelVector kv, RTableSettings *s, 
                 char *name );

EXT_DECL void M_DECL RTableSetIndFn( MKernelVector kv, RTableSettings *s, 
                ALGEB indfn );

EXT_DECL void M_DECL RTableSetType( MKernelVector kv, RTableSettings *s, 
                 M_INT id, char *name );

EXT_DECL void M_DECL RTableSparseCompact( MKernelVector kv, ALGEB rt );

EXT_DECL NAG_INT* M_DECL RTableSparseIndexRow( MKernelVector kv, ALGEB rt, 
                M_INT dim );

EXT_DECL ALGEB M_DECL RTableSparseIndexSort( MKernelVector kv, ALGEB rt, 
                M_INT by_dim );

EXT_DECL void M_DECL RTableSparseSetNumElems( MKernelVector kv, ALGEB rt, 
                M_INT num ); 

EXT_DECL M_INT M_DECL RTableSparseSize( MKernelVector kv, ALGEB rt );

EXT_DECL void M_DECL RTableSparseResize( MKernelVector kv, ALGEB rt, M_INT size );

EXT_DECL ALGEB M_DECL RTableZipReIm( MKernelVector kv, RTableSettings *s, 
                            ALGEB rt_re, ALGEB rt_im );

/* Wrapper Helper Functions */
/* creates a foreign maple array referencing 'ptr' */
EXT_DECL void* M_CDECL ArrayToMaple( MKernelVector kv, void *ptr, 
	    M_INT data_type, M_INT storage, M_INT p1, M_INT p2, 
	    M_INT order, ALGEB indfn, M_INT subtype,
	    ALGEB *ref, M_INT num_dims, M_INT offset, 
	    /* low-bound 1, upp-bound 1, low-bound 2, upp-bound 2 */ ... );

/* creates a C array from 'rt', or just points to it's data if params match */
EXT_DECL void* M_CDECL MapleToArray( MKernelVector kv, ALGEB rt,
	    M_INT data_type, M_INT storage, M_INT p1, M_INT p2, 
	    M_INT order, ALGEB indfn, M_INT subtype,
	    ALGEB *ref, M_INT num_dims, M_INT offset, 
	    /* low-bound 1, upp-bound 1, low-bound 2, upp-bound 2 */ ... ); 

/* allocate a new array and return it's data pointer */
EXT_DECL void* M_CDECL MapleAllocArray( MKernelVector kv, ALGEB rt,
	    M_INT data_type, M_INT storage, M_INT p1, M_INT p2, 
	    M_INT order, ALGEB indfn, M_INT subtype,
	    ALGEB *ref, M_INT num_dims, M_INT offset, 
	    /* low-bound 1, upp-bound 1, low-bound 2, upp-bound 2 */ ... ); 

/* Data Selection */
EXT_DECL ALGEB M_DECL MapleSelectImaginaryPart( MKernelVector kv, ALGEB s );

EXT_DECL ALGEB M_DECL MapleSelectIndexed( MKernelVector kv, ALGEB s, M_INT dim, M_INT *ind );

EXT_DECL ALGEB M_DECL MapleSelectRealPart( MKernelVector kv, ALGEB s );

/* Data Uniquification */
EXT_DECL ALGEB M_DECL MapleUnique( MKernelVector kv,ALGEB s );

/* ---------------------------------------------------------------------- */
/*                             Table Manipulation                         */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_DECL MapleTableAlloc( MKernelVector kv );

EXT_DECL void M_DECL MapleTableAssign( MKernelVector kv, ALGEB table, ALGEB ind, ALGEB val );

EXT_DECL void M_DECL MapleTableDelete( MKernelVector kv, ALGEB table, ALGEB ind );

EXT_DECL M_BOOL M_DECL MapleTableHasEntry( MKernelVector kv, ALGEB table, ALGEB ind );

EXT_DECL ALGEB M_DECL MapleTableSelect( MKernelVector kv, ALGEB table, ALGEB ind );

/* ---------------------------------------------------------------------- */
/*                        Conversion to Maple Objects                     */
/* ---------------------------------------------------------------------- */

EXT_DECL ALGEB M_DECL ToMapleBoolean( MKernelVector kv, M_INT b );

EXT_DECL ALGEB M_DECL ToMapleChar( MKernelVector kv, M_INT c );

EXT_DECL ALGEB M_DECL ToMapleComplex( MKernelVector kv, double re, double im );

EXT_DECL ALGEB M_DECL ToMapleComplexFloat( MKernelVector kv, ALGEB re, ALGEB im );

EXT_DECL ALGEB M_CDECL ToMapleExpressionSequence( MKernelVector kv, int nargs, /* ALGEB arg1, ALGEB arg2, */ ... );

EXT_DECL ALGEB M_DECL ToMapleInteger( MKernelVector kv, M_INT i ); 

EXT_DECL ALGEB M_DECL ToMapleInteger64( MKernelVector kv, INTEGER64 i ); 

EXT_DECL ALGEB M_DECL GMPIntegerToMaple( MKernelVector kv, mpz_ptr g );

EXT_DECL ALGEB M_DECL ToMapleFloat( MKernelVector kv, double f ); 

EXT_DECL ALGEB M_DECL ToMapleHFloat( MKernelVector kv, double f ); 

EXT_DECL ALGEB M_CDECL ToMapleFunction( MKernelVector kv, ALGEB fn, int nargs, ... );
EXT_DECL ALGEB M_DECL ToMapleFunc( MKernelVector kv, ALGEB fn, ALGEB expseq );

EXT_DECL ALGEB M_DECL ToMapleName( MKernelVector kv, const char *n, M_BOOL is_global ); 

EXT_DECL ALGEB M_DECL ToMapleNULL( MKernelVector kv ); 

EXT_DECL ALGEB M_DECL ToMapleNULLPointer( MKernelVector kv ); 

EXT_DECL ALGEB M_DECL ToMaplePointer( MKernelVector kv, void *v, M_INT type );

EXT_DECL ALGEB M_DECL ToMapleRelation( MKernelVector kv, const char *rel, 
                          ALGEB lhs, ALGEB rhs );

EXT_DECL ALGEB M_DECL ToMapleString( MKernelVector kv, const char *s ); 

EXT_DECL ALGEB M_DECL ToMapleUneval( MKernelVector kv, ALGEB s ); 

EXT_DECL ALGEB M_DECL ToMaplehfData( MKernelVector kv, hfdata *d ); 

/* ---------------------------------------------------------------------- */
/*                       Memory Management                                */
/* ---------------------------------------------------------------------- */

/* Allocate Using Maple's Allocator */
EXT_DECL void* M_DECL MapleAlloc( MKernelVector kv, M_INT nbytes );

/* Allocate a new Maple Object */
EXT_DECL ALGEB M_DECL MapleNew( MKernelVector kv, MapleID id, M_INT len );
EXT_DECL ALGEB M_DECL MapleCreate( MKernelVector kv, MapleID id, M_INT len,
    ALGEB *elems );

/* Free Memory Allocated By Maple's Allocator */ 
EXT_DECL void M_DECL MapleDispose( MKernelVector kv, ALGEB s );

/* Allow a to be garbage collected */
EXT_DECL void M_DECL MapleGcAllow( MKernelVector kv, ALGEB a );

/* Prevent a from being garbage collected */
EXT_DECL void M_DECL MapleGcProtect( MKernelVector kv, ALGEB a );

/* Check if an object is protected from being garbage collected */
EXT_DECL M_BOOL M_DECL MapleGcIsProtected( MKernelVector kv, ALGEB a );

/* Apply to all stored Maple objects during a gc sweep */
EXT_DECL void M_DECL MapleGcMark( MKernelVector kv, ALGEB a );

/* ---------------------------------------------------------------------- */
/*                 GMP Memory Management Functions                        */
/* ---------------------------------------------------------------------- */

/*
 * Maple sets the gmp memory allocation functions, if one wants to use GMP in
 * external call with different memory management functions (which is proabably
 * a good idea) then use the functions below to set and unset the memory functions.  
 * Don't use the mp_set_memory_functions.
 *
 * These functions maintain a stack of allocations functions.  Setting functions this
 * way is thread safe.  i.e. you are setting memory functions for the current thread
 * only.
 */

/* push a new set of memory functions on the stack */
EXT_DECL void M_DECL MaplePushGMPAllocators( MKernelVector kv,
                                void *(GMP_DECL *malloc)( size_t ),
                                void *(GMP_DECL *realloc)( void *, size_t, size_t ),
                                void (GMP_DECL *free)( void *, size_t ) );

/* pop the current set of memory functions off the stack */
EXT_DECL void M_DECL MaplePopGMPAllocators( MKernelVector kv );


/* ---------------------------------------------------------------------- */
/*                Foreign Object (MaplePointer) Management                */
/* ---------------------------------------------------------------------- */

/* query the user-supplied pointer type marker */
EXT_DECL M_INT M_DECL MaplePointerType( MKernelVector kv, ALGEB a );

/* set the pointer type marker */
EXT_DECL void M_DECL MaplePointerSetType( MKernelVector kv, ALGEB a, M_INT type );

/* Set the function to be called during a gc mark sweep.
   MapleGcMark can then be called on all Maple objects contained
   in the foreign data-structure so they won't be collected. 
*/
EXT_DECL void M_DECL MaplePointerSetMarkFunction( MKernelVector kv, ALGEB a,
    void (M_DECL *markfn)( ALGEB a ) );

/* Set a function to call when a Pointer object is about to be
   garbage collected.
*/
EXT_DECL void M_DECL MaplePointerSetDisposeFunction( MKernelVector kv, ALGEB a,
    void (M_DECL *disposefn)( ALGEB a ) );

/* set the function to be called in order to convert a Pointer object
   into a printable Maple object during printing */
EXT_DECL void M_DECL MaplePointerSetPrintFunction( MKernelVector kv, ALGEB a,
    ALGEB (M_DECL *printfn)( ALGEB a ) );

/* ---------------------------------------------------------------------- */
/*                Extra Call Backs                                        */
/* ---------------------------------------------------------------------- */

/* set a function to be called after a restart */
EXT_DECL void M_DECL RegisterRestartCallBack( MKernelVector kv,
    void (M_DECL *restartCB)( void *data ) );


/* ------------------------------------------------------------------------ */
/*                Access to the Random Number Generator                     */
/* ------------------------------------------------------------------------ */

EXT_DECL INTEGER32 M_DECL MapleRandomInt32( MKernelVector kv );
EXT_DECL INTEGER64 M_DECL MapleRandomInt64( MKernelVector kv );
EXT_DECL M_INT M_DECL MapleRandomM_INT( MKernelVector kv );
EXT_DECL double M_DECL MapleRandomDouble01( MKernelVector kv );
EXT_DECL ALGEB M_DECL MapleRandomSoftwareFloat01( MKernelVector kv );
EXT_DECL ALGEB M_DECL MapleRandomSoftwareInteger( MKernelVector kv, M_INT bits );
EXT_DECL M_INT M_DECL MapleRandomCalcBits( MKernelVector kv, ALGEB range );
EXT_DECL INTEGER32 M_DECL MapleRandomRangeInt32( MKernelVector kv, INTEGER32 range, M_INT bits );
EXT_DECL ALGEB M_DECL MapleRandomRangeSoftwareInt( MKernelVector kv, ALGEB range, M_INT bits );


/* ---------------------------------------------------------------------- */
/*                       Access the Help System                           */
/* ---------------------------------------------------------------------- */

/* Possible values of the attribute passed to the writeAttrib call back. */
#define FN_NORM  0 	/* normal text mode  */
#define FN_ITAL  1	/* italic text mode */
#define FN_BOLD  3	/* boldfaced text mode */
#define FN_UNDER 4	/* underlined text mode */

/* 
    The HelpLookUpText function searches for and retrieves a help page or 
    a section of a help page, based on the topic passed to it. The results 
    are passed as a stream of characters and attributes to the specified 
    call-back functions. 

    topic:  Specifies the help page retrieved

    section:  Indicates which section of the page to display. If this
	     is passed as "" or NULL, the entire page is displayed. 
             To restrict display to a particular section of the page, 
             one of the following values can be passed: 

	      "usage" 
		    Shows just the function name (one-line description) and 
		    calling sequence information. 
	      "description" 
		    Shows the detailed description of the function. 
	      "examples" 
		    Shows examples of the function's usage. 
	      "seealso" 
		    Shows a list of alternate topics that may be related to 
		    this function. 

    writeChar: Function to which output is sent.  The writeChar
	       function can terminate rendering by returning TRUE. 

    writeAttrib: Function to which attribute information is passed.
	         Each given attribute applies to all subsequent characters
                 sent to writeChar until a new attribute is given.
                 Possible attribute values are describe above (FN_*). 
		 The writeAttrib function can be omitted by passing NULL 
		 for the writeAttrib parameter. 

    width:  Indicates the width, in characters, to which the help
	    information should be formatted. 

    data: The data parameter given to StartMaple.

    MapleHelp returns NULL if successful, or it returns a pointer 
    to an error message if unsuccessful.

*/

EXT_DECL char * M_DECL MapleHelp( 
    MKernelVector kv, 
    char *topic, 
    char *section,
    M_BOOL (M_DECL *writechar) ( void *data, int c ),
    M_BOOL (M_DECL *writeattrib) ( void *data, int a ),
    int width, 
    void *data );

/* ------------------------------------------------------------------------ */
/*                           Misc Functions                                 */
/* ------------------------------------------------------------------------ */

/* Periodically call in safe spot to allow user interrupt via non-trappable
 * Maple exception.  Will not return if the user has pressed the Stop
 * button or hit Ctrl-C prior to call.
 */
EXT_DECL void M_DECL MapleCheckInterrupt( MKernelVector kv );

/*
 * Check to see if an interrupt has occured without actually jumping away.  This gives
 * the programmer a chance to clean up before actually calling MapleCheckInterrupt
 * to handle the interrupt.
 */
EXT_DECL M_BOOL M_DECL MapleGetInterruptValue( MKernelVector kv );

/* Override default Maple version to launch in Windows OpenMaple.
 * Unix and Linux look at the MAPLE environment variable to find
 * the installed version of Maple, and won't respect a call to SetMapleVersion.
 *
 * Call prior to StartMaple() with arguments like "10" or "11". 
 * Only use the major release number (ie. from a version that has its own
 * directory under "C:\Program Files").  
 */
EXT_DECL void M_DECL SetMapleVersion( char *version );

/* ------------------------------------------------------------------------ */
/*                           Internal Use Only                              */
/* ------------------------------------------------------------------------ */

EXT_DECL void M_DECL RegisterPlotCallBack( MKernelVector kv,
    void (M_DECL *callback) ( void *data, void* ) );
EXT_DECL void M_DECL RegisterStreamCallBack( MKernelVector kv, char *name, 
    ALGEB (M_DECL *callback) ( void *data, ALGEB ) );
EXT_DECL void M_DECL RegisterGUIStreamCallBack( MKernelVector kv, char *name, 
    ALGEB (M_DECL *callback) ( void *data, ALGEB ) );
EXT_DECL void M_DECL RegisterPostRestartCallBack( MKernelVector kv, 
    void (M_DECL *callback) ( void *data ) );
EXT_DECL char* M_DECL MapleAppendFeature( MKernelVector kv, char*, char*, char*, char* );

/* ------------------------------------------------------------------------ */
/*                         Multi-threaded support                           */
/* ------------------------------------------------------------------------ */

EXT_DECL ALGEB M_DECL MapleMutexCreate( MKernelVector kv, ALGEB options );
EXT_DECL void M_DECL MapleMutexDestroy( MKernelVector kv, ALGEB id );
EXT_DECL void M_DECL MapleMutexLock( MKernelVector kv, ALGEB id );
EXT_DECL void M_DECL MapleMutexUnlock( MKernelVector kv, ALGEB id );

/* ------------------------------------------------------------------------ */
/*                         Task Programming support                         */
/* ------------------------------------------------------------------------ */

#define MAPLE_ROOT_TASK 0
EXT_DECL int M_DECL MapleStartRootTask( MKernelVector kv, void *root, 
        int (*TaskFunction)( void *parent, int arg_number, void *self ), 
        void *args, void (*MarkTaskFunction)( void *self ), void **retValue, 
        M_INT options );
EXT_DECL void M_DECL MapleCreateContinuationTask( MKernelVector kv, 
        int (*TaskFunction)( void *parent, int arg_number, void *self ), 
        void *args, void (*MarkTaskFunction)( void *self ) );
EXT_DECL void M_DECL MapleStartChildTask( MKernelVector kv, int arg_number, 
        int (*TaskFunction)( void *parent, int arg_number, void *self ), 
        void *args, void (*MarkTaskFunction)( void *self ) );
EXT_DECL int M_DECL MapleTaskReturn( MKernelVector kv, void *value,
        void (*MarkFunction)( void *self ) );

EXT_DECL M_INT M_DECL MapleRegisterThread( MKernelVector kv, M_INT options );
EXT_DECL M_INT M_DECL MapleUnregisterThread( MKernelVector kv );

#ifdef __cplusplus
}
#endif

#endif /* _MAPLE_EXTERNAL_C_H_ */
