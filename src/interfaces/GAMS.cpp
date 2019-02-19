/* HiGHS link code */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"
#include "optcc.h"

/* HiGHS API */
#include "Highs.h"

#if defined(_WIN32)
#if !defined(STDCALL)
#define STDCALL __stdcall
#endif
#if !defined(DllExport)
#define DllExport __declspec(dllexport)
#endif
#else
#if !defined(STDCALL)
#define STDCALL
#endif
#if !defined(DllExport)
#define DllExport
#endif
#endif

struct gamshighs_s
{
   gmoHandle_t gmo;
   gevHandle_t gev;
   optHandle_t opt;
   int         debug;

   Highs*      highs;
   HighsLp*    lp;
};
typedef struct gamshighs_s gamshighs_t;

static
int setupProblem(
   gamshighs_t* gh
)
{
   int numCol;
   int numRow;
   int numNz;
   int i;
   int rc = 1;

   assert(gh != NULL);
   assert(gh->highs == NULL);
   assert(gh->lp == NULL);

   HighsOptions options;
   gh->highs = new Highs(options);

   numCol = gmoN(gh->gmo);
   numRow = gmoM(gh->gmo);
   numNz = gmoNZ(gh->gmo);

   gh->lp = new HighsLp();

   gh->lp->numRow_ = numRow;
   gh->lp->numCol_ = numCol;
   gh->lp->nnz_ = numNz;

   /* columns */
   gh->lp->colUpper_.resize(numCol);
   gh->lp->colLower_.resize(numCol);
   gmoGetVarLower(gh->gmo, &gh->lp->colLower_[0]);
   gmoGetVarUpper(gh->gmo, &gh->lp->colUpper_[0]);

   /* objective */
   gh->lp->colCost_.resize(numCol);
   gmoGetObjVector(gh->gmo, &gh->lp->colCost_[0], NULL);
   if( gmoSense(gh->gmo) == gmoObj_Min )
      gh->lp->sense_ = 1;
   else
      gh->lp->sense_ = -1;

   /* row left- and right-hand-side */
   gh->lp->rowLower_.resize(numRow);
   gh->lp->rowUpper_.resize(numRow);
   for( i = 0; i < numRow; ++i )
   {
      switch( gmoGetEquTypeOne(gh->gmo, i) )
      {
         case gmoequ_E :
            gh->lp->rowLower_[i] = gh->lp->rowUpper_[i] = gmoGetRhsOne(gh->gmo, i);
            break;

         case gmoequ_G :
            gh->lp->rowLower_[i] = gmoGetRhsOne(gh->gmo, i);
            gh->lp->rowUpper_[i] = HIGHS_CONST_INF;
            break;

         case gmoequ_L :
            gh->lp->rowLower_[i] = -HIGHS_CONST_INF;
            gh->lp->rowUpper_[i] = gmoGetRhsOne(gh->gmo, i);
            break;

         case gmoequ_N:
         case gmoequ_X:
         case gmoequ_C:
         case gmoequ_B:
            /* these should not occur */
            goto TERMINATE;
      }
   }

   /* coefficients matrix */
   gh->lp->Astart_.resize(numCol + 1);
   gh->lp->Aindex_.resize(numNz);
   gh->lp->Avalue_.resize(numNz);
   gmoGetMatrixCol(gh->gmo, &gh->lp->Astart_[0], &gh->lp->Aindex_[0], &gh->lp->Avalue_[0], NULL);

   gh->highs->initializeLp(*gh->lp);

   rc = 0;
TERMINATE:

   return rc;
}

static
int processSolve(
   gamshighs_t* gh
)
{
   double* x;
   double* pi;

   assert(gh != NULL);
   assert(gh->highs != NULL);
   assert(gh->lp != NULL);

   gmoSetHeadnTail(gh->gmo, gmoHresused, gevTimeDiffStart(gh->gev));
   gmoSetHeadnTail(gh->gmo, gmoHiterused, 42 /* FIXME number of iterations */);

   /* FIXME below assumes best possible outcome: solved to optimality, have optimal primal and dual solution */

   x = new double[gmoN(gh->gmo)];
   pi = new double[gmoM(gh->gmo)];

   /* TODO fill x with primal solution
    * TODO fill pi with dual solution (w.r.t. rows)
    * TODO probably should use gmoSetSolution or gmoSetSolution8
    */

   gmoSetSolution2(gh->gmo, x, pi);
   gmoCompleteSolution(gh->gmo);

   gmoModelStatSet(gh->gmo, gmoModelStat_OptimalGlobal);
   gmoSolveStatSet(gh->gmo, gmoSolveStat_Normal);

   delete[] x;
   delete[] pi;

   return 0;
}


extern "C"
{


void his_Initialize(void)
{
   gmoInitMutexes();
   gevInitMutexes();
   optInitMutexes();
}

void his_Finalize(void)
{
   gmoFiniMutexes();
   gevFiniMutexes();
   optFiniMutexes();
}

DllExport void STDCALL hisXCreate(void** Cptr)
{
   assert(Cptr != NULL);

   *Cptr = calloc(1, sizeof(gamshighs_t));
}

DllExport int STDCALL hiscreate(void** Cptr, char* msgBuf, int msgBufLen)
{
   assert(Cptr != NULL);
   assert(msgBufLen > 0);
   assert(msgBuf != NULL);

   *Cptr = calloc(1, sizeof(gamshighs_t));

   msgBuf[0] = 0;

   return 1;
}

DllExport void STDCALL hisXFree(void** Cptr)
{
   assert(Cptr != NULL);
   assert(*Cptr != NULL);

   free(*Cptr);
   *Cptr = NULL;

   gmoLibraryUnload();
   gevLibraryUnload();
}

DllExport int STDCALL hisfree(void** Cptr)
{
   hisXFree(Cptr);

   return 1;
}

/* comp returns the compatibility mode:
           0: client is too old for the DLL, no compatibility
           1: client version and DLL version are the same, full compatibility
           2: client is older than DLL, but defined as compatible, backward compatibility
           3: client is newer than DLL, forward compatibility
           FIXME: for now, we just claim full compatibility
 */
DllExport int STDCALL C__hisXAPIVersion(int api, char* Msg, int* comp)
{
   *comp = 1;
   return 1;
}

DllExport int STDCALL D__hisXAPIVersion(int api, char* Msg, int* comp)
{
   *comp = 1;
   return 1;
}

DllExport int STDCALL C__hisXCheck(const char* funcn, int ClNrArg, int Clsign[], char* Msg) { return 1; }

DllExport int STDCALL D__hisXCheck(const char* funcn, int ClNrArg, int Clsign[], char* Msg) { return 1; }

DllExport int STDCALL C__hisReadyAPI(void* Cptr, gmoHandle_t Gptr, optHandle_t Optr)
{
   gamshighs_t* gh;

   assert(Cptr != NULL);
   assert(Gptr != NULL);

   char msg[256];
   if(!gmoGetReady(msg, sizeof(msg)))
      return 1;
   if(!gevGetReady(msg, sizeof(msg)))
      return 1;

   gh = (gamshighs_t*)Cptr;
   gh->gmo = Gptr;
   gh->gev = (gevHandle_t)gmoEnvironment(gh->gmo);
   gh->opt = Optr;

   return 0;
}

DllExport int STDCALL C__hisCallSolver(void* Cptr)
{
   int rc = 1;
   char buffer[1024];
   gamshighs_t* gh;

   gh = (gamshighs_t*)Cptr;
   assert(gh->gmo != NULL);
   assert(gh->gev != NULL);

   gevLogStat(gh->gev, "This is the GAMS link to HiGHS.");

   gmoModelStatSet(gh->gmo, gmoModelStat_NoSolutionReturned);
   gmoSolveStatSet(gh->gmo, gmoSolveStat_SystemErr);

   /*
   if( dooptions(highs) )
      goto TERMINATE;
    */


   /* get the problem into a normal form */
   gmoObjStyleSet(gh->gmo, gmoObjType_Fun);
   gmoObjReformSet(gh->gmo, 1);
   gmoIndexBaseSet(gh->gmo, 0);
   gmoSetNRowPerm(gh->gmo); /* hide =N= rows */
   gmoMinfSet(gh->gmo, -HIGHS_CONST_INF);
   gmoPinfSet(gh->gmo,  HIGHS_CONST_INF);

   if( !setupProblem(gh) )
      goto TERMINATE;

   /* set timelimit */
//   gh->model->dblOption[DBLOPT_TIME_LIMIT] = gevGetDblOpt(gh->gev, gevResLim);

   gevTimeSetStart(gh->gev);

   /* TODO solve the problem here */

   /* pass solution, status, etc back to GMO here */

   /* process solve outcome */
   if( !processSolve(gh) )
      goto TERMINATE;

   rc = 0;
TERMINATE:

   if( gh->opt != NULL )
      optFree(&gh->opt);

   delete gh->lp;
   gh->lp = NULL;

   delete gh->highs;
   gh->highs= NULL;

   return rc;
}

DllExport int STDCALL C__hisHaveModifyProblem(void* Cptr) { return 0; }

DllExport int STDCALL C__hisModifyProblem(void* Cptr)
{
   assert(Cptr != NULL);
   return 1;
}

void his_Initialize(void);
void his_Finalize(void);

} // extern "C"
