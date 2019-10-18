/* HiGHS link code */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

/* GAMS API */
#include "gmomcc.h"
#include "gevmcc.h"

typedef struct optRec* optHandle_t;

/* HiGHS API */
#include "Highs.h"
#include "io/LoadOptions.h"  /* for loadOptionsFromFile */
#include "io/FilereaderLp.h"
#include "io/FilereaderMps.h"

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
   int         debug;

   Highs*      highs;
   HighsLp*    lp;
   HighsOptions* options;
};
typedef struct gamshighs_s gamshighs_t;

static
void gevprint(
   int           level,
   const char*   msg,
   void*         msgcb_data)
{
   gevHandle_t gev = (gevHandle_t) msgcb_data;
   gevLogPChar(gev, msg);
}

static
void gevlog(
   HighsMessageType type,
   const char*      msg,
   void*            msgcb_data)
{
   gevHandle_t gev = (gevHandle_t) msgcb_data;
   if( type == HighsMessageType::INFO )
      gevLogPChar(gev, msg);
   else
      gevLogStatPChar(gev, msg);
}

static
int setupOptions(
   gamshighs_t* gh
)
{
   assert(gh != NULL);
   assert(gh->options == NULL);

   gh->options = new HighsOptions;

   gh->options->time_limit = gevGetDblOpt(gh->gev, gevResLim);
   gh->options->simplex_iteration_limit = gevGetIntOpt(gh->gev, gevIterLim);

   if( gevGetIntOpt(gh->gev, gevUseCutOff) )
      gh->options->dual_objective_value_upper_bound = gevGetDblOpt(gh->gev, gevCutOff);

   if( gmoOptFile(gh->gmo) > 0 )
   {
      char optfilename[GMS_SSSIZE];
      gmoNameOptFile(gh->gmo, optfilename);
      gh->options->options_file = optfilename;
      if( !loadOptionsFromFile(*gh->options) )
         return 1;
   }

   gh->options->printmsgcb = gevprint;
   gh->options->logmsgcb = gevlog;
   gh->options->msgcb_data = (void*)gh->gev;

   return 0;
}

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
   assert(gh->options != NULL);
   assert(gh->highs == NULL);
   assert(gh->lp == NULL);

   gh->highs = new Highs(*gh->options);

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

   //FilereaderLp().writeModelToFile("highs.lp", *gh->lp);
   //FilereaderMps().writeModelToFile("highs.mps", *gh->lp);

   rc = 0;
TERMINATE:

   return rc;
}

static
int processSolve(
   gamshighs_t* gh
)
{
   assert(gh != NULL);
   assert(gh->highs != NULL);
   assert(gh->lp != NULL);

   Highs* highs = gh->highs;

   gmoSetHeadnTail(gh->gmo, gmoHresused, gevTimeDiffStart(gh->gev));
   gmoSetHeadnTail(gh->gmo, gmoHiterused, highs->getIterationCount());

   HighsSolution sol = gh->highs->getSolution();
   bool havesol = sol.col_value.size() > 0;

   if( havesol )
   {
      /* TODO probably should use gmoSetSolution or gmoSetSolution8 */
      gmoSetSolution2(gh->gmo, &sol.col_value[0], &sol.row_dual[0]);
      gmoCompleteSolution(gh->gmo);
   }

   gmoSolveStatSet(gh->gmo, gmoSolveStat_Normal);

   switch( gh->highs->getModelStatus() )
   {
      case HighsModelStatus::NOTSET:
      case HighsModelStatus::LOAD_ERROR:
      case HighsModelStatus::MODEL_ERROR:
      case HighsModelStatus::MODEL_EMPTY:
      case HighsModelStatus::PRESOLVE_ERROR:
      case HighsModelStatus::SOLVE_ERROR:
      case HighsModelStatus::POSTSOLVE_ERROR:
         gmoModelStatSet(gh->gmo, gmoModelStat_ErrorNoSolution);
         gmoSolveStatSet(gh->gmo, gmoSolveStat_SolverErr);
         break;

      case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
         /* TODO is solution feasible? */
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_InfeasibleIntermed : gmoModelStat_NoSolutionReturned);
         gmoSolveStatSet(gh->gmo, gmoSolveStat_Solver);
         break;

      case HighsModelStatus::PRIMAL_UNBOUNDED:
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_Unbounded : gmoModelStat_UnboundedNoSolution);
         break;

      case HighsModelStatus::PRIMAL_INFEASIBLE:
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_InfeasibleGlobal : gmoModelStat_InfeasibleNoSolution);
         break;

      case HighsModelStatus::PRIMAL_FEASIBLE:
         assert(havesol);
         gmoModelStatSet(gh->gmo, gmoModelStat_Feasible);
         break;

      case HighsModelStatus::DUAL_FEASIBLE:
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_InfeasibleIntermed : gmoModelStat_NoSolutionReturned);
         gmoSolveStatSet(gh->gmo, gmoSolveStat_Solver);
         break;

      case HighsModelStatus::OPTIMAL:
         assert(havesol);
         gmoModelStatSet(gh->gmo, gmoModelStat_OptimalGlobal);
         break;

      case HighsModelStatus::REACHED_TIME_LIMIT:
         /* TODO is solution feasible? */
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_InfeasibleIntermed : gmoModelStat_NoSolutionReturned);
         gmoSolveStatSet(gh->gmo, gmoSolveStat_Resource);
         break;

      case HighsModelStatus::REACHED_ITERATION_LIMIT:
         /* TODO is solution feasible? */
         gmoModelStatSet(gh->gmo, havesol ? gmoModelStat_InfeasibleIntermed : gmoModelStat_NoSolutionReturned);
         gmoSolveStatSet(gh->gmo, gmoSolveStat_Iteration);
         break;
   }

   return 0;
}


extern "C"
{


void his_Initialize(void)
{
   gmoInitMutexes();
   gevInitMutexes();
}

void his_Finalize(void)
{
   gmoFiniMutexes();
   gevFiniMutexes();
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
   assert(Optr == NULL);

   char msg[256];
   if(!gmoGetReady(msg, sizeof(msg)))
      return 1;
   if(!gevGetReady(msg, sizeof(msg)))
      return 1;

   gh = (gamshighs_t*)Cptr;
   gh->gmo = Gptr;
   gh->gev = (gevHandle_t)gmoEnvironment(gh->gmo);

   return 0;
}

#define XQUOTE(x) QUOTE(x)
#define QUOTE(x) #x

DllExport int STDCALL C__hisCallSolver(void* Cptr)
{
   int rc = 1;
   gamshighs_t* gh;
   HighsStatus status;

   gh = (gamshighs_t*)Cptr;
   assert(gh->gmo != NULL);
   assert(gh->gev != NULL);

   gevLogStatPChar(gh->gev, "HiGHS " XQUOTE(HIGHS_VERSION_MAJOR) "." XQUOTE(HIGHS_VERSION_MINOR) "." XQUOTE(HIGHS_VERSION_PATCH) " [date: " HIGHS_COMPILATION_DATE ", git hash: " HIGHS_GITHASH "]\n");
   gevLogStatPChar(gh->gev, "Copyright (c) 2019 ERGO-Code under MIT licence terms.\n");

   gmoModelStatSet(gh->gmo, gmoModelStat_NoSolutionReturned);
   gmoSolveStatSet(gh->gmo, gmoSolveStat_SystemErr);

   /* get the problem into a normal form */
   gmoObjStyleSet(gh->gmo, gmoObjType_Fun);
   gmoObjReformSet(gh->gmo, 1);
   gmoIndexBaseSet(gh->gmo, 0);
   gmoSetNRowPerm(gh->gmo); /* hide =N= rows */
   gmoMinfSet(gh->gmo, -HIGHS_CONST_INF);
   gmoPinfSet(gh->gmo,  HIGHS_CONST_INF);

   if( setupOptions(gh) )
      goto TERMINATE;

   if( setupProblem(gh) )
      goto TERMINATE;

   gevTimeSetStart(gh->gev);

   /* solve the problem */
   status = gh->highs->run();
   if( status != HighsStatus::OK )
      goto TERMINATE;

   /* pass solution, status, etc back to GMO */
   if( processSolve(gh) )
      goto TERMINATE;

   rc = 0;
TERMINATE:

   delete gh->lp;
   gh->lp = NULL;

   delete gh->highs;
   gh->highs= NULL;

   delete gh->options;
   gh->options = NULL;

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
