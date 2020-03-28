# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, allocator, make_unique
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdio cimport FILE

cdef extern from "HConst.h" nogil:
    cdef enum HighsPrintMessageLevel:
        ML_MIN = 0
        ML_NONE = ML_MIN
        ML_VERBOSE = 1
        ML_DETAILED = 2
        ML_MINIMAL = 4
        ML_ALWAYS = ML_VERBOSE | ML_DETAILED | ML_MINIMAL
        ML_MAX = ML_ALWAYS

cdef extern from "HighsOptions.h" nogil:
    cdef cppclass HighsOptions:
        FILE * output
        int message_level
        bool mip
        string solution_file
        bool write_solution_to_file
        bool write_solution_pretty

cdef extern from "HighsRuntimeOptions.h" nogil:
    bool loadOptions(int argc, char** argv, HighsOptions& options)

cdef extern from "HighsIO.h" nogil:
    void HighsPrintMessage(FILE* pass_output, const int message_level, const int level, const char* format, ...)

cdef extern from "HighsLp.h" nogil:
    cdef cppclass HighsLp:
        int numCol_
        int numRow_
        int numInt_
        vector[double] Avalue_
        string model_name_

    cdef enum HighsModelStatus:
        HighsModelStatusNOTSET "HighsModelStatus::NOTSET"
        HighsModelStatusLOAD_ERROR "HighsModelStatus::LOAD_ERROR"
        HighsModelStatusMODEL_ERROR "HighsModelStatus::MODEL_ERROR"
        HighsModelStatusMODEL_EMPTY "HighsModelStatus::MODEL_EMPTY"
        HighsModelStatusPRESOLVE_ERROR "HighsModelStatus::PRESOLVE_ERROR"
        HighsModelStatusSOLVE_ERROR "HighsModelStatus::SOLVE_ERROR"
        HighsModelStatusPOSTSOLVE_ERROR "HighsModelStatus::POSTSOLVE_ERROR"
        HighsModelStatusPRIMAL_INFEASIBLE "HighsModelStatus::PRIMAL_INFEASIBLE"
        HighsModelStatusPRIMAL_UNBOUNDED "HighsModelStatus::PRIMAL_UNBOUNDED"
        HighsModelStatusOPTIMAL "HighsModelStatus::OPTIMAL"
        HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND "HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND"
        HighsModelStatusREACHED_TIME_LIMIT "HighsModelStatus::REACHED_TIME_LIMIT"
        HighsModelStatusREACHED_ITERATION_LIMIT "HighsModelStatus::REACHED_ITERATION_LIMIT"

cdef extern from "HighsStatus.h" nogil:
    cdef enum HighsStatus:
        HighsStatusOK "HighsStatus::OK"
        HighsStatusWarning "HighsStatus::Warning"
        HighsStatusError "HighsStatus::Error"

    string HighsStatusToString(HighsStatus status)

cdef extern from "LoadProblem.h" nogil:
    HighsStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp)

cdef void reportLpStatsOrError(FILE* output, int message_level, const HighsStatus read_status, const HighsLp& lp):
    if read_status == HighsStatusError:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n")
    else:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "LP       : %s\n", lp.model_name_.c_str())
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Rows     : %d\n", lp.numRow_)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Cols     : %d\n", lp.numCol_)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Nonzeros : %d\n", lp.Avalue_.size())
        if lp.numInt_:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "Integer  : %d\n", lp.numInt_)

cdef extern from "HighsInfo.h" nogil:
    cdef cppclass HighsInfo:
        int simplex_iteration_count
        int ipm_iteration_count
        int crossover_iteration_count
        int primal_status
        int dual_status
        double objective_function_value
        int num_primal_infeasibilities
        double max_primal_infeasibility
        double sum_primal_infeasibilities
        int num_dual_infeasibilities
        double max_dual_infeasibility
        double sum_dual_infeasibilities

cdef extern from "Highs.h":
    cdef cppclass Highs:
        HighsStatus passHighsOptions(const HighsOptions& options)
        HighsStatus passModel(const HighsLp& lp)
        HighsStatus run()
        HighsStatus setHighsLogfile(FILE* logfile)
        HighsStatus setHighsOutput(FILE* output)
        HighsStatus writeHighsOptions(const string filename, const bool report_only_non_default_values = true)

        # split up for cython below
        #const HighsModelStatus& getModelStatus(const bool scaled_model = False) const
        const HighsModelStatus & getModelStatus() const
        const HighsModelStatus & getModelStatus(const bool scaled_model) const

        const HighsInfo& getHighsInfo() const
        string highsModelStatusToString(const HighsModelStatus model_status) const
        #HighsStatus getHighsInfoValue(const string& info, int& value)
        HighsStatus getHighsInfoValue(const string& info, double& value) const
        const HighsOptions& getHighsOptions() const

        # split up for cython
        #HighsStatus writeSolution(const string filename, const bool pretty = false) const
        HighsStatus writeSolution(const string filename) const
        HighsStatus writeSolution(const string filename, const bool pretty) const

cdef void reportSolvedLpStats(FILE* output, int message_level, const HighsStatus run_status, const Highs& highs):
    cdef string statusname
    cdef HighsModelStatus model_status
    cdef HighsModelStatus scaled_model_status
    cdef HighsInfo highs_info
    cdef double objective_function_value = 0 # initialized but written over for cython to stop complaining
    cdef const HighsOptions * options

    if run_status == HighsStatusError:
        statusname = HighsStatusToString(run_status)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n", statusname.c_str())
    else:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "\n")
        model_status = highs.getModelStatus()
        scaled_model_status = highs.getModelStatus(True)
        highs_info = highs.getHighsInfo()
        if model_status != scaled_model_status:
            if scaled_model_status == HighsModelStatusOPTIMAL:
                HighsPrintMessage(output, message_level, ML_ALWAYS,
                                  "Primal infeasibility: %10.3e (%d)\n",
                                  highs_info.max_primal_infeasibility,
                                  highs_info.num_primal_infeasibilities);
                HighsPrintMessage(output, message_level, ML_ALWAYS,
                                  "Dual   infeasibility: %10.3e (%d)\n",
                                  highs_info.max_dual_infeasibility,
                                  highs_info.num_dual_infeasibilities);
                model_status = scaled_model_status;

        HighsPrintMessage(output, message_level, ML_ALWAYS, "Model   status      : %s\n", highs.highsModelStatusToString(model_status).c_str())
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Simplex   iterations: %d\n", highs_info.simplex_iteration_count);
        if highs_info.ipm_iteration_count:
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "IPM       iterations: %d\n",
                              highs_info.ipm_iteration_count)
        if highs_info.crossover_iteration_count:
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "Crossover iterations: %d\n",
                              highs_info.crossover_iteration_count)
        if model_status == HighsModelStatusOPTIMAL:
            highs.getHighsInfoValue("objective_function_value".encode(), objective_function_value)
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "Objective value     : %13.6e\n",
                              objective_function_value)

        # Possibly write the solution to a file
        options = &highs.getHighsOptions()
        if options.write_solution_to_file:
            highs.writeSolution(options.solution_file, options.write_solution_pretty)


cdef HighsStatus callLpSolver(const HighsOptions& options, const HighsLp& lp, FILE* output, int message_level, bool run_quiet):
    # Solve LP case.
    cdef Highs highs
    cdef HighsStatus return_status = highs.passHighsOptions(options)
    if return_status != HighsStatusOK:
        if return_status == HighsStatusWarning:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "HighsStatus::Warning return from passHighsOptions\n")
        else:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "In main: fail return from passHighsOptions\n")
        return return_status

    if run_quiet:
        highs.setHighsLogfile(NULL)
        highs.setHighsOutput(NULL)

    cdef HighsStatus init_status = highs.passModel(lp)
    if init_status != HighsStatusOK:
        if init_status == HighsStatusWarning:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "HighsStatus::Warning return setting HighsLp\n")
        else:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "Error setting HighsLp\n")
        return HighsStatusError

    highs.writeHighsOptions("".encode())

    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Before calling highs.run()\n")

    # Run HiGHS.
    cdef HighsStatus run_status = highs.run()

    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "After calling highs.run()\n")

    reportSolvedLpStats(output, message_level, run_status, highs)
    return run_status

cdef extern from "HighsMipSolver.h" nogil:
    cdef enum HighsMipStatus:
        HighsMipStatuskOptimal "HighsMipStatus::kOptimal"
        HighsMipStatuskTimeout "HighsMipStatus::kTimeout"
        HighsMipStatuskError "HighsMipStatus::kError"
        HighsMipStatuskNodeOptimal "HighsMipStatus::kNodeOptimal"
        HighsMipStatuskNodeInfeasible "HighsMipStatus::kNodeInfeasible"
        HighsMipStatuskNodeUnbounded "HighsMipStatus::kNodeUnbounded"
        HighsMipStatuskNodeNotOptimal "HighsMipStatus::kNodeNotOptimal"
        HighsMipStatuskNodeError "HighsMipStatus::kNodeError"
        HighsMipStatuskRootNodeOptimal "HighsMipStatus::kRootNodeOptimal"
        HighsMipStatuskRootNodeNotOptimal "HighsMipStatus::kRootNodeNotOptimal"
        HighsMipStatuskRootNodeError "HighsMipStatus::kRootNodeError"
        HighsMipStatuskMaxNodeReached "HighsMipStatus::kMaxNodeReached"
        HighsMipStatuskUnderDevelopment "HighsMipStatus::kUnderDevelopment"
        HighsMipStatuskTreeExhausted "HighsMipStatus::kTreeExhausted"

    cdef cppclass HighsMipSolver:
        HighsMipSolver(const HighsOptions& options, const HighsLp& lp) except +
        HighsMipStatus runMipSolver()

cdef HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp, FILE* output, int message_level, bool run_quiet):
    #cdef HighsMipSolver solver(options, lp)
    cdef unique_ptr[HighsMipSolver] solver = make_unique[HighsMipSolver](options, lp)
    cdef HighsMipStatus status = solver.get().runMipSolver()
    if status == HighsMipStatuskOptimal:
        return HighsStatusOK
    return HighsStatusError

def linprog(model_file, solver, bool run_quiet=True):

    # Parse the inputs and put into char**
    args = {
        b'--model_file': model_file.encode(),
        b'--solver': solver.encode(),
    }
    cdef allocator[char *] ptr_al
    cdef unique_ptr[char *] argv
    argv.reset(ptr_al.allocate(len(args)*2+1))
    argv.get()[0] = 'highs' # name of program in argv[0]
    for ii, (k, v) in enumerate(args.items()):
        argv.get()[2*ii+1] = k
        argv.get()[2*ii+2] = v

    # Load user options.
    cdef HighsOptions options
    cdef bool options_ok = loadOptions(len(args)*2+1, argv.get(), options)
    if not options_ok:
        return 0

    # Set message level.
    cdef FILE* output = options.output
    cdef int message_level = options.message_level

    #cdef bool run_quiet = True #False
    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "In main: running highs.run() quietly\n")
    output = options.output
    message_level = options.message_level

    # Load problem.
    cdef HighsLp lp
    cdef HighsStatus read_status = loadLpFromFile(options, lp)
    reportLpStatsOrError(output, message_level, read_status, lp)
    if read_status == HighsStatusError:
        return <int>HighsStatusError

    # Run LP or MIP solver.
    cdef HighsStatus run_status = HighsStatusError
    if not options.mip:
        run_status = callLpSolver(options, lp, output, message_level, run_quiet)
    else:
        run_status = callMipSolver(options, lp, output, message_level, run_quiet)

    return <int>run_status
