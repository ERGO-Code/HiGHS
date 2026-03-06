/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HipoExtrasCApi.cpp
 * @brief C-style API implementation for the HiPO library dependencies
 */

#include "HipoExtrasCApi.h"

extern "C" {

HIPO_EXTRAS_API int hipo_extras_get_abi_version(void) {
  return HIPO_EXTRAS_ABI_VERSION;
}

HIPO_EXTRAS_API const char* hipo_extras_get_version(void) {
  return HIPO_EXTRAS_VERSION;
}

}  // extern "C"

// Note: We use extern "C" here to disable C++ name mangling, allowing
// GetProcAddress/dlsym to find the function by its simple name "hipo_solve_lp".
// While extern "C" is typically for C-compatible interfaces, using C++ types
// (references, classes) works here because:
// 1. Both highspy and highspy-extras are built with the same C++ compiler/ABI
// 2. The ABI version check ensures struct layouts match between builds
// 3. We only need C linkage for symbol lookup, not C type compatibility

extern "C" HIPO_EXTRAS_API int hipo_extras_metis_set_default_options(
    idx_t* options) {
  return Highs_METIS_SetDefaultOptions(options);
}

extern "C" HIPO_EXTRAS_API int hipo_extras_metis_nodend(
    idx_t* nvtxs, const idx_t* xadj, const idx_t* adjncy, idx_t* vwgt,
    idx_t* options, idx_t* perm, idx_t* iperm) {
  return Highs_METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

extern "C" HIPO_EXTRAS_API void hipo_extras_amd_defaults(double Control[]) {
  Highs_amd_defaults(Control);
}

extern "C" HIPO_EXTRAS_API int hipo_extras_amd_order(
    amd_int n, const amd_int Ap[], const amd_int Ai[], amd_int P[],
    double Control[], double Info[]) {
  HighsInt status = Highs_amd_order(n, Ap, Ai, P, Control, Info);
  return status;
}

extern "C" HIPO_EXTRAS_API int hipo_extras_genrcm(HighsInt node_num,
                                                  HighsInt adj_num,
                                                  const HighsInt adj_row[],
                                                  const HighsInt adj[],
                                                  HighsInt perm[]) {
  HighsInt status = Highs_genrcm(node_num, adj_num, adj_row, adj, perm);
  return status;
}

// extern "C" HIPO_EXTRAS_API HighsStatus hipo_solve_lp(
//     const HighsOptions& options,
//     HighsTimer& timer,
//     const HighsLp& lp,
//     HighsBasis& highs_basis,
//     HighsSolution& highs_solution,
//     HighsModelStatus& model_status,
//     HighsInfo& highs_info,
//     HighsCallback& callback) {
//   // Call the actual HiPO solver implementation
//   // Note: solveLpHipo is defined in IpxWrapper.cpp when HIPO is defined
//   return solveLpHipo(options, timer, lp, highs_basis, highs_solution,
//                      model_status, highs_info, callback);
// }

// HighsStatus solveLpHipo(HighsLpSolverObject& solver_object) {
//   return solveLpHipo(solver_object.options_, solver_object.timer_,
//                      solver_object.lp_, solver_object.basis_,
//                      solver_object.solution_, solver_object.model_status_,
//                      solver_object.highs_info_, solver_object.callback_);
// }

#ifdef HIPO_USES_OPENBLAS
// function to set number of threads of openblas
extern "C" {
void openblas_set_num_threads(int num_threads);
}
#endif

// HighsStatus solveLpHipo(const HighsOptions& options, HighsTimer& timer,
//                         const HighsLp& lp, HighsBasis& highs_basis,
//                         HighsSolution& highs_solution,
//                         HighsModelStatus& model_status, HighsInfo&
//                         highs_info, HighsCallback& callback) {
//   // Use HiPO
//   //
//   // Can return HighsModelStatus (HighsStatus) values:
//   //
//   // 1. kSolveError (kError) if various unlikely solution errors occur
//   //
//   // 2. kTimeLimit (kWarning) if time limit is reached
//   //
//   // 3. kIterationLimit (kWarning) if iteration limit is reached
//   //
//   // 4. kUnknown (kWarning) if HiPO makes no progress or if
//   // IPM/crossover are imprecise
//   //
//   // 5. kInfeasible (kOk) if HiPO identifies primal infeasibility
//   //
//   // 6. kUnboundedOrInfeasible (kOk) if HiPO identifies dual
//   // infeasibility
//   //
//   // kOptimal (kOk) if HiPO/crossover identify optimality
//   //
//   // With a non-error return, if just HiPO has been run then a
//   // non-vertex primal solution is obtained; if crossover has been run
//   // then a basis and primal+dual solution are obtained.
//   //
//   //
//   // Indicate that there is no valid primal solution, dual solution or basis
//   highs_basis.valid = false;
//   highs_solution.value_valid = false;
//   highs_solution.dual_valid = false;
//   // Indicate that no imprecise solution has (yet) been found
//   resetModelStatusAndHighsInfo(model_status, highs_info);

// #ifdef HIPO_USES_OPENBLAS
//   // force openblas to run in serial, for determinism and better performance
//   openblas_set_num_threads(1);
// #endif

//   // Create solver instance
//   hipo::Solver hipo{};
//   // This creates ipx::LpSolver ipx_lps_, in case HiPO has to switch
//   // to IPX, so use the current HiGHS time as an offset for the
//   // ipx_lps.control_ elapsed time
//   hipo.setIpxTimerOffset(timer.read());

//   hipo::Options hipo_options{};

//   hipo_options.display = true;
//   if (!options.output_flag | !options.log_to_console)
//     hipo_options.display = false;

//   hipo_options.log_options = &options.log_options;

//   // Debug option is already considered through log_options.log_dev_level in
//   // hipo::LogHighs::debug

//   hipo_options.timeless_log = options.timeless_log;
//   hipo_options.feasibility_tol =
//   std::min(options.primal_feasibility_tolerance,
//                                           options.dual_feasibility_tolerance);
//   hipo_options.optimality_tol = options.ipm_optimality_tolerance;
//   hipo_options.crossover_tol = options.start_crossover_tolerance;

//   if (options.kkt_tolerance != kDefaultKktTolerance) {
//     hipo_options.feasibility_tol = options.kkt_tolerance;
//     hipo_options.optimality_tol = 1e-1 * options.kkt_tolerance;
//     hipo_options.crossover_tol = 1e-1 * options.kkt_tolerance;
//     highsLogUser(options.log_options, HighsLogType::kInfo,
//                  "IpxWrapper: feasibility_tol = %g; optimality_tol = %g; "
//                  "crossover_tol = %g\n",
//                  hipo_options.feasibility_tol, hipo_options.optimality_tol,
//                  hipo_options.crossover_tol);
//   }

//   // hipo uses same timer as highs, so it is fine to pass the same time limit
//   hipo_options.time_limit = options.time_limit;

//   hipo_options.max_iter =
//       options.ipm_iteration_limit - highs_info.ipm_iteration_count;

//   if (options.run_crossover == kHighsOnString)
//     hipo_options.crossover = hipo::kOptionCrossoverOn;
//   else if (options.run_crossover == kHighsOffString)
//     hipo_options.crossover = hipo::kOptionCrossoverOff;
//   else {
//     assert(options.run_crossover == kHighsChooseString);
//     hipo_options.crossover = hipo::kOptionCrossoverChoose;
//   }

//   // Potentially control if ipx is used for refinement and if it is displayed
//   // hipo_options.refine_with_ipx = true;
//   hipo_options.display_ipx = true;

//   // if option parallel is on, it can be refined by option hipo_parallel_type
//   if (options.parallel == kHighsOnString) {
//     if (options.hipo_parallel_type == kHipoTreeString)
//       hipo_options.parallel = hipo::kOptionParallelTreeOnly;
//     else if (options.hipo_parallel_type == kHipoNodeString)
//       hipo_options.parallel = hipo::kOptionParallelNodeOnly;
//     else if (options.hipo_parallel_type == kHipoBothString)
//       hipo_options.parallel = hipo::kOptionParallelOn;
//     else {
//       highsLogUser(options.log_options, HighsLogType::kError,
//                    "Unknown value of option %s\n",
//                    kHipoParallelString.c_str());
//       model_status = HighsModelStatus::kSolveError;
//       return HighsStatus::kError;
//     }
//   }
//   // otherwise, option hipo_parallel_type is ignored
//   else if (options.parallel == kHighsOffString)
//     hipo_options.parallel = hipo::kOptionParallelOff;
//   else {
//     assert(options.parallel == kHighsChooseString);
//     hipo_options.parallel = hipo::kOptionParallelChoose;
//   }

//   // Parse hipo_system option
//   if (options.hipo_system == kHipoAugmentedString) {
//     hipo_options.nla = hipo::kOptionNlaAugmented;
//   } else if (options.hipo_system == kHipoNormalEqString) {
//     hipo_options.nla = hipo::kOptionNlaNormEq;
//   } else if (options.hipo_system == kHighsChooseString) {
//     hipo_options.nla = hipo::kOptionNlaChoose;
//   } else {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Unknown value of option %s\n", kHipoSystemString.c_str());
//     model_status = HighsModelStatus::kSolveError;
//     return HighsStatus::kError;
//   }

//   // Reordering heuristic
//   if (options.hipo_ordering != kHipoMetisString &&
//       options.hipo_ordering != kHipoAmdString &&
//       options.hipo_ordering != kHipoRcmString &&
//       options.hipo_ordering != kHighsChooseString) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Unknown value of option %s\n",
//                  kHipoOrderingString.c_str());
//     model_status = HighsModelStatus::kSolveError;
//     return HighsStatus::kError;
//   }
//   hipo_options.ordering = options.hipo_ordering;

//   // block size option
//   hipo_options.block_size = options.hipo_block_size;

//   hipo.setOptions(hipo_options);
//   hipo.setTimer(timer);
//   hipo.setCallback(callback);

//   // Transform problem to correct formulation
//   hipo::Int num_col, num_row;
//   std::vector<double> obj, rhs, lower, upper, Aval;
//   std::vector<hipo::Int> Aptr, Aind;
//   std::vector<char> constraints;
//   double offset;
//   fillInIpxData(lp, num_col, num_row, offset, obj, lower, upper, Aptr, Aind,
//                 Aval, rhs, constraints);
//   highsLogUser(options.log_options, HighsLogType::kInfo,
//                "HiPO model has %" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
//                " columns and %" HIGHSINT_FORMAT " nonzeros\n",
//                num_row, num_col, Aptr[num_col]);

//   // Load the problem
//   hipo::Int load_status = hipo.load(
//       num_col, num_row, obj.data(), rhs.data(), lower.data(), upper.data(),
//       Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset);

//   if (load_status) {
//     model_status = HighsModelStatus::kSolveError;
//     return HighsStatus::kError;
//   }

//   hipo.solve();

//   // const bool report_solve_data =
//   //    kHighsAnalysisLevelSolverSummaryData & options.highs_analysis_level;

//   // Differently from IPX, HiPO returns a single status. So, dealing with
//   // statuses is a bit different.
//   // hipo.solved(), hipo.stopped(), hipo.failed() can be used to query if the
//   // status belongs to the solved, stopped or failed group.
//   // If primal-dual feasible solution is found (non-vertex solution), then
//   the
//   // status is kStatusPDfeas.
//   // If crossover is successful, then the status is kStatusBasic.
//   // Otherwise, the specific crossover status can be accessed through the
//   // ipx_info stored in hipo_info.

//   // Get solver and solution information.
//   const hipo::Info hipo_info = hipo.getInfo();
//   hipo::Status solve_status = hipo_info.status;
//   highs_info.ipm_iteration_count +=
//       hipo_info.ipm_iter + hipo_info.ipx_info.iter;
//   highs_info.crossover_iteration_count +=
//   hipo_info.ipx_info.updates_crossover;

//   // Report hipo status
//   const HighsStatus solve_return_status =
//       reportHipoStatus(options, solve_status, hipo);
//   if (solve_return_status == HighsStatus::kError) {
//     model_status = HighsModelStatus::kSolveError;
//     return HighsStatus::kError;
//   }

//   // Report crossover status
//   const HighsStatus crossover_return_status =
//       reportHipoCrossoverStatus(options,
//       hipo_info.ipx_info.status_crossover);
//   if (crossover_return_status == HighsStatus::kError) {
//     model_status = HighsModelStatus::kSolveError;
//     return HighsStatus::kError;
//   }

//   // Failures should have been handled. Status should be stopper or solved.
//   if (ipxStatusError(!hipo.solved() && !hipo.stopped(), options, "Hipo",
//                      "status should be solved or stopped but value is",
//                      solve_status))
//     return HighsStatus::kError;

//   if (hipo.stopped()) {
//     const HighsModelStatus local_model_status = HighsModelStatus::kUnknown;

//     getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
//                              hipo, local_model_status, highs_solution);

//     // For crossover
//     if (illegalIpxStoppedCrossoverStatus(hipo_info.ipx_info, options))
//       return HighsStatus::kError;
//     // Can stop and reach time limit
//     if (hipo_info.ipx_info.status_crossover == IPX_STATUS_time_limit) {
//       model_status = HighsModelStatus::kTimeLimit;
//       return HighsStatus::kWarning;
//     }

//     // if crossover didn't time out, then solver can only stop as follows
//     if (solve_status == hipo::kStatusUserInterrupt) {
//       model_status = HighsModelStatus::kInterrupt;
//       return HighsStatus::kWarning;
//     } else if (solve_status == hipo::kStatusTimeLimit) {
//       model_status = HighsModelStatus::kTimeLimit;
//       return HighsStatus::kWarning;
//     } else if (solve_status == hipo::kStatusMaxIter) {
//       model_status = HighsModelStatus::kIterationLimit;
//       return HighsStatus::kWarning;
//     } else if (solve_status == hipo::kStatusNoProgress) {
//       reportHipoNoProgress(options, hipo_info);
//       model_status = HighsModelStatus::kUnknown;
//       return HighsStatus::kWarning;
//     } else {
//       assert(1 == 0);
//     }
//   }

//   // Stopped status should have been handled. Status should be solved.
//   if (ipxStatusError(!hipo.solved(), options, "Hipo",
//                      "status should be solved but value is", solve_status))
//     return HighsStatus::kError;

//   // primal/dual infeasible
//   if (solve_status == hipo::kStatusPrimalInfeasible ||
//       solve_status == hipo::kStatusDualInfeasible) {
//     if (solve_status == hipo::kStatusPrimalInfeasible)
//       model_status = HighsModelStatus::kInfeasible;
//     else
//       model_status = HighsModelStatus::kUnboundedOrInfeasible;

//     getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
//                              hipo, model_status, highs_solution);

//     return HighsStatus::kOk;
//   }

//   // Status should be optimal or imprecise
//   if (ipxStatusError(solve_status != hipo::kStatusPDFeas &&
//                          solve_status != hipo::kStatusBasic &&
//                          solve_status != hipo::kStatusImprecise,
//                      options, "Hipo",
//                      "status should be optimal or imprecise but value is",
//                      solve_status))
//     return HighsStatus::kError;

//   const bool have_basic_solution =
//       hipo_info.ipx_used &&
//       hipo_info.ipx_info.status_crossover != IPX_STATUS_not_run;

//   const bool imprecise_solution =
//       hipo_info.status == hipo::kStatusImprecise ||
//       hipo_info.ipx_info.status_crossover == IPX_STATUS_imprecise;

//   if (have_basic_solution) {
//     IpxSolution hipo_solution;
//     hipo_solution.num_col = num_col;
//     hipo_solution.num_row = num_row;
//     hipo_solution.ipx_col_value.resize(num_col);
//     hipo_solution.ipx_row_value.resize(num_row);
//     hipo_solution.ipx_col_dual.resize(num_col);
//     hipo_solution.ipx_row_dual.resize(num_row);
//     hipo_solution.ipx_row_status.resize(num_row);
//     hipo_solution.ipx_col_status.resize(num_col);

//     hipo::Int errflag = hipo.getBasicSolution(
//         hipo_solution.ipx_col_value, hipo_solution.ipx_row_value,
//         hipo_solution.ipx_row_dual, hipo_solution.ipx_col_dual,
//         hipo_solution.ipx_row_status.data(),
//         hipo_solution.ipx_col_status.data());

//     if (errflag) {
//       highsLogUser(options.log_options, HighsLogType::kError,
//                    "IPX crossover getting basic solution: flag = %d\n",
//                    (int)errflag);
//       return HighsStatus::kError;
//     }
//     // Convert the IPX basic solution to a HiGHS basic solution
//     HighsStatus status = ipxBasicSolutionToHighsBasicSolution(
//         options.log_options, lp, rhs, constraints, hipo_solution,
//         highs_basis, highs_solution);
//     if (status != HighsStatus::kOk) {
//       highsLogUser(
//           options.log_options, HighsLogType::kError,
//           "Failed to convert IPX basic solution to Highs basic solution\n");
//       return HighsStatus::kError;
//     }
//   } else {
//     const HighsModelStatus local_model_status =
//         imprecise_solution ? HighsModelStatus::kUnknown
//                            : HighsModelStatus::kOptimal;
//     getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
//                              hipo, local_model_status, highs_solution);
//     assert(!highs_basis.valid);
//   }

//   highs_info.basis_validity =
//       highs_basis.valid ? kBasisValidityValid : kBasisValidityInvalid;

//   HighsStatus return_status;
//   if (imprecise_solution) {
//     model_status = HighsModelStatus::kUnknown;
//     return_status = HighsStatus::kWarning;
//   } else {
//     model_status = HighsModelStatus::kOptimal;
//     return_status = HighsStatus::kOk;
//   }
//   return return_status;
// }

// void reportHipoNoProgress(const HighsOptions& options,
//                           const hipo::Info& hipo_info) {
//   highsLogUser(options.log_options, HighsLogType::kWarning,
//                "No progress: primal objective value       = %11.4g\n",
//                hipo_info.p_obj);
//   highsLogUser(options.log_options, HighsLogType::kWarning,
//                "No progress: max absolute primal residual = %11.4g\n",
//                hipo_info.p_res_abs);
//   highsLogUser(options.log_options, HighsLogType::kWarning,
//                "No progress: max absolute   dual residual = %11.4g\n",
//                hipo_info.d_res_abs);
// }

// void getHipoNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
//                               const hipo::Int num_col, const hipo::Int
//                               num_row, const std::vector<double>& rhs, const
//                               std::vector<char>& constraint_type, const
//                               hipo::Solver& hipo, const HighsModelStatus
//                               model_status, HighsSolution& highs_solution) {
//   std::vector<double> x(num_col);
//   std::vector<double> xl(num_col);
//   std::vector<double> xu(num_col);
//   std::vector<double> zl(num_col);
//   std::vector<double> zu(num_col);
//   std::vector<double> slack(num_row);
//   std::vector<double> y(num_row);

//   hipo.getInteriorSolution(x, xl, xu, slack, y, zl, zu);
//   ipxSolutionToHighsSolution(options, lp, rhs, constraint_type, num_col,
//                              num_row, x, slack, y, zl, zu, highs_solution);
// }

// HighsStatus reportHipoStatus(const HighsOptions& options,
//                              const hipo::Int status, const hipo::Solver&
//                              hipo) {
//   if (hipo.solved()) {
//     highsLogUser(options.log_options, HighsLogType::kInfo, "Hipo: Solved\n");
//     return HighsStatus::kOk;
//   }

//   // these are warnings
//   else if (status == hipo::kStatusTimeLimit) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Time limit\n");
//     return HighsStatus::kWarning;
//   } else if (status == hipo::kStatusUserInterrupt) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: User interrupt\n");
//     return HighsStatus::kWarning;
//   } else if (status == hipo::kStatusMaxIter) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Reached maximum iterations\n");
//     return HighsStatus::kWarning;
//   } else if (status == hipo::kStatusNoProgress) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: No progress\n");
//     return HighsStatus::kWarning;
//   } else if (status == hipo::kStatusImprecise) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Imprecise solution\n");
//     return HighsStatus::kWarning;
//   }

//   // these are errors
//   else if (status == hipo::kStatusError) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Internal error\n");
//   } else if (status == hipo::kStatusOverflow) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Integer overflow\n");
//   } else if (status == hipo::kStatusErrorAnalyse) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Error in analyse phase\n");
//   } else if (status == hipo::kStatusErrorFactorise) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Error in factorise phase\n");
//   } else if (status == hipo::kStatusErrorSolve) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Error in solve phase\n");
//   } else if (status == hipo::kStatusBadModel) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Invalid model\n");
//   } else {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Unrecognized status\n");
//   }
//   return HighsStatus::kError;
// }

// HighsStatus reportHipoCrossoverStatus(const HighsOptions& options,
//                                       const ipx::Int status) {
//   if (status == IPX_STATUS_not_run) {
//     if (options.run_crossover == kHighsOnString) {
//       // Warn if crossover not run and run_crossover option is "on"
//       highsLogUser(options.log_options, HighsLogType::kWarning,
//                    "Hipo: Crossover not run\n");
//       return HighsStatus::kWarning;
//     }
//     return HighsStatus::kOk;
//   } else if (status == IPX_STATUS_optimal) {
//     highsLogUser(options.log_options, HighsLogType::kInfo,
//                  "Hipo: Crossover optimal\n");
//     return HighsStatus::kOk;
//   } else if (status == IPX_STATUS_imprecise) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover imprecise\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_primal_infeas) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover primal infeasible\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_dual_infeas) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover dual infeasible\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_user_interrupt) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover user interrupt\n");
//     return HighsStatus::kOk;
//   } else if (status == IPX_STATUS_time_limit) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover reached time limit\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_iter_limit) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover reached iteration limit\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_no_progress) {
//     highsLogUser(options.log_options, HighsLogType::kWarning,
//                  "Hipo: Crossover no progress\n");
//     return HighsStatus::kWarning;
//   } else if (status == IPX_STATUS_failed) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Crossover failed\n");
//     return HighsStatus::kError;
//   } else if (status == IPX_STATUS_debug) {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Crossover debug\n");
//     return HighsStatus::kError;
//   } else {
//     highsLogUser(options.log_options, HighsLogType::kError,
//                  "Hipo: Crossover unrecognised status\n");
//     return HighsStatus::kError;
//   }
//   return HighsStatus::kError;
// }