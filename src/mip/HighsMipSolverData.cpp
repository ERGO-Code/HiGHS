/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolverData.h"

#include <random>

#include "lp_data/HighsLpUtils.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsRedcostFixing.h"
#include "pdqsort/pdqsort.h"
#include "presolve/HAggregator.h"
#include "presolve/HPresolve.h"
#include "util/HighsIntegers.h"

bool HighsMipSolverData::checkSolution(const std::vector<double>& solution) {
  for (HighsInt i = 0; i != mipsolver.model_->num_col_; ++i) {
    if (solution[i] < mipsolver.model_->col_lower_[i] - feastol) return false;
    if (solution[i] > mipsolver.model_->col_upper_[i] + feastol) return false;
    if (mipsolver.variableType(i) == HighsVarType::kInteger &&
        std::abs(solution[i] - std::floor(solution[i] + 0.5)) > feastol)
      return false;
  }

  for (HighsInt i = 0; i != mipsolver.model_->num_row_; ++i) {
    double rowactivity = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      rowactivity += solution[ARindex_[j]] * ARvalue_[j];

    if (rowactivity > mipsolver.rowUpper(i) + feastol) return false;
    if (rowactivity < mipsolver.rowLower(i) - feastol) return false;
  }

  return true;
}

bool HighsMipSolverData::trySolution(const std::vector<double>& solution,
                                     char source) {
  if (int(solution.size()) != mipsolver.model_->num_col_) return false;

  HighsCDouble obj = 0;

  for (HighsInt i = 0; i != mipsolver.model_->num_col_; ++i) {
    if (solution[i] < mipsolver.model_->col_lower_[i] - feastol) return false;
    if (solution[i] > mipsolver.model_->col_upper_[i] + feastol) return false;
    if (mipsolver.variableType(i) == HighsVarType::kInteger &&
        std::abs(solution[i] - std::floor(solution[i] + 0.5)) > feastol)
      return false;

    obj += mipsolver.colCost(i) * solution[i];
  }

  for (HighsInt i = 0; i != mipsolver.model_->num_row_; ++i) {
    double rowactivity = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      rowactivity += solution[ARindex_[j]] * ARvalue_[j];

    if (rowactivity > mipsolver.rowUpper(i) + feastol) return false;
    if (rowactivity < mipsolver.rowLower(i) - feastol) return false;
  }

  return addIncumbent(solution, double(obj), source);
}

bool HighsMipSolverData::moreHeuristicsAllowed() {
  // in the beginning of the search and in sub-MIP heuristics we only allow
  // what is proportionally for the currently spent effort plus an initial
  // offset. This is because in a sub-MIP we usually do a truncated search and
  // therefore should not extrapolate the time we spent for heuristics as in
  // the other case. Moreover, since we estimate the total effort for
  // exploring the tree based on the weight of the already pruned nodes, the
  // estimated effort the is not expected to be a good prediction in the
  // beginning.
  if (mipsolver.submip) {
    return heuristic_lp_iterations < total_lp_iterations * heuristic_effort;
  } else if (pruned_treeweight < 1e-3 &&
             num_leaves - num_leaves_before_run < 10) {
    // in the main MIP solver allow an initial offset of 10000 heuristic LP
    // iterations
    if (heuristic_lp_iterations <
        total_lp_iterations * heuristic_effort + 10000)
      return true;
  } else if (heuristic_lp_iterations <
             100000 + ((total_lp_iterations - heuristic_lp_iterations -
                        sb_lp_iterations) >>
                       1)) {
    // compute the node LP iterations in the current run as only those should be
    // used when estimating the total required LP iterations to complete the
    // search
    int64_t heur_iters_curr_run =
        heuristic_lp_iterations - heuristic_lp_iterations_before_run;
    int64_t sb_iters_curr_run = sb_lp_iterations - sb_lp_iterations_before_run;
    int64_t node_iters_curr_run = total_lp_iterations -
                                  total_lp_iterations_before_run -
                                  heur_iters_curr_run - sb_iters_curr_run;
    // now estimate the total fraction of LP iterations that we have spent on
    // heuristics by assuming the node iterations of the current run will
    // grow proportional to the pruned weight of the current tree and the
    // iterations spent for anything else are just added as an offset
    double total_heuristic_effort_estim =
        heuristic_lp_iterations /
        ((total_lp_iterations - node_iters_curr_run) +
         node_iters_curr_run / std::max(1e-3, double(pruned_treeweight)));
    // since heuristics help most in the beginning of the search, we want to
    // spent the time we have for heuristics in the first 80% of the tree
    // exploration. Additionally we want to spent the proportional effort
    // of heuristics that is allowed in the the first 30% of tree exploration as
    // fast as possible, which is why we have the max(0.3/0.8,...).
    // Hence, in the first 30% of the tree exploration we allow to spent all
    // effort available for heuristics in that part of the search as early as
    // possible, whereas after that we allow the part that is proportionally
    // adequate when we want to spent all available time in the first 80%.
    if (total_heuristic_effort_estim <
        std::max(0.3 / 0.8, std::min(double(pruned_treeweight), 0.8) / 0.8) *
            heuristic_effort) {
      // printf(
      //     "heuristic lp iterations: %ld, total_lp_iterations: %ld, "
      //     "total_heur_effort_estim = %.3f%%\n",
      //     heuristic_lp_iterations, total_lp_iterations,
      //     total_heuristic_effort_estim);
      return true;
    }
  }

  return false;
}

void HighsMipSolverData::removeFixedIndices() {
  integral_cols.erase(
      std::remove_if(integral_cols.begin(), integral_cols.end(),
                     [&](HighsInt col) { return domain.isFixed(col); }),
      integral_cols.end());
  integer_cols.erase(
      std::remove_if(integer_cols.begin(), integer_cols.end(),
                     [&](HighsInt col) { return domain.isFixed(col); }),
      integer_cols.end());
  implint_cols.erase(
      std::remove_if(implint_cols.begin(), implint_cols.end(),
                     [&](HighsInt col) { return domain.isFixed(col); }),
      implint_cols.end());
  continuous_cols.erase(
      std::remove_if(continuous_cols.begin(), continuous_cols.end(),
                     [&](HighsInt col) { return domain.isFixed(col); }),
      continuous_cols.end());
}

void HighsMipSolverData::init() {
  postSolveStack.initializeIndexMaps(mipsolver.model_->num_row_,
                                     mipsolver.model_->num_col_);
  mipsolver.orig_model_ = mipsolver.model_;
  if (mipsolver.clqtableinit) cliquetable.buildFrom(*mipsolver.clqtableinit);
  if (mipsolver.implicinit) implications.buildFrom(*mipsolver.implicinit);
  feastol = mipsolver.options_mip_->mip_feasibility_tolerance;
  epsilon = mipsolver.options_mip_->small_matrix_value;
  heuristic_effort = mipsolver.options_mip_->mip_heuristic_effort;
  detectSymmetries = mipsolver.options_mip_->mip_detect_symmetry;

  firstlpsolobj = -kHighsInf;
  rootlpsolobj = -kHighsInf;
  analyticCenterComputed = false;
  numRestarts = 0;
  numRestartsRoot = 0;
  numImprovingSols = 0;
  pruned_treeweight = 0;
  avgrootlpiters = 0;
  num_nodes = 0;
  num_nodes_before_run = 0;
  num_leaves = 0;
  num_leaves_before_run = 0;
  total_lp_iterations = 0;
  heuristic_lp_iterations = 0;
  sepa_lp_iterations = 0;
  sb_lp_iterations = 0;
  total_lp_iterations_before_run = 0;
  heuristic_lp_iterations_before_run = 0;
  sepa_lp_iterations_before_run = 0;
  sb_lp_iterations_before_run = 0;
  num_disp_lines = 0;
  cliquesExtracted = false;
  rowMatrixSet = false;
  lower_bound = -kHighsInf;
  upper_bound = kHighsInf;
  upper_limit = mipsolver.options_mip_->objective_bound;

  if (mipsolver.options_mip_->mip_report_level == 0)
    dispfreq = 0;
  else if (mipsolver.options_mip_->mip_report_level == 1)
    dispfreq = 2000;
  else
    dispfreq = 100;
}

void HighsMipSolverData::runPresolve() {
#ifdef HIGHS_DEBUGSOL
  bool debugSolActive = false;
  std::swap(debugSolution.debugSolActive, debugSolActive);
#endif

  mipsolver.timer_.start(mipsolver.timer_.presolve_clock);
  presolve::HPresolve presolve;
  presolve.setInput(mipsolver);
  mipsolver.modelstatus_ = presolve.run(postSolveStack);
  mipsolver.timer_.stop(mipsolver.timer_.presolve_clock);

#ifdef HIGHS_DEBUGSOL
  debugSolution.debugSolActive = debugSolActive;
  if (debugSolution.debugSolActive) debugSolution.registerDomain(domain);
  assert(!debugSolution.debugSolActive ||
         checkSolution(debugSolution.debugSolution));
#endif
}

void HighsMipSolverData::runSetup() {
  const HighsLp& model = *mipsolver.model_;

  last_disptime = -kHighsInf;

  // transform the objective limit to the current model
  upper_limit -= mipsolver.model_->offset_;
  lower_bound -= mipsolver.model_->offset_;
  upper_bound -= mipsolver.model_->offset_;

  if (mipsolver.numCol() == 0) addIncumbent(std::vector<double>(), 0, 'P');

  redcostfixing = HighsRedcostFixing();
  pseudocost = HighsPseudocost(mipsolver);
  nodequeue.setNumCol(mipsolver.numCol());

  continuous_cols.clear();
  integer_cols.clear();
  implint_cols.clear();
  integral_cols.clear();

  rowMatrixSet = false;
  if (!rowMatrixSet) {
    rowMatrixSet = true;
    highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                         model.a_matrix_.index_, model.a_matrix_.value_,
                         ARstart_, ARindex_, ARvalue_);
    uplocks.resize(model.num_col_);
    downlocks.resize(model.num_col_);
    for (HighsInt i = 0; i != model.num_col_; ++i) {
      HighsInt start = model.a_matrix_.start_[i];
      HighsInt end = model.a_matrix_.start_[i + 1];
      for (HighsInt j = start; j != end; ++j) {
        HighsInt row = model.a_matrix_.index_[j];

        if (model.row_lower_[row] != -kHighsInf) {
          if (model.a_matrix_.value_[j] < 0)
            ++uplocks[i];
          else
            ++downlocks[i];
        }
        if (model.row_upper_[row] != kHighsInf) {
          if (model.a_matrix_.value_[j] < 0)
            ++downlocks[i];
          else
            ++uplocks[i];
        }
      }
    }
  }

  rowintegral.resize(mipsolver.model_->num_row_);

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.model_->num_row_);
  for (HighsInt i = 0; i != mipsolver.model_->num_row_; ++i) {
    double maxabsval = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];
    bool integral = true;
    for (HighsInt j = start; j != end; ++j) {
      if (integral) {
        if (mipsolver.variableType(ARindex_[j]) == HighsVarType::kContinuous)
          integral = false;
        else {
          double intval = std::floor(ARvalue_[j] + 0.5);
          if (std::abs(ARvalue_[j] - intval) > epsilon) integral = false;
        }
      }

      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));
    }

    if (integral) {
      if (presolvedModel.row_lower_[i] != -kHighsInf)
        presolvedModel.row_lower_[i] =
            std::ceil(presolvedModel.row_lower_[i] - feastol);

      if (presolvedModel.row_upper_[i] != kHighsInf)
        presolvedModel.row_upper_[i] =
            std::floor(presolvedModel.row_upper_[i] + feastol);
    }

    rowintegral[i] = integral;
    maxAbsRowCoef[i] = maxabsval;
  }

  // compute row activities and propagate all rows once
  domain.computeRowActivities();
  domain.propagate();
  if (domain.infeasible()) {
    mipsolver.modelstatus_ = HighsModelStatus::kInfeasible;
    lower_bound = kHighsInf;
    pruned_treeweight = 1.0;
    return;
  }

  if (model.num_col_ == 0) {
    mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    return;
  }

  if (checkLimits()) return;
  // extract cliques if they have not been extracted before

  for (HighsInt col : domain.getChangedCols())
    implications.cleanupVarbounds(col);
  domain.clearChangedCols();

  lp.getLpSolver().setOptionValue("presolve", "off");
  // lp.getLpSolver().setOptionValue("dual_simplex_cost_perturbation_multiplier",
  // 0.0); lp.getLpSolver().setOptionValue("parallel", "on");
  lp.getLpSolver().setOptionValue("simplex_initial_condition_check", false);

  checkObjIntegrality();
  basisTransfer();
  rootlpsol.clear();
  firstlpsol.clear();
  HighsInt numBin = 0;

  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    switch (mipsolver.variableType(i)) {
      case HighsVarType::kContinuous:
        continuous_cols.push_back(i);
        break;
      case HighsVarType::kImplicitInteger:
        implint_cols.push_back(i);
        integral_cols.push_back(i);
        break;
      case HighsVarType::kInteger:
        integer_cols.push_back(i);
        integral_cols.push_back(i);
        numBin += ((mipsolver.model_->col_lower_[i] == 0.0) &
                   (mipsolver.model_->col_upper_[i] == 1.0));
    }
  }
  numintegercols = integer_cols.size();
  detectSymmetries = detectSymmetries && numBin > 0;

  if (numRestarts == 0) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 // clang-format off
               "\nSolving MIP model with:\n"
               "   %" HIGHSINT_FORMAT " rows\n"
               "   %" HIGHSINT_FORMAT " cols (%" HIGHSINT_FORMAT" binary, %" HIGHSINT_FORMAT " integer, %" HIGHSINT_FORMAT" implied int., %" HIGHSINT_FORMAT " continuous)\n"
               "   %" HIGHSINT_FORMAT " nonzeros\n",
                 // clang-format on
                 mipsolver.numRow(), mipsolver.numCol(), numBin,
                 numintegercols - numBin, (HighsInt)implint_cols.size(),
                 (HighsInt)continuous_cols.size(), mipsolver.numNonzero());
  } else {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "Model after restart has %" HIGHSINT_FORMAT
                 " rows, %" HIGHSINT_FORMAT " cols (%" HIGHSINT_FORMAT
                 " bin., %" HIGHSINT_FORMAT " int., %" HIGHSINT_FORMAT
                 " impl., %" HIGHSINT_FORMAT " cont.), and %" HIGHSINT_FORMAT
                 " nonzeros\n",
                 mipsolver.numRow(), mipsolver.numCol(), numBin,
                 numintegercols - numBin, (HighsInt)implint_cols.size(),
                 (HighsInt)continuous_cols.size(), mipsolver.numNonzero());
  }

  heuristics.setupIntCols();

#ifdef HIGHS_DEBUGSOL
  if (numRestarts == 0) {
    debugSolution.activate();
    assert(!debugSolution.debugSolActive ||
           checkSolution(debugSolution.debugSolution));
  }
#endif

  symmetries.clear();

  if (detectSymmetries) {
    if (numRestarts == 0)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "\n");

    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "(%4.1fs) Starting symmetry detection\n",
                 mipsolver.timer_.read(mipsolver.timer_.solve_clock));
    HighsSymmetryDetection symDetection;
    symDetection.loadModelAsGraph(mipsolver.mipdata_->presolvedModel,
                                  mipsolver.options_mip_->small_matrix_value);
    symDetection.run(symmetries);
    if (symmetries.numGenerators == 0) {
      detectSymmetries = false;
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "(%4.1fs) No symmetry present\n",
                   mipsolver.timer_.read(mipsolver.timer_.solve_clock));
    } else if (symmetries.orbitopes.size() == 0) {
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "(%4.1fs) Found %" HIGHSINT_FORMAT " generators\n",
                   mipsolver.timer_.read(mipsolver.timer_.solve_clock),
                   symmetries.numGenerators);

    } else {
      if (symmetries.numPerms != 0) {
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                     "(%4.1fs) Found %" HIGHSINT_FORMAT
                     " generators and %" HIGHSINT_FORMAT
                     " full orbitope(s) acting on %" HIGHSINT_FORMAT
                     " columns\n",
                     mipsolver.timer_.read(mipsolver.timer_.solve_clock),
                     symmetries.numPerms, (HighsInt)symmetries.orbitopes.size(),
                     (HighsInt)symmetries.columnToOrbitope.size());
      } else {
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                     "(%4.1fs) Found %" HIGHSINT_FORMAT
                     " full orbitope(s) acting on %" HIGHSINT_FORMAT
                     " columns\n",
                     mipsolver.timer_.read(mipsolver.timer_.solve_clock),
                     (HighsInt)symmetries.orbitopes.size(),
                     (HighsInt)symmetries.columnToOrbitope.size());
      }
      for (HighsOrbitopeMatrix& orbitope : symmetries.orbitopes)
        orbitope.determineOrbitopeType(cliquetable, domain);

      if (!domain.getChangedCols().empty()) {
        domain.propagate();
        if (domain.infeasible()) {
          mipsolver.modelstatus_ = HighsModelStatus::kInfeasible;
          lower_bound = kHighsInf;
          pruned_treeweight = 1.0;
          return;
        }
        domain.clearChangedCols();
      }
    }
  }

  if (numRestarts != 0)
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "\n");
}

double HighsMipSolverData::transformNewIncumbent(
    const std::vector<double>& sol) {
  HighsSolution solution;
  solution.col_value = sol;
  calculateRowValues(*mipsolver.model_, solution);

  postSolveStack.undoPrimal(*mipsolver.options_mip_, solution);
  calculateRowValues(*mipsolver.orig_model_, solution);
  bool allow_try_again = true;
try_again:

  // compute the objective value in the original space
  double bound_violation_ = 0;
  double row_violation_ = 0;
  double integrality_violation_ = 0;

  HighsCDouble obj = mipsolver.orig_model_->offset_;
  assert((HighsInt)solution.col_value.size() ==
         mipsolver.orig_model_->num_col_);
  for (HighsInt i = 0; i != mipsolver.orig_model_->num_col_; ++i) {
    obj += mipsolver.orig_model_->col_cost_[i] * solution.col_value[i];

    bound_violation_ =
        std::max(bound_violation_,
                 mipsolver.orig_model_->col_lower_[i] - solution.col_value[i]);
    bound_violation_ =
        std::max(bound_violation_,
                 solution.col_value[i] - mipsolver.orig_model_->col_upper_[i]);

    if (mipsolver.orig_model_->integrality_[i] == HighsVarType::kInteger) {
      double intval = std::floor(solution.col_value[i] + 0.5);
      integrality_violation_ = std::max(
          std::abs(intval - solution.col_value[i]), integrality_violation_);
    }
  }

  for (HighsInt i = 0; i != mipsolver.orig_model_->num_row_; ++i) {
    row_violation_ =
        std::max(row_violation_,
                 mipsolver.orig_model_->row_lower_[i] - solution.row_value[i]);
    row_violation_ =
        std::max(row_violation_,
                 solution.row_value[i] - mipsolver.orig_model_->row_upper_[i]);
  }

  bool feasible =
      bound_violation_ <= mipsolver.options_mip_->mip_feasibility_tolerance &&
      integrality_violation_ <=
          mipsolver.options_mip_->mip_feasibility_tolerance &&
      row_violation_ <=
          mipsolver.options_mip_->mip_feasibility_tolerance + kHighsTiny;

  if (!feasible && allow_try_again) {
    // printf(
    //     "trying to repair sol that is violated by %.12g bounds, %.12g "
    //     "integrality, %.12g rows\n",
    //     bound_violation_, integrality_violation_, row_violation_);
    HighsLp fixedModel = *mipsolver.orig_model_;
    fixedModel.integrality_.clear();
    for (HighsInt i = 0; i != mipsolver.orig_model_->num_col_; ++i) {
      if (mipsolver.orig_model_->integrality_[i] == HighsVarType::kInteger) {
        double solval = std::round(solution.col_value[i]);
        fixedModel.col_lower_[i] = std::max(fixedModel.col_lower_[i], solval);
        fixedModel.col_upper_[i] = std::min(fixedModel.col_upper_[i], solval);
      }
    }
    Highs tmpSolver;
    tmpSolver.setOptionValue("output_flag", false);
    tmpSolver.setOptionValue("simplex_scale_strategy", 0);
    tmpSolver.setOptionValue("presolve", "off");
    tmpSolver.setOptionValue("primal_feasibility_tolerance",
                             mipsolver.options_mip_->mip_feasibility_tolerance);
    tmpSolver.passModel(std::move(fixedModel));
    tmpSolver.run();

    if (tmpSolver.getInfo().primal_solution_status == 2) {
      solution = tmpSolver.getSolution();
      allow_try_again = false;
      goto try_again;
    }
  }
  // store the solution as incumbent in the original space if there is no
  // solution or if it is feasible
  if (feasible) {
    // if (!allow_try_again)
    //   printf("repaired solution with value %g\n", double(obj));
    // store
    mipsolver.row_violation_ = row_violation_;
    mipsolver.bound_violation_ = bound_violation_;
    mipsolver.integrality_violation_ = integrality_violation_;
    mipsolver.solution_ = std::move(solution.col_value);
    mipsolver.solution_objective_ = double(obj);
  } else {
    bool currentFeasible =
        mipsolver.solution_objective_ != kHighsInf &&
        mipsolver.bound_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.integrality_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.row_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance + kHighsTiny;
    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kWarning,
        "Untransformed solution with objective %g is violated by %.12g for the "
        "original model\n",
        double(obj),
        std::max({bound_violation_, integrality_violation_, row_violation_}));
    if (!currentFeasible) {
      // if the current incumbent is non existent or also not feasible we still
      // store the new one
      mipsolver.row_violation_ = row_violation_;
      mipsolver.bound_violation_ = bound_violation_;
      mipsolver.integrality_violation_ = integrality_violation_;
      mipsolver.solution_ = std::move(solution.col_value);
      mipsolver.solution_objective_ = double(obj);
    }

    // return infinity so that it is not used for bounding
    return kHighsInf;
  }

  // return the objective value in the transformed space
  if (mipsolver.orig_model_->sense_ == ObjSense::kMaximize)
    return -double(obj + mipsolver.model_->offset_);

  return double(obj - mipsolver.model_->offset_);
}

double HighsMipSolverData::percentageInactiveIntegers() const {
  return 100.0 * (1.0 - double(integer_cols.size() -
                               cliquetable.getSubstitutions().size()) /
                            numintegercols);
}

void HighsMipSolverData::performRestart() {
  HighsBasis root_basis;
  HighsPseudocostInitialization pscostinit(
      pseudocost, mipsolver.options_mip_->mip_pscost_minreliable,
      postSolveStack);

  mipsolver.pscostinit = &pscostinit;
  ++numRestarts;
  num_leaves_before_run = num_leaves;
  num_nodes_before_run = num_nodes;
  num_nodes_before_run = num_nodes;
  total_lp_iterations_before_run = total_lp_iterations;
  heuristic_lp_iterations_before_run = heuristic_lp_iterations;
  sepa_lp_iterations_before_run = sepa_lp_iterations;
  sb_lp_iterations_before_run = sb_lp_iterations;
  HighsInt numLpRows = lp.getLp().num_row_;
  HighsInt numModelRows = mipsolver.numRow();
  HighsInt numCuts = numLpRows - numModelRows;
  if (numCuts > 0) postSolveStack.appendCutsToModel(numCuts);
  auto integrality = std::move(presolvedModel.integrality_);
  double offset = presolvedModel.offset_;
  presolvedModel = lp.getLp();
  presolvedModel.offset_ = offset;
  presolvedModel.integrality_ = std::move(integrality);
  const HighsBasis& basis = lp.getLpSolver().getBasis();
  if (basis.valid) {
    // if we have a basis after solving the root LP, we expand it to the
    // original space so that it can be used for constructing a starting basis
    // for the presolved model after the restart
    root_basis.col_status.resize(postSolveStack.getOrigNumCol());
    root_basis.row_status.resize(postSolveStack.getOrigNumRow());
    root_basis.valid = true;

    for (HighsInt i = 0; i != mipsolver.model_->num_col_; ++i)
      root_basis.col_status[postSolveStack.getOrigColIndex(i)] =
          basis.col_status[i];

    for (HighsInt i = 0; i != mipsolver.model_->num_row_; ++i)
      root_basis.row_status[postSolveStack.getOrigRowIndex(i)] =
          basis.row_status[i];

    mipsolver.rootbasis = &root_basis;
  }

  // transform the objective upper bound into the original space, as it is
  // expected during presolve
  upper_limit += mipsolver.model_->offset_;
  upper_bound += mipsolver.model_->offset_;
  lower_bound += mipsolver.model_->offset_;

  // remove the current incumbent. Any incumbent is already transformed into the
  // original space and kept there
  incumbent.clear();
  pruned_treeweight = 0;
  nodequeue.clear();
  globalOrbits.reset();

  runPresolve();

  if (mipsolver.modelstatus_ != HighsModelStatus::kNotset) {
    if (mipsolver.solution_objective_ != kHighsInf &&
        mipsolver.modelstatus_ == HighsModelStatus::kInfeasible)
      mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    // transform the objective limit to the current model
    upper_limit -= mipsolver.model_->offset_;
    upper_bound -= mipsolver.model_->offset_;
    lower_bound = upper_bound;
    return;
  }
  runSetup();

  postSolveStack.removeCutsFromModel(numCuts);

  // HighsNodeQueue oldNodeQueue;
  // std::swap(nodequeue, oldNodeQueue);

  // remove the pointer into the stack-space of this function
  if (mipsolver.rootbasis == &root_basis) mipsolver.rootbasis = nullptr;
  mipsolver.pscostinit = nullptr;
}

void HighsMipSolverData::basisTransfer() {
  // if a root basis is given, construct a basis for the root LP from
  // in the reduced problem space after presolving
  if (mipsolver.rootbasis) {
    HighsInt numRow = mipsolver.numRow() + cutpool.getNumCuts();
    firstrootbasis.col_status.assign(mipsolver.numCol(),
                                     HighsBasisStatus::kNonbasic);
    firstrootbasis.row_status.assign(numRow, HighsBasisStatus::kNonbasic);
    firstrootbasis.valid = true;
    HighsInt missingbasic = numRow;

    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      HighsBasisStatus status =
          mipsolver.rootbasis->col_status[postSolveStack.getOrigColIndex(i)];

      if (status == HighsBasisStatus::kBasic) {
        --missingbasic;
        firstrootbasis.col_status[i] = status;

        if (missingbasic == 0) break;
      }
    }

    if (missingbasic != 0) {
      for (HighsInt i = 0; i != numRow; ++i) {
        HighsBasisStatus status =
            mipsolver.rootbasis->row_status[postSolveStack.getOrigRowIndex(i)];

        if (status == HighsBasisStatus::kBasic) {
          --missingbasic;
          firstrootbasis.row_status[i] = status;
          if (missingbasic == 0) break;
        }
      }
    }

    const HighsLp& model = *mipsolver.model_;

    // there are missing basic variables; first add the sparsest nonbasic
    // structural columns to the basis whenever the column does not contain
    // any basic row. Then proceed by adding logical columns of rows which
    // contain no basic variables until the basis is complete
    if (missingbasic != 0) {
      std::vector<HighsInt> nonbasiccols;
      nonbasiccols.reserve(model.num_col_);
      for (HighsInt i = 0; i != model.num_col_; ++i) {
        if (firstrootbasis.col_status[i] != HighsBasisStatus::kBasic)
          nonbasiccols.push_back(i);
      }
      pdqsort(nonbasiccols.begin(), nonbasiccols.end(),
              [&](HighsInt col1, HighsInt col2) {
                HighsInt len1 = model.a_matrix_.start_[col1 + 1] -
                                model.a_matrix_.start_[col1];
                HighsInt len2 = model.a_matrix_.start_[col2 + 1] -
                                model.a_matrix_.start_[col2];
                return std::make_pair(len1, col1) < std::make_pair(len2, col2);
              });
      nonbasiccols.resize(std::min(nonbasiccols.size(), size_t(missingbasic)));
      for (HighsInt i : nonbasiccols) {
        const HighsInt start = model.a_matrix_.start_[i];
        const HighsInt end = model.a_matrix_.start_[i + 1];

        bool hasbasic = false;
        for (HighsInt j = start; j != end; ++j) {
          if (firstrootbasis.row_status[model.a_matrix_.index_[j]] ==
              HighsBasisStatus::kBasic) {
            hasbasic = true;
            break;
          }
        }

        if (!hasbasic) {
          firstrootbasis.col_status[i] = HighsBasisStatus::kBasic;
          --missingbasic;
          if (missingbasic == 0) break;
        }
      }

      if (missingbasic != 0) {
        std::vector<std::pair<HighsInt, int>> nonbasicrows;

        for (HighsInt i = 0; i != model.num_row_; ++i) {
          if (firstrootbasis.row_status[i] == HighsBasisStatus::kBasic)
            continue;

          const HighsInt start = ARstart_[i];
          const HighsInt end = ARstart_[i + 1];

          HighsInt nbasic = 0;
          for (HighsInt j = start; j != end; ++j) {
            if (firstrootbasis.col_status[ARindex_[j]] ==
                HighsBasisStatus::kBasic) {
              ++nbasic;
            }
          }

          if (nbasic == 0) {
            firstrootbasis.row_status[i] = HighsBasisStatus::kBasic;
            --missingbasic;
            if (missingbasic == 0) break;
          } else {
            nonbasicrows.emplace_back(nbasic, i);
          }
        }

        pdqsort(nonbasicrows.begin(), nonbasicrows.end());
        nonbasicrows.resize(missingbasic);

        for (std::pair<HighsInt, int> nonbasicrow : nonbasicrows)
          firstrootbasis.row_status[nonbasicrow.second] =
              HighsBasisStatus::kBasic;
      }
    }
  }
}

const std::vector<double>& HighsMipSolverData::getSolution() const {
  return incumbent;
}

bool HighsMipSolverData::addIncumbent(const std::vector<double>& sol,
                                      double solobj, char source) {
  if (solobj < upper_bound) {
    if (solobj <= upper_limit) {
      solobj = transformNewIncumbent(sol);
      if (solobj >= upper_bound) return false;
    }
    upper_bound = solobj;
    incumbent = sol;
    double new_upper_limit;
    if (objintscale != 0.0) {
      new_upper_limit =
          (std::floor(objintscale * solobj - 0.5) / objintscale) + feastol;
    } else {
      new_upper_limit = solobj - feastol;
    }
    if (new_upper_limit < upper_limit) {
      ++numImprovingSols;
      upper_limit = new_upper_limit;
      debugSolution.newIncumbentFound();
      redcostfixing.propagateRootRedcost(mipsolver);
      if (domain.infeasible()) {
        pruned_treeweight = 1.0;
        nodequeue.clear();
        return true;
      }
      cliquetable.extractObjCliques(mipsolver);
      if (domain.infeasible()) {
        pruned_treeweight = 1.0;
        nodequeue.clear();
        return true;
      }
      pruned_treeweight += nodequeue.performBounding(upper_limit);
      printDisplayLine(source);
    }
  } else if (incumbent.empty())
    incumbent = sol;

  return true;
}

static std::array<char, 16> convertToPrintString(int64_t val) {
  double l = std::log10(std::max(1.0, double(val)));
  std::array<char, 16> printString;
  switch (int(l)) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      std::snprintf(printString.data(), 16, "%ld", val);
      break;
    case 6:
    case 7:
    case 8:
      std::snprintf(printString.data(), 16, "%ldk", val / 1000);
      break;
    default:
      std::snprintf(printString.data(), 16, "%ldm", val / 1000000);
  }

  return printString;
}

static std::array<char, 32> convertToPrintString(double val) {
  std::array<char, 32> printString;
  double l = std::abs(val) == kHighsInf ? 0.0 : std::log10(std::max(1e-6, val));
  switch (int(l)) {
    case 0:
    case 1:
    case 2:
    case 3:
      std::snprintf(printString.data(), 32, "%.10g", val);
      break;
    case 4:
      std::snprintf(printString.data(), 32, "%.11g", val);
      break;
    case 5:
      std::snprintf(printString.data(), 32, "%.12g", val);
      break;
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
      std::snprintf(printString.data(), 32, "%.13g", val);
      break;
    default:
      std::snprintf(printString.data(), 32, "%.9g", val);
  }

  return printString;
}

void HighsMipSolverData::printDisplayLine(char first) {
  double time = mipsolver.timer_.read(mipsolver.timer_.solve_clock);
  if (first == ' ' && time - last_disptime < 5.) return;

  last_disptime = time;

  double offset = mipsolver.model_->offset_;
  if (num_disp_lines % 20 == 0) {
    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
        "\n        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints |       Work      \n"
          "     Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | LpIters     Time\n\n"
        // clang-format on
    );

    //"   %7s | %10s | %10s | %10s | %10s | %-15s | %-15s | %7s | %7s "
    //"| %8s | %8s\n",
    //"time", "open nodes", "nodes", "leaves", "lpiters", "dual bound",
    //"primal bound", "cutpool", "confl.", "gap", "explored");
  }

  ++num_disp_lines;

  std::array<char, 16> print_nodes = convertToPrintString(num_nodes);
  std::array<char, 16> queue_nodes = convertToPrintString(nodequeue.numNodes());
  std::array<char, 16> print_leaves =
      convertToPrintString(num_leaves - num_leaves_before_run);

  double explored = 100 * double(pruned_treeweight);

  double lb = lower_bound + offset;
  if (std::abs(lb) <= epsilon) lb = 0;
  double ub = kHighsInf;
  double gap = kHighsInf;

  std::array<char, 16> print_lp_iters =
      convertToPrintString(total_lp_iterations);
  if (upper_bound != kHighsInf) {
    ub = upper_bound + offset;
    if (std::abs(ub) <= epsilon) ub = 0;
    lb = std::min(ub, lb);
    gap = std::min(9999., 100 * (ub - lb) / std::max(1.0, std::abs(ub)));

    std::array<char, 32> lb_string = convertToPrintString(lb);
    std::array<char, 32> ub_string = convertToPrintString(ub);

    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
                 " %c %7s %7s   %7s %6.2f%%   %-15s %-15s %7.2f%%   %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT "   %7s %7.1fs\n",
        // clang-format on
        first, print_nodes.data(), queue_nodes.data(), print_leaves.data(),
        explored, lb_string.data(), ub_string.data(), gap, cutpool.getNumCuts(),
        lp.numRows() - lp.getNumModelRows(), conflictPool.getNumConflicts(),
        print_lp_iters.data(), time);
  } else {
    std::array<char, 32> lb_string = convertToPrintString(lb);
    std::array<char, 32> ub_string = convertToPrintString(ub);

    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::kInfo,
        // clang-format off
        " %c %7s %7s   %7s %6.2f%%   %-15s %-15s %8.2f   %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT " %6" HIGHSINT_FORMAT "   %7s %7.1fs\n",
        // clang-format on
        first, print_nodes.data(), queue_nodes.data(), print_leaves.data(),
        explored, lb_string.data(), ub_string.data(), gap, cutpool.getNumCuts(),
        lp.numRows() - lp.getNumModelRows(), conflictPool.getNumConflicts(),
        print_lp_iters.data(), time);
  }
}

bool HighsMipSolverData::rootSeparationRound(
    HighsSeparation& sepa, HighsInt& ncuts, HighsLpRelaxation::Status& status) {
  int64_t tmpLpIters = -lp.getNumLpIterations();
  ncuts = sepa.separationRound(domain, status);
  tmpLpIters += lp.getNumLpIterations();
  avgrootlpiters = lp.getAvgSolveIters();
  total_lp_iterations += tmpLpIters;
  sepa_lp_iterations += tmpLpIters;

  status = evaluateRootLp();
  if (status == HighsLpRelaxation::Status::kInfeasible) return true;

  const std::vector<double>& solvals = lp.getLpSolver().getSolution().col_value;

  if (mipsolver.submip || incumbent.empty()) {
    heuristics.randomizedRounding(solvals);
    heuristics.flushStatistics();
    status = evaluateRootLp();
    if (status == HighsLpRelaxation::Status::kInfeasible) return true;
  }

  return false;
}

HighsLpRelaxation::Status HighsMipSolverData::evaluateRootLp() {
  do {
    domain.propagate();

    if (globalOrbits && !domain.infeasible())
      globalOrbits->orbitalFixing(domain);

    if (domain.infeasible()) {
      lower_bound = std::min(kHighsInf, upper_bound);
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return HighsLpRelaxation::Status::kInfeasible;
    }

    bool lpBoundsChanged = false;
    if (!domain.getChangedCols().empty()) {
      lpBoundsChanged = true;
      removeFixedIndices();
      lp.flushDomain(domain);
    }

    bool lpWasSolved = false;
    HighsLpRelaxation::Status status;
    if (lpBoundsChanged ||
        lp.getLpSolver().getModelStatus(true) == HighsModelStatus::kNotset) {
      int64_t lpIters = -lp.getNumLpIterations();
      status = lp.resolveLp(&domain);
      lpIters += lp.getNumLpIterations();
      total_lp_iterations += lpIters;
      avgrootlpiters = lp.getAvgSolveIters();
      lpWasSolved = true;
      if (status == HighsLpRelaxation::Status::kOptimal &&
          lp.getFractionalIntegers().empty() &&
          addIncumbent(lp.getLpSolver().getSolution().col_value,
                       lp.getObjective(), 'T')) {
        mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
        lower_bound = upper_bound;
        pruned_treeweight = 1.0;
        num_nodes += 1;
        num_leaves += 1;
        return HighsLpRelaxation::Status::kInfeasible;
      }
    } else
      status = lp.getStatus();

    if (status == HighsLpRelaxation::Status::kInfeasible) {
      lower_bound = std::min(kHighsInf, upper_bound);
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return status;
    }

    if (lp.unscaledDualFeasible(lp.getStatus())) {
      lower_bound = std::max(lp.getObjective(), lower_bound);
      if (lpWasSolved) {
        redcostfixing.addRootRedcost(mipsolver,
                                     lp.getLpSolver().getSolution().col_dual,
                                     lp.getObjective());
        if (upper_limit != kHighsInf)
          redcostfixing.propagateRootRedcost(mipsolver);
      }
    }

    if (lower_bound > upper_limit) {
      lower_bound = std::min(kHighsInf, upper_bound);
      pruned_treeweight = 1.0;
      num_nodes += 1;
      num_leaves += 1;
      return HighsLpRelaxation::Status::kInfeasible;
    }

    if (domain.getChangedCols().empty()) return status;
  } while (true);
}

void HighsMipSolverData::evaluateRootNode() {
  HighsInt maxSepaRounds = mipsolver.submip ? 5 : kHighsIInf;
restart:
  // lp.getLpSolver().setOptionValue(
  //     "dual_simplex_cost_perturbation_multiplier", 10.0);
  lp.setIterationLimit();
  lp.loadModel();
  domain.clearChangedCols();
  lp.setObjectiveLimit(upper_limit);

  if (symmetries.numPerms != 0)
    globalOrbits = symmetries.computeStabilizerOrbits(domain);

  // add all cuts again after restart
  if (cutpool.getNumCuts() != 0) {
    highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                "\nAdding %" HIGHSINT_FORMAT
                " cuts to the LP after performing a restart\n",
                cutpool.getNumCuts());
    assert(numRestarts != 0);
    HighsCutSet cutset;
    cutpool.separateLpCutsAfterRestart(cutset);
#ifdef HIGHS_DEBUGSOL
    for (HighsInt i = 0; i < cutset.numCuts(); ++i) {
      debugSolution.checkCut(cutset.ARindex_.data() + cutset.ARstart_[i],
                             cutset.ARvalue_.data() + cutset.ARstart_[i],
                             cutset.ARstart_[i + 1] - cutset.ARstart_[i],
                             cutset.upper_[i]);
    }
#endif
    lp.addCuts(cutset);
    // solve the first root lp
    highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                "Solving root node LP relaxation\n");
  } else if (numRestarts == 0) {
    // solve the first root lp
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                 "\nSolving root node LP relaxation\n");
  }

  if (firstrootbasis.valid)
    lp.getLpSolver().setBasis(firstrootbasis);
  else
    lp.getLpSolver().setOptionValue("presolve", "on");
  if (mipsolver.options_mip_->highs_debug_level)
    lp.getLpSolver().setOptionValue("output_flag",
                                    mipsolver.options_mip_->output_flag);
  //  lp.getLpSolver().setOptionValue("log_dev_level", kHighsLogDevLevelInfo);
  //  lp.getLpSolver().setOptionValue("log_file",
  //  mipsolver.options_mip_->log_file);
  HighsLpRelaxation::Status status = evaluateRootLp();
  if (numRestarts == 0) firstrootlpiters = total_lp_iterations;

  lp.getLpSolver().setOptionValue("output_flag", false);
  lp.getLpSolver().setOptionValue("presolve", "off");
  lp.setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));
  lp.getLpSolver().setOptionValue("parallel", "off");

  if (status == HighsLpRelaxation::Status::kInfeasible) return;

  heuristics.randomizedRounding(firstlpsol);
  heuristics.flushStatistics();

  status = evaluateRootLp();
  if (status == HighsLpRelaxation::Status::kInfeasible) return;

  firstlpsol = lp.getSolution().col_value;
  firstlpsolobj = lp.getObjective();

  if (lp.getLpSolver().getBasis().valid && lp.numRows() == mipsolver.numRow())
    firstrootbasis = lp.getLpSolver().getBasis();
  else {
    // the root basis is later expected to be consistent for the model without
    // cuts so set it to the slack basis if the current basis already includes
    // cuts, e.g. due to a restart
    firstrootbasis.col_status.assign(mipsolver.numCol(),
                                     HighsBasisStatus::kNonbasic);
    firstrootbasis.row_status.assign(mipsolver.numRow(),
                                     HighsBasisStatus::kBasic);
    firstrootbasis.valid = true;
  }
  rootlpsolobj = firstlpsolobj;

  // begin separation
  std::vector<double> avgdirection;
  std::vector<double> curdirection;
  avgdirection.resize(mipsolver.numCol());
  curdirection.resize(mipsolver.numCol());

  HighsInt stall = 0;
  double smoothprogress = 0.0;
  HighsInt nseparounds = 0;
  HighsSeparation sepa(mipsolver);
  sepa.setLpRelaxation(&lp);

  while (lp.scaledOptimal(status) && !lp.getFractionalIntegers().empty() &&
         stall < 3) {
    printDisplayLine();

    if (checkLimits()) return;

    if (nseparounds == maxSepaRounds) break;

    removeFixedIndices();

    if (!mipsolver.submip &&
        mipsolver.options_mip_->presolve != kHighsOffString) {
      double fixingRate = percentageInactiveIntegers();
      if (fixingRate >= 10.0) {
        stall = -1;
        break;
      }
    }

    ++nseparounds;

    HighsInt ncuts;
    if (rootSeparationRound(sepa, ncuts, status)) return;
    if (nseparounds >= 5 && !mipsolver.submip && !analyticCenterComputed) {
      analyticCenterComputed = true;
      heuristics.centralRounding();
      heuristics.flushStatistics();

      status = evaluateRootLp();
      if (status == HighsLpRelaxation::Status::kInfeasible) return;
    }

    HighsCDouble sqrnorm = 0.0;
    const auto& solvals = lp.getSolution().col_value;

    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      curdirection[i] = firstlpsol[i] - solvals[i];

      // if (mip.integrality_[i] == 2 && lp.getObjective() > firstobj &&
      //    std::abs(curdirection[i]) > 1e-6)
      //  pseudocost.addObservation(i, -curdirection[i],
      //                            lp.getObjective() - firstobj);

      sqrnorm += curdirection[i] * curdirection[i];
    }
#if 1
    double scale = double(1.0 / sqrt(sqrnorm));
    sqrnorm = 0.0;
    HighsCDouble dotproduct = 0.0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      avgdirection[i] =
          (scale * curdirection[i] - avgdirection[i]) / nseparounds;
      sqrnorm += avgdirection[i] * avgdirection[i];
      dotproduct += avgdirection[i] * curdirection[i];
    }
#endif

    double progress = double(dotproduct / sqrt(sqrnorm));

    if (nseparounds == 1) {
      smoothprogress = progress;
    } else {
      double alpha = 1.0 / 3.0;
      double nextprogress = (1.0 - alpha) * smoothprogress + alpha * progress;

      if (nextprogress < smoothprogress * 1.01 &&
          (lp.getObjective() - firstlpsolobj) <=
              (rootlpsolobj - firstlpsolobj) * 1.001)
        ++stall;
      else {
        stall = 0;
      }
      smoothprogress = nextprogress;
    }

    rootlpsolobj = lp.getObjective();
    lp.setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));
    if (ncuts == 0) break;
  }

  lp.setIterationLimit();
  status = evaluateRootLp();
  if (status == HighsLpRelaxation::Status::kInfeasible) return;

  rootlpsol = lp.getLpSolver().getSolution().col_value;
  rootlpsolobj = lp.getObjective();
  lp.setIterationLimit(std::max(10000, int(10 * avgrootlpiters)));

  if (!analyticCenterComputed || upper_limit == kHighsInf) {
    analyticCenterComputed = true;
    heuristics.centralRounding();
    heuristics.flushStatistics();

    // if there are new global bound changes we reevaluate the LP and do one
    // more separation round
    bool separate = !domain.getChangedCols().empty();
    status = evaluateRootLp();
    if (status == HighsLpRelaxation::Status::kInfeasible) return;
    if (separate && lp.scaledOptimal(status)) {
      HighsInt ncuts;
      if (rootSeparationRound(sepa, ncuts, status)) return;
      ++nseparounds;
      printDisplayLine();
    }
  }

  printDisplayLine();

  do {
    if (rootlpsol.empty()) break;
    if (upper_limit != kHighsInf && !moreHeuristicsAllowed()) break;

    double oldLimit = upper_limit;
    heuristics.rootReducedCost();
    heuristics.flushStatistics();

    // if there are new global bound changes we reevaluate the LP and do one
    // more separation round
    bool separate = !domain.getChangedCols().empty();
    status = evaluateRootLp();
    if (status == HighsLpRelaxation::Status::kInfeasible) return;
    if (separate && lp.scaledOptimal(status)) {
      HighsInt ncuts;
      if (rootSeparationRound(sepa, ncuts, status)) return;

      ++nseparounds;
      printDisplayLine();
    }

    if (upper_limit != kHighsInf && !moreHeuristicsAllowed()) break;

    heuristics.RENS(rootlpsol);
    heuristics.flushStatistics();

    // if there are new global bound changes we reevaluate the LP and do one
    // more separation round
    separate = !domain.getChangedCols().empty();
    status = evaluateRootLp();
    if (status == HighsLpRelaxation::Status::kInfeasible) return;
    if (separate && lp.scaledOptimal(status)) {
      HighsInt ncuts;
      if (rootSeparationRound(sepa, ncuts, status)) return;

      ++nseparounds;

      printDisplayLine();
    }

    if (upper_limit != kHighsInf || mipsolver.submip) break;
    heuristics.feasibilityPump();
    heuristics.flushStatistics();

    status = evaluateRootLp();
    if (status == HighsLpRelaxation::Status::kInfeasible) return;
  } while (false);

  if (lower_bound > upper_limit) {
    mipsolver.modelstatus_ = HighsModelStatus::kOptimal;
    pruned_treeweight = 1.0;
    num_nodes += 1;
    num_leaves += 1;
    return;
  }

  // if there are new global bound changes we reevaluate the LP and do one
  // more separation round
  bool separate = !domain.getChangedCols().empty();
  status = evaluateRootLp();
  if (status == HighsLpRelaxation::Status::kInfeasible) return;
  if (separate && lp.scaledOptimal(status)) {
    HighsInt ncuts;
    if (rootSeparationRound(sepa, ncuts, status)) return;

    ++nseparounds;
    printDisplayLine();
  }

  removeFixedIndices();
  if (lp.getLpSolver().getBasis().valid) lp.removeObsoleteRows();
  rootlpsolobj = lp.getObjective();

  printDisplayLine();

  if (lower_bound <= upper_limit) {
    if (!mipsolver.submip &&
        mipsolver.options_mip_->presolve != kHighsOffString) {
      double fixingRate = percentageInactiveIntegers();
      if (fixingRate >= 2.5 + 7.5 * mipsolver.submip ||
          (!mipsolver.submip && fixingRate > 0 && numRestarts == 0)) {
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                     "\n%.1f%% inactive integer columns, restarting\n",
                     fixingRate);
        if (stall != -1) maxSepaRounds = std::min(maxSepaRounds, nseparounds);
        performRestart();
        ++numRestartsRoot;
        if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) goto restart;

        return;
      }
    }
    // add the root node to the nodequeue to initialize the search
    nodequeue.emplaceNode(std::vector<HighsDomainChange>(),
                          std::vector<HighsInt>(), lower_bound,
                          lp.computeBestEstimate(pseudocost), 1);
  }
}

bool HighsMipSolverData::checkLimits() const {
  const HighsOptions& options = *mipsolver.options_mip_;
  if (options.mip_max_nodes != kHighsIInf &&
      num_nodes >= options.mip_max_nodes) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "reached node limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kIterationLimit;
    }
    return true;
  }
  if (options.mip_max_leaves != kHighsIInf &&
      num_leaves >= options.mip_max_leaves) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "reached leave node limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kIterationLimit;
    }
    return true;
  }
  if (mipsolver.timer_.read(mipsolver.timer_.solve_clock) >=
      options.time_limit) {
    if (mipsolver.modelstatus_ == HighsModelStatus::kNotset) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "reached time limit\n");
      mipsolver.modelstatus_ = HighsModelStatus::kTimeLimit;
    }
    return true;
  }

  return false;
}

void HighsMipSolverData::checkObjIntegrality() {
  objintscale = 600.0;

  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.colCost(i) == 0.0) continue;

    if (mipsolver.variableType(i) == HighsVarType::kContinuous) {
      objintscale = 0.0;
      break;
    }

    double cost = mipsolver.colCost(i);
    double intcost = std::floor(objintscale * cost + 0.5) / objintscale;
    if (std::abs(cost - intcost) > epsilon) {
      objintscale = 0.0;
      break;
    }
  }

  if (objintscale != 0.0) {
    int64_t currgcd = 0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.colCost(i) == 0.0) continue;
      int64_t intval = std::floor(mipsolver.colCost(i) * objintscale + 0.5);
      if (currgcd == 0) {
        currgcd = intval < 0 ? -intval : intval;
        continue;
      }
      currgcd = HighsIntegers::gcd(intval, currgcd);
      if (currgcd == 1) break;
    }

    if (currgcd != 0) objintscale /= currgcd;

    if (numRestarts == 0)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "Objective function is integral with scale %g\n",
                   objintscale);
  }
}

void HighsMipSolverData::setupDomainPropagation() {
  const HighsLp& model = *mipsolver.model_;
  highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                       model.a_matrix_.index_, model.a_matrix_.value_, ARstart_,
                       ARindex_, ARvalue_);

  pseudocost = HighsPseudocost(mipsolver);

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.model_->num_row_);
  for (HighsInt i = 0; i != mipsolver.model_->num_row_; ++i) {
    double maxabsval = 0.0;

    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];
    for (HighsInt j = start; j != end; ++j)
      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));

    maxAbsRowCoef[i] = maxabsval;
  }

  domain = HighsDomain(mipsolver);
  domain.computeRowActivities();
}
