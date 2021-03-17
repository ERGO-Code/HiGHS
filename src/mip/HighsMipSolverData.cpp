/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolverData.h"

#include <random>

#include "lp_data/HighsLpUtils.h"
#include "presolve/HAggregator.h"
#include "presolve/HPresolve.h"
#include "util/HighsIntegers.h"

bool HighsMipSolverData::trySolution(const std::vector<double>& solution,
                                     char source) {
  if (int(solution.size()) != mipsolver.model_->numCol_) return false;

  HighsCDouble obj = 0;

  for (int i = 0; i != mipsolver.model_->numCol_; ++i) {
    if (solution[i] < mipsolver.model_->colLower_[i] - feastol) return false;
    if (solution[i] > mipsolver.model_->colUpper_[i] + feastol) return false;
    if (mipsolver.variableType(i) == HighsVarType::INTEGER &&
        std::abs(solution[i] - std::floor(solution[i] + 0.5)) > feastol)
      return false;

    obj += mipsolver.colCost(i) * solution[i];
  }

  for (int i = 0; i != mipsolver.model_->numRow_; ++i) {
    double rowactivity = 0.0;

    int start = ARstart_[i];
    int end = ARstart_[i + 1];

    for (int j = start; j != end; ++j)
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
  } else if (pruned_treeweight < 1e-3 && num_leaves < 10) {
    // in the main MIP solver allow an initial offset of 10000 heuristic LP
    // iterations
    if (heuristic_lp_iterations <
        total_lp_iterations * heuristic_effort + 10000)
      return true;
  } else {
    double total_heuristic_effort_estim =
        heuristic_lp_iterations /
        (heuristic_lp_iterations + sb_lp_iterations +
         (total_lp_iterations - heuristic_lp_iterations - sb_lp_iterations) /
             std::max(1e-3, double(pruned_treeweight)));
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
        std::max(0.3 / 0.8, double(pruned_treeweight) / 0.8) *
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
                     [&](int col) { return domain.isFixed(col); }),
      integral_cols.end());
  integer_cols.erase(
      std::remove_if(integer_cols.begin(), integer_cols.end(),
                     [&](int col) { return domain.isFixed(col); }),
      integer_cols.end());
  implint_cols.erase(
      std::remove_if(implint_cols.begin(), implint_cols.end(),
                     [&](int col) { return domain.isFixed(col); }),
      implint_cols.end());
  continuous_cols.erase(
      std::remove_if(continuous_cols.begin(), continuous_cols.end(),
                     [&](int col) { return domain.isFixed(col); }),
      continuous_cols.end());
}

void HighsMipSolverData::init() {
  postSolveStack.initializeIndexMaps(mipsolver.model_->numRow_,
                                     mipsolver.model_->numCol_);
  mipsolver.orig_model_ = mipsolver.model_;
  feastol = mipsolver.options_mip_->mip_feasibility_tolerance;
  epsilon = mipsolver.options_mip_->mip_epsilon;
  heuristic_effort = mipsolver.options_mip_->mip_heuristic_effort;

  firstlpsolobj = -HIGHS_CONST_INF;
  rootlpsolobj = -HIGHS_CONST_INF;

  numRestarts = 0;
  numImprovingSols = 0;
  pruned_treeweight = 0;
  maxrootlpiters = 0;
  num_nodes = 0;
  num_leaves = 0;
  total_lp_iterations = 0;
  heuristic_lp_iterations = 0;
  sepa_lp_iterations = 0;
  sb_lp_iterations = 0;
  num_disp_lines = 0;
  last_displeave = 0;
  cliquesExtracted = false;
  rowMatrixSet = false;
  lower_bound = -HIGHS_CONST_INF;
  upper_bound = HIGHS_CONST_INF;
  upper_limit = mipsolver.options_mip_->dual_objective_value_upper_bound;

  if (mipsolver.options_mip_->mip_report_level == 0)
    dispfreq = 0;
  else if (mipsolver.options_mip_->mip_report_level == 1)
    dispfreq = 100;
  else
    dispfreq = 1;
}

void HighsMipSolverData::runPresolve() {
  presolve::HPresolve presolve;

  presolve.setInput(mipsolver);

  mipsolver.modelstatus_ = presolve.run(postSolveStack);
}

void HighsMipSolverData::runSetup() {
  const HighsLp& model = *mipsolver.model_;

  // transform the objective limit to the current model
  upper_limit -= mipsolver.model_->offset_;
  lower_bound -= mipsolver.model_->offset_;
  upper_bound -= mipsolver.model_->offset_;

  redcostfixing = HighsRedcostFixing();
  mipsolver.mipdata_->pseudocost = HighsPseudocost(mipsolver.model_->numCol_);
  continuous_cols.clear();
  integer_cols.clear();
  implint_cols.clear();
  integral_cols.clear();
  if (mipsolver.submip) pseudocost.setMinReliable(0);

  rowMatrixSet = false;
  if (!rowMatrixSet) {
    rowMatrixSet = true;
    highsSparseTranspose(model.numRow_, model.numCol_, model.Astart_,
                         model.Aindex_, model.Avalue_, ARstart_, ARindex_,
                         ARvalue_);
    uplocks.resize(model.numCol_);
    downlocks.resize(model.numCol_);
    for (int i = 0; i != model.numCol_; ++i) {
      int start = model.Astart_[i];
      int end = model.Astart_[i + 1];
      for (int j = start; j != end; ++j) {
        int row = model.Aindex_[j];

        if (model.rowLower_[row] != -HIGHS_CONST_INF) {
          if (model.Avalue_[j] < 0)
            ++uplocks[i];
          else
            ++downlocks[i];
        }
        if (model.rowUpper_[row] != HIGHS_CONST_INF) {
          if (model.Avalue_[j] < 0)
            ++downlocks[i];
          else
            ++uplocks[i];
        }
      }
    }
  }

  rowintegral.resize(mipsolver.model_->numRow_);

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.model_->numRow_);
  for (int i = 0; i != mipsolver.model_->numRow_; ++i) {
    double maxabsval = 0.0;

    int start = ARstart_[i];
    int end = ARstart_[i + 1];
    bool integral = true;
    for (int j = start; j != end; ++j) {
      if (integral) {
        if (mipsolver.variableType(ARindex_[j]) == HighsVarType::CONTINUOUS)
          integral = false;
        else {
          double intval = std::floor(ARvalue_[j] + 0.5);
          if (std::abs(ARvalue_[j] - intval) > epsilon) integral = false;
        }
      }

      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));
    }

    rowintegral[i] = integral;

    maxAbsRowCoef[i] = maxabsval;
  }

  // compute row activities and propagate all rows once
  domain.computeRowActivities();
  domain.propagate();
  if (domain.infeasible()) {
    mipsolver.modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
    lower_bound = HIGHS_CONST_INF;
    pruned_treeweight = 1.0;
    return;
  }

  if (model.numCol_ == 0) {
    mipsolver.modelstatus_ = HighsModelStatus::OPTIMAL;
    return;
  }

  if (checkLimits()) return;
  // extract cliques if they have not been extracted before

  for (int col : domain.getChangedCols()) implications.cleanupVarbounds(col);
  domain.clearChangedCols();

  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  // lp.getLpSolver().setHighsOptionValue("dual_simplex_cleanup_strategy", 0);
  // lp.getLpSolver().setHighsOptionValue("dual_simplex_cost_perturbation_multiplier",
  // 0.0); lp.getLpSolver().setHighsOptionValue("parallel", "on");
  lp.getLpSolver().setHighsOptionValue("simplex_initial_condition_check",
                                       false);

  checkObjIntegrality();
  basisTransfer();

  for (int i = 0; i != mipsolver.numCol(); ++i) {
    switch (mipsolver.variableType(i)) {
      case HighsVarType::CONTINUOUS:
        continuous_cols.push_back(i);
        break;
      case HighsVarType::IMPLICIT_INTEGER:
        implint_cols.push_back(i);
        integral_cols.push_back(i);
        break;
      case HighsVarType::INTEGER:
        integer_cols.push_back(i);
        integral_cols.push_back(i);
    }
  }
  nodequeue.setNumCol(mipsolver.numCol());
  numintegercols = integer_cols.size();

  debugSolution.activate();
}

double HighsMipSolverData::transformNewIncumbent(
    const std::vector<double>& sol) {
  HighsSolution solution;
  solution.col_value = sol;
  calculateRowValues(*mipsolver.model_, solution);

  // todo, add interface to postSolveStack that does not expect a basis
  HighsBasis basis;

  postSolveStack.undo(*mipsolver.options_mip_, solution, basis);
  calculateRowValues(*mipsolver.orig_model_, solution);

  // compute the objective value in the original space
  double bound_violation_ = 0;
  double row_violation_ = 0;
  double integrality_violation_ = 0;

  HighsCDouble obj = mipsolver.orig_model_->offset_;
  for (int i = 0; i != mipsolver.orig_model_->numCol_; ++i) {
    obj += mipsolver.orig_model_->colCost_[i] * solution.col_value[i];

    bound_violation_ =
        std::max(bound_violation_,
                 mipsolver.orig_model_->colLower_[i] - solution.col_value[i]);
    bound_violation_ =
        std::max(bound_violation_,
                 solution.col_value[i] - mipsolver.orig_model_->colUpper_[i]);

    if (mipsolver.orig_model_->integrality_[i] == HighsVarType::INTEGER) {
      double intval = std::floor(solution.col_value[i] + 0.5);
      integrality_violation_ = std::max(
          std::abs(intval - solution.col_value[i]), integrality_violation_);
    }
  }

  for (int i = 0; i != mipsolver.orig_model_->numRow_; ++i) {
    row_violation_ =
        std::max(row_violation_,
                 mipsolver.orig_model_->rowLower_[i] - solution.row_value[i]);
    row_violation_ =
        std::max(row_violation_,
                 solution.row_value[i] - mipsolver.orig_model_->rowUpper_[i]);
  }

  bool feasible =
      bound_violation_ <= mipsolver.options_mip_->mip_feasibility_tolerance &&
      integrality_violation_ <=
          mipsolver.options_mip_->mip_feasibility_tolerance &&
      row_violation_ <= mipsolver.options_mip_->mip_feasibility_tolerance;
  assert(feasible);
  // store the solution as incumbent in the original space if there is no
  // solution or if it is feasible
  if (feasible) {
    // store
    mipsolver.row_violation_ = row_violation_;
    mipsolver.bound_violation_ = bound_violation_;
    mipsolver.integrality_violation_ = integrality_violation_;
    mipsolver.solution_ = std::move(solution.col_value);
    mipsolver.solution_objective_ = double(obj);
  } else {
    bool currentFeasible =
        mipsolver.solution_objective_ != HIGHS_CONST_INF &&
        mipsolver.bound_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.integrality_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance &&
        mipsolver.row_violation_ <=
            mipsolver.options_mip_->mip_feasibility_tolerance;
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
    return HIGHS_CONST_INF;
  }

  // return the objective value in the transformed space
  if (mipsolver.orig_model_->sense_ == ObjSense::MAXIMIZE)
    return -double(obj + mipsolver.model_->offset_);

  return double(obj - mipsolver.model_->offset_);
}

void HighsMipSolverData::performRestart() {
  HighsBasis rootBasis;
  ++numRestarts;
  int numLpRows = lp.getLp().numRow_;
  int numModelRows = mipsolver.numRow();
  int numCuts = numLpRows - numModelRows;
  if (numCuts > 0) postSolveStack.appendCutsToModel(numCuts);
  auto integrality = std::move(presolvedModel.integrality_);
  presolvedModel = lp.getLp();
  presolvedModel.integrality_ = std::move(integrality);
  const HighsBasis& basis = lp.getLpSolver().getBasis();
  if (basis.valid_) {
    // if we have a basis after solving the root LP, we expand it to the
    // original space so that it can be used for constructing a starting basis
    // for the presolved model after the restart
    rootBasis.col_status.resize(postSolveStack.getOrigNumCol());
    rootBasis.row_status.resize(postSolveStack.getOrigNumRow());
    rootBasis.valid_ = true;

    for (int i = 0; i != mipsolver.model_->numCol_; ++i)
      rootBasis.col_status[postSolveStack.getOrigColIndex(i)] =
          basis.col_status[i];

    for (int i = 0; i != mipsolver.model_->numRow_; ++i)
      rootBasis.row_status[postSolveStack.getOrigRowIndex(i)] =
          basis.row_status[i];

    mipsolver.rootbasis = &rootBasis;
  }

  // transform the objective upper bound into the original space, as it is
  // expected during presolve
  upper_limit += mipsolver.model_->offset_;
  upper_bound += mipsolver.model_->offset_;
  lower_bound += mipsolver.model_->offset_;

  // remove the current incumbent. Any incumbent is already transformed into the
  // original space and kept there
  incumbent.clear();

  runPresolve();

  if (mipsolver.modelstatus_ != HighsModelStatus::NOTSET) {
    if (mipsolver.solution_objective_ != HIGHS_CONST_INF &&
        mipsolver.modelstatus_ == HighsModelStatus::PRIMAL_INFEASIBLE)
      mipsolver.modelstatus_ = HighsModelStatus::OPTIMAL;
    return;
  }
  runSetup();

  pruned_treeweight = 0;
  nodequeue.clear();
  // HighsNodeQueue oldNodeQueue;
  // std::swap(nodequeue, oldNodeQueue);

  // remove the pointer into the stack-space of this function
  if (mipsolver.rootbasis == &rootBasis) mipsolver.rootbasis = nullptr;
}

void HighsMipSolverData::basisTransfer() {
  // if a root basis is given, construct a basis for the root LP from
  // in the reduced problem space after presolving
  if (mipsolver.rootbasis) {
    int numRow = mipsolver.numRow() + cutpool.getNumCuts();
    firstrootbasis.col_status.assign(mipsolver.numCol(),
                                     HighsBasisStatus::NONBASIC);
    firstrootbasis.row_status.assign(numRow, HighsBasisStatus::NONBASIC);
    firstrootbasis.valid_ = true;
    int missingbasic = numRow;

    for (int i = 0; i != mipsolver.numCol(); ++i) {
      HighsBasisStatus status =
          mipsolver.rootbasis->col_status[postSolveStack.getOrigColIndex(i)];

      if (status == HighsBasisStatus::BASIC) {
        --missingbasic;
        firstrootbasis.col_status[i] = status;

        if (missingbasic == 0) break;
      }
    }

    if (missingbasic != 0) {
      for (int i = 0; i != numRow; ++i) {
        HighsBasisStatus status =
            mipsolver.rootbasis->row_status[postSolveStack.getOrigRowIndex(i)];

        if (status == HighsBasisStatus::BASIC) {
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
      std::vector<int> nonbasiccols;
      nonbasiccols.reserve(model.numCol_);
      for (int i = 0; i != model.numCol_; ++i) {
        if (firstrootbasis.col_status[i] != HighsBasisStatus::BASIC)
          nonbasiccols.push_back(i);
      }
      std::sort(nonbasiccols.begin(), nonbasiccols.end(),
                [&](int col1, int col2) {
                  int len1 = model.Astart_[col1 + 1] - model.Astart_[col1];
                  int len2 = model.Astart_[col2 + 1] - model.Astart_[col2];
                  return len1 < len2;
                });
      nonbasiccols.resize(std::min(nonbasiccols.size(), size_t(missingbasic)));
      for (int i : nonbasiccols) {
        const int start = model.Astart_[i];
        const int end = model.Astart_[i + 1];

        bool hasbasic = false;
        for (int j = start; j != end; ++j) {
          if (firstrootbasis.row_status[model.Aindex_[j]] ==
              HighsBasisStatus::BASIC) {
            hasbasic = true;
            break;
          }
        }

        if (!hasbasic) {
          firstrootbasis.col_status[i] = HighsBasisStatus::BASIC;
          --missingbasic;
          if (missingbasic == 0) break;
        }
      }

      if (missingbasic != 0) {
        std::vector<std::pair<int, int>> nonbasicrows;

        for (int i = 0; i != model.numRow_; ++i) {
          if (firstrootbasis.row_status[i] == HighsBasisStatus::BASIC) continue;

          const int start = ARstart_[i];
          const int end = ARstart_[i + 1];

          int nbasic = 0;
          for (int j = start; j != end; ++j) {
            if (firstrootbasis.col_status[ARindex_[j]] ==
                HighsBasisStatus::BASIC) {
              ++nbasic;
            }
          }

          if (nbasic == 0) {
            firstrootbasis.row_status[i] = HighsBasisStatus::BASIC;
            --missingbasic;
            if (missingbasic == 0) break;
          } else {
            nonbasicrows.emplace_back(nbasic, i);
          }
        }

        std::sort(nonbasicrows.begin(), nonbasicrows.end());
        nonbasicrows.resize(missingbasic);

        for (std::pair<int, int> nonbasicrow : nonbasicrows)
          firstrootbasis.row_status[nonbasicrow.second] =
              HighsBasisStatus::BASIC;
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
      debugSolution.newIncumbentFound();
      upper_limit = new_upper_limit;
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
  }

  return true;
}

void HighsMipSolverData::printDisplayLine(char first) {
  double offset = mipsolver.model_->offset_;
  if (num_disp_lines % 20 == 0) {
    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::INFO,
        "   %7s | %10s | %10s | %10s | %10s | %-14s | %-14s | %7s | %7s "
        "| %8s | %8s\n",
        "time", "open nodes", "nodes", "leaves", "lpiters", "dual bound",
        "primal bound", "cutpool", "lpcuts", "gap", "explored");
  }

  ++num_disp_lines;
  last_displeave = num_leaves;

  double lb = mipsolver.mipdata_->lower_bound + offset;
  double ub = HIGHS_CONST_INF;
  double gap = HIGHS_CONST_INF;
  int lpcuts = mipsolver.mipdata_->lp.numRows() - mipsolver.model_->numRow_;

  if (upper_bound != HIGHS_CONST_INF) {
    ub = upper_bound + offset;
    lb = std::min(ub, lb);
    gap = 100 * (ub - lb) / std::max(1.0, std::abs(ub));

    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::INFO,
        " %c %6.1fs | %10lu | %10lu | %10lu | %10lu | %-14.9g | %-14.9g | "
        "%7d | %7d | %7.2f%% | %7.2f%%\n",
        first, mipsolver.timer_.read(mipsolver.timer_.solve_clock),
        nodequeue.numNodes(), num_nodes, num_leaves, total_lp_iterations, lb,
        ub, mipsolver.mipdata_->cutpool.getNumCuts(), lpcuts, gap,
        100 * double(pruned_treeweight));
  } else {
    highsLogUser(
        mipsolver.options_mip_->log_options, HighsLogType::INFO,
        " %c %6.1fs | %10lu | %10lu | %10lu | %10lu | %-14.9g | %-14.9g | "
        "%7d | %7d | %8.2f | %7.2f%%\n",
        first, mipsolver.timer_.read(mipsolver.timer_.solve_clock),
        nodequeue.numNodes(), num_nodes, num_leaves, total_lp_iterations, lb,
        ub, mipsolver.mipdata_->cutpool.getNumCuts(), lpcuts, gap,
        100 * double(pruned_treeweight));
  }
}

bool HighsMipSolverData::rootSeparationRound(
    HighsSeparation& sepa, int& ncuts, HighsLpRelaxation::Status& status) {
  size_t tmpLpIters = lp.getNumLpIterations();
  ncuts = sepa.separationRound(domain, status);
  maxrootlpiters =
      std::max(maxrootlpiters, lp.getNumLpIterations() - tmpLpIters);

  total_lp_iterations = lp.getNumLpIterations();
  sepa_lp_iterations = total_lp_iterations - firstrootlpiters;

  if (status == HighsLpRelaxation::Status::Infeasible) {
    pruned_treeweight = 1.0;
    lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
    num_nodes = 1;
    num_leaves = 1;
    return true;
  }

  const std::vector<double>& solvals = lp.getLpSolver().getSolution().col_value;

  if (incumbent.empty()) {
    heuristics.randomizedRounding(solvals);
    heuristics.flushStatistics();
  }

  if (lp.unscaledDualFeasible(status)) {
    lower_bound = lp.getObjective();
    redcostfixing.addRootRedcost(
        mipsolver, lp.getLpSolver().getSolution().col_dual, lower_bound);
    if (upper_limit != HIGHS_CONST_INF) {
      redcostfixing.propagateRootRedcost(mipsolver);

      if (domain.infeasible())
        status = HighsLpRelaxation::Status::Infeasible;
      else if (!domain.getChangedCols().empty())
        status = lp.resolveLp(&domain);

      if (status == HighsLpRelaxation::Status::Infeasible) {
        pruned_treeweight = 1.0;
        lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
        total_lp_iterations = lp.getNumLpIterations();
        sepa_lp_iterations = total_lp_iterations - firstrootlpiters;
        num_nodes = 1;
        num_leaves = 1;
        return true;
      }
    }
  }

  if (mipsolver.mipdata_->lower_bound > mipsolver.mipdata_->upper_limit) {
    lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
    total_lp_iterations = lp.getNumLpIterations();
    pruned_treeweight = 1.0;
    num_nodes = 1;
    num_leaves = 1;
    return true;
  }

  return false;
}

void HighsMipSolverData::evaluateRootNode() {
  int maxSepaRounds = mipsolver.submip ? 5 : HIGHS_CONST_I_INF;
restart:
  // solve the first root lp
  highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::INFO,
               "\nSolving root node LP relaxation\n");
  lp.loadModel();

  // add all cuts again after restart
  if (cutpool.getNumCuts() != 0) {
    highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                 "Adding %d cuts to LP after restart\n", cutpool.getNumCuts());
    assert(numRestarts != 0);
    HighsCutSet cutset;
    cutpool.separateLpCutsAfterRestart(cutset);
    lp.addCuts(cutset);
  }

  if (firstrootbasis.valid_) lp.getLpSolver().setBasis(firstrootbasis);
  lp.getLpSolver().setHighsOptionValue("presolve", "on");

  lp.getLpSolver().setHighsOptionValue("output_flag",
                                       mipsolver.options_mip_->output_flag);
  //  lp.getLpSolver().setHighsOptionValue("log_dev_level", LOG_DEV_LEVEL_INFO);
  //  lp.getLpSolver().setHighsOptionValue("log_file",
  //  mipsolver.options_mip_->log_file);
  HighsLpRelaxation::Status status = lp.resolveLp();

  lp.getLpSolver().setHighsOptionValue("output_flag", false);

  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  maxrootlpiters = std::max(lp.getNumLpIterations(), maxrootlpiters);
  if (numRestarts == 0) firstrootlpiters = lp.getNumLpIterations();

  lp.setIterationLimit(std::max(10000, int(50 * maxrootlpiters)));
  //  lp.getLpSolver().setHighsOptionValue("output_flag", false);
  //  lp.getLpSolver().setHighsOptionValue("log_dev_level", 0);
  lp.getLpSolver().setHighsOptionValue("parallel", "off");

  firstlpsol = lp.getLpSolver().getSolution().col_value;
  firstlpsolobj = lp.getObjective();
  firstrootbasis = lp.getLpSolver().getBasis();
  rootlpsolobj = firstlpsolobj;

  if (lp.unscaledDualFeasible(lp.getStatus())) {
    lower_bound = lp.getObjective();
    redcostfixing.addRootRedcost(
        mipsolver, lp.getLpSolver().getSolution().col_dual, lower_bound);
    if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF)
      redcostfixing.propagateRootRedcost(mipsolver);
  }

  if (!domain.infeasible()) {
    heuristics.randomizedRounding(firstlpsol);
    heuristics.flushStatistics();
  }

  if (status == HighsLpRelaxation::Status::Infeasible ||
      mipsolver.mipdata_->domain.infeasible() ||
      mipsolver.mipdata_->lower_bound > mipsolver.mipdata_->upper_limit) {
    lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
    total_lp_iterations = lp.getNumLpIterations();
    sepa_lp_iterations =
        total_lp_iterations - firstrootlpiters - heuristic_lp_iterations;
    pruned_treeweight = 1.0;
    num_nodes = 1;
    num_leaves = 1;
    return;
  }

  // begin separation
  std::vector<double> avgdirection;
  std::vector<double> curdirection;
  avgdirection.resize(mipsolver.numCol());
  curdirection.resize(mipsolver.numCol());

  int stall = 0;
  double smoothprogress = 0.0;
  int nseparounds = 0;

  HighsSeparation sepa(mipsolver);
  sepa.setLpRelaxation(&lp);

  total_lp_iterations = lp.getNumLpIterations();
  sepa_lp_iterations =
      total_lp_iterations - firstrootlpiters - heuristic_lp_iterations;

  while (lp.scaledOptimal(status) && !lp.getFractionalIntegers().empty() &&
         stall < 3) {
    printDisplayLine();
    if (checkLimits()) return;

    removeFixedIndices();

    if (mipsolver.options_mip_->presolve != off_string) {
      double fixingRate =
          100.0 * (1.0 - double(integer_cols.size() +
                                cliquetable.getSubstitutions().size()) /
                             numintegercols);
      if (fixingRate >= 10.0) {
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                     "%.1f%% inactive integer columns, restarting\n",
                     fixingRate);
        performRestart();
        if (mipsolver.modelstatus_ == HighsModelStatus::NOTSET) goto restart;

        return;
      }
    }

    ++nseparounds;

    int ncuts;
    if (rootSeparationRound(sepa, ncuts, status)) return;

    HighsCDouble sqrnorm = 0.0;
    const auto& solvals = lp.getSolution().col_value;

    for (int i = 0; i != mipsolver.numCol(); ++i) {
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
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      avgdirection[i] += scale * curdirection[i];
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
    if (lp.unscaledDualFeasible(status)) lower_bound = lp.getObjective();

    total_lp_iterations = lp.getNumLpIterations();
    sepa_lp_iterations =
        total_lp_iterations - firstrootlpiters - heuristic_lp_iterations;

    lp.setIterationLimit(std::max(10000, int(50 * maxrootlpiters)));
    if (ncuts == 0) break;
    if (nseparounds == maxSepaRounds) break;
  }

  lp.setIterationLimit();
  status = lp.resolveLp(&domain);
  total_lp_iterations = lp.getNumLpIterations();
  sepa_lp_iterations =
      total_lp_iterations - firstrootlpiters - heuristic_lp_iterations;
  if (status == HighsLpRelaxation::Status::Optimal &&
      lp.getFractionalIntegers().empty()) {
    mipsolver.modelstatus_ = HighsModelStatus::OPTIMAL;
    pruned_treeweight = 1.0;
    num_nodes = 1;
    num_leaves = 1;
    addIncumbent(lp.getLpSolver().getSolution().col_value, lp.getObjective(),
                 'T');
    return;
  } else {
    rootlpsol = lp.getLpSolver().getSolution().col_value;
    rootlpsolobj = lp.getObjective();
    lp.setIterationLimit(std::max(10000, int(50 * maxrootlpiters)));

    if ((!mipsolver.submip && numRestarts == 0) ||
        upper_limit == HIGHS_CONST_INF) {
      heuristics.centralRounding();
      heuristics.flushStatistics();
    }

    if (moreHeuristicsAllowed() || upper_limit == HIGHS_CONST_INF) {
      heuristics.RENS(rootlpsol);
      heuristics.flushStatistics();

      if (upper_limit == HIGHS_CONST_INF && !mipsolver.submip) {
        heuristics.feasibilityPump();
        heuristics.flushStatistics();
      }
    }
  }

  // if global propagation found bound changes, we update the local domain
  if (!domain.getChangedCols().empty()) {
    int ncuts;
    if (rootSeparationRound(sepa, ncuts, status)) {
      return;
    }

    if (lp.unscaledDualFeasible(status)) lower_bound = lp.getObjective();

    printDisplayLine();
  }

  removeFixedIndices();
  lp.removeObsoleteRows();
  rootlpsolobj = lp.getObjective();

  if (lower_bound <= upper_limit) {
    if (mipsolver.options_mip_->presolve != off_string) {
      double fixingRate =
          100.0 * (1.0 - double(integer_cols.size() +
                                cliquetable.getSubstitutions().size()) /
                             numintegercols);
      if (fixingRate >= 2.5 ||
          (!mipsolver.submip && fixingRate > 0 && numRestarts == 0)) {
        highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                     "%.1f%% inactive integer columns, restarting\n",
                     fixingRate);
        maxSepaRounds = std::min(maxSepaRounds, nseparounds);
        performRestart();
        if (mipsolver.modelstatus_ == HighsModelStatus::NOTSET) goto restart;

        return;
      }
    }
    // add the root node to the nodequeue to initialize the search
    nodequeue.emplaceNode(std::vector<HighsDomainChange>(), lower_bound,
                          lp.getObjective(), lp.getObjective(), 1);
  }
}

bool HighsMipSolverData::checkLimits() const {
  const HighsOptions& options = *mipsolver.options_mip_;
  if (options.mip_max_nodes != HIGHS_CONST_I_INF &&
      num_nodes >= size_t(options.mip_max_nodes)) {
    highsLogDev(options.log_options, HighsLogType::INFO,
                "reached node limit\n");
    mipsolver.modelstatus_ = HighsModelStatus::REACHED_ITERATION_LIMIT;
    return true;
  }
  if (options.mip_max_leaves != HIGHS_CONST_I_INF &&
      num_leaves >= size_t(options.mip_max_leaves)) {
    highsLogDev(options.log_options, HighsLogType::INFO,
                "reached leave node limit\n");
    mipsolver.modelstatus_ = HighsModelStatus::REACHED_ITERATION_LIMIT;
    return true;
  }
  if (mipsolver.timer_.read(mipsolver.timer_.solve_clock) >=
      options.time_limit) {
    highsLogDev(options.log_options, HighsLogType::INFO,
                "reached time limit\n");
    mipsolver.modelstatus_ = HighsModelStatus::REACHED_TIME_LIMIT;
    return true;
  }

  return false;
}

void HighsMipSolverData::checkObjIntegrality() {
  objintscale = 600.0;

  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.colCost(i) == 0.0) continue;

    if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS) {
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
    for (int i = 0; i != mipsolver.numCol(); ++i) {
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

    highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                "objective is always integral with scale %g\n", objintscale);
  }
}

void HighsMipSolverData::setupDomainPropagation() {
  const HighsLp& model = *mipsolver.model_;
  highsSparseTranspose(model.numRow_, model.numCol_, model.Astart_,
                       model.Aindex_, model.Avalue_, ARstart_, ARindex_,
                       ARvalue_);

  // compute the maximal absolute coefficients to filter propagation
  maxAbsRowCoef.resize(mipsolver.model_->numRow_);
  for (int i = 0; i != mipsolver.model_->numRow_; ++i) {
    double maxabsval = 0.0;

    int start = ARstart_[i];
    int end = ARstart_[i + 1];
    for (int j = start; j != end; ++j)
      maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));

    maxAbsRowCoef[i] = maxabsval;
  }

  domain = HighsDomain(mipsolver);
  domain.computeRowActivities();
}
