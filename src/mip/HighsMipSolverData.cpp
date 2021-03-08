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

  addIncumbent(solution, double(obj), source);
  return true;
}

HighsMipSolverData::ModelCleanup::ModelCleanup(HighsMipSolver& mipsolver) {
  origmodel = mipsolver.model_;
  HighsLp model;

  // we initialize everything except the matrix and column bounds,
  // as they are modified first via HAggregator and HighsDomain
  // and the copy would be thrown away
  model.colCost_ = origmodel->colCost_;
  model.integrality_ = origmodel->integrality_;
  model.numCol_ = origmodel->numCol_;
  model.numRow_ = origmodel->numRow_;
  model.rowLower_ = origmodel->rowLower_;
  model.rowUpper_ = origmodel->rowUpper_;
  model.sense_ = origmodel->sense_;
  model.lp_name_ = origmodel->lp_name_;
  model.model_name_ = origmodel->model_name_;
  model.offset_ = origmodel->offset_;
  model.orientation_ = origmodel->orientation_;

  std::vector<double>& colLower = mipsolver.mipdata_->domain.colLower_;
  std::vector<double>& colUpper = mipsolver.mipdata_->domain.colUpper_;

  std::vector<int> coldeleted(model.numCol_);
  std::vector<uint8_t> rowdeleted(model.numRow_);

  presolve::HAggregator aggregator(model.rowLower_, model.rowUpper_,
                                   model.colCost_, model.offset_,
                                   model.integrality_, colLower, colUpper);
  aggregator.fromCSC(origmodel->Avalue_, origmodel->Aindex_,
                     origmodel->Astart_);

  origsol.resize(model.numCol_);
  // printf("before model cleanup: %d nonzeros\n", aggregator.numNonzeros());
  int nfixed = 0;
  for (int i = 0; i != model.numCol_; ++i) {
    if (colLower[i] != colUpper[i]) continue;
    aggregator.removeFixedCol(i);
    ++nfixed;
    origsol[i] = colLower[i];
    coldeleted[i] = -1;
  }

  // printf("removed %d fixed columns: %d nonzeros\n", nfixed,
  //       aggregator.numNonzeros());

  for (int delrow : mipsolver.mipdata_->cliquetable.getDeletedRows()) {
    if (!rowdeleted[delrow]) {
      aggregator.removeRow(delrow);
      rowdeleted[delrow] = true;
    }
  }
  mipsolver.mipdata_->cliquetable.getDeletedRows().clear();

  aggregator.removeRedundantRows(rowdeleted);

  int nremoved = 0;
  for (int i = 0; i != model.numRow_; ++i) {
    if (rowdeleted[i]) ++nremoved;
  }

  // printf("removed %d redundant rows: %d nonzeros\n", nremoved,
  //       aggregator.numNonzeros());
  std::vector<std::pair<int, HighsCliqueTable::CliqueVar>>& extensionvars =
      mipsolver.mipdata_->cliquetable.getCliqueExtensions();
  int addednnz = extensionvars.size();
  for (std::pair<int, HighsCliqueTable::CliqueVar> cliqueextension :
       extensionvars) {
    if (rowdeleted[cliqueextension.first]) {
      --addednnz;
      continue;
    }
    double val;
    if (cliqueextension.second.val == 0) {
      model.rowLower_[cliqueextension.first] -= 1;
      model.rowUpper_[cliqueextension.first] -= 1;
      val = -1.0;
    } else
      val = 1.0;
    aggregator.addNonzero(cliqueextension.first, cliqueextension.second.col,
                          val);
  }
  extensionvars.clear();
  // printf("clique extension added %d nonzeros\n", addednnz);

  int ndelrows = 0;

  for (int delrow : mipsolver.mipdata_->cliquetable.getDeletedRows()) {
    if (!rowdeleted[delrow]) {
      ++ndelrows;
      aggregator.removeRow(delrow);
      rowdeleted[delrow] = true;
    }
  }

  nremoved += ndelrows;
  mipsolver.mipdata_->cliquetable.getDeletedRows().clear();

  // run clique extension before performing substitutions

  // printf("after clique extension: %d fixed cols, %d nonzeros\n", nfixed,
  //       aggregator.numNonzeros());

  // clean up fixings
  mipsolver.mipdata_->cliquetable.cleanupFixed(mipsolver.mipdata_->domain);
  for (int i = 0; i != model.numCol_; ++i) {
    if (coldeleted[i] != 0 ||
        std::abs(colLower[i] - colUpper[i]) > HIGHS_CONST_TINY)
      continue;
    aggregator.removeFixedCol(i);

    ++nfixed;
    origsol[i] = colLower[i];
    coldeleted[i] = -1;
  }

  // printf("removed fixed columns: %d fixed cols, %d nonzeros\n", nfixed,
  //       aggregator.numNonzeros());
  int nsubst = 0;

  for (const HighsSubstitution& substitution :
       mipsolver.mipdata_->implications.substitutions) {
    if (coldeleted[substitution.substcol] == -1 ||
        coldeleted[substitution.staycol] == -1)
      continue;
    aggregator.substitute(substitution.substcol, substitution.staycol,
                          substitution.offset, substitution.scale);

    substitutionStack.push_back(substitution);
    coldeleted[substitution.substcol] = substitutionStack.size();
    ++nsubst;
  }

  // printf("after %d non-binary probing substitutions: %d nonzeros\n", nsubst,
  //       aggregator.numNonzeros());
  nsubst = 0;

  for (HighsCliqueTable::Substitution subst :
       mipsolver.mipdata_->cliquetable.getSubstitutions()) {
    if (coldeleted[subst.substcol] == -1) continue;
    assert(coldeleted[subst.substcol] == 0);

    double scale;
    double offset;

    if (subst.replace.val == 0) {
      scale = -1.0;
      offset = 1.0;
    } else {
      scale = 1.0;
      offset = 0.0;
    }

    HighsSubstitution substitution{subst.substcol, int(subst.replace.col),
                                   scale, double(offset)};

    aggregator.substitute(substitution.substcol, substitution.staycol,
                          substitution.offset, substitution.scale);

    substitutionStack.push_back(substitution);
    coldeleted[substitution.substcol] = substitutionStack.size();
    ++nsubst;
  }

  // printf("after %d cliquetable substitutions: %d nonzeros\n", nsubst,
  //       aggregator.numNonzeros());

  aggregator.removeRedundantRows(rowdeleted);
  nremoved = 0;
  for (int i = 0; i != model.numRow_; ++i) {
    if (rowdeleted[i]) ++nremoved;
  }

  int numstrengthened = aggregator.strengthenInequalities();

  if (numstrengthened != 0)
    highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                "strengthened %d coefficients\n", numstrengthened);

  // printf("removed redundant rows: %d removed rows, %d nonzeros\n", nremoved,
  //       aggregator.numNonzeros());

  aggregator.toCSC(model.Avalue_, model.Aindex_, model.Astart_);

  model.colLower_ = std::move(colLower);
  model.colUpper_ = std::move(colUpper);
  cleanedUpModel = std::move(model);
  mipsolver.model_ = &cleanedUpModel;

  cIndex.assign(cleanedUpModel.numCol_, -1);
  rIndex.assign(cleanedUpModel.numRow_, -1);

  int numreducedcol = 0;
  for (int i = 0; i != cleanedUpModel.numCol_; ++i) {
    if (coldeleted[i]) continue;
    cIndex[i] = numreducedcol++;
  }

  int numreducedrow = 0;
  for (int i = 0; i != cleanedUpModel.numRow_; ++i) {
    if (rowdeleted[i]) continue;
    rIndex[i] = numreducedrow++;
  }

  for (int i = 0; i != cleanedUpModel.numCol_; ++i) {
    int newindex = cIndex[i];
    if (newindex == -1) continue;
    cleanedUpModel.Astart_[newindex] = cleanedUpModel.Astart_[i];
    cleanedUpModel.colCost_[newindex] = cleanedUpModel.colCost_[i];
    cleanedUpModel.integrality_[newindex] = cleanedUpModel.integrality_[i];
    cleanedUpModel.colLower_[newindex] = cleanedUpModel.colLower_[i];
    cleanedUpModel.colUpper_[newindex] = cleanedUpModel.colUpper_[i];
  }
  cleanedUpModel.Astart_[numreducedcol] = cleanedUpModel.Avalue_.size();
  cleanedUpModel.Astart_.resize(numreducedcol + 1);
  cleanedUpModel.colCost_.resize(numreducedcol);
  cleanedUpModel.integrality_.resize(numreducedcol);
  cleanedUpModel.colLower_.resize(numreducedcol);
  cleanedUpModel.colUpper_.resize(numreducedcol);

  for (int i = 0; i != cleanedUpModel.numRow_; ++i) {
    int newindex = rIndex[i];
    if (newindex == -1) continue;
    cleanedUpModel.rowLower_[newindex] = cleanedUpModel.rowLower_[i];
    cleanedUpModel.rowUpper_[newindex] = cleanedUpModel.rowUpper_[i];
  }
  cleanedUpModel.rowLower_.resize(numreducedrow);
  cleanedUpModel.rowUpper_.resize(numreducedrow);

  int numnnz = cleanedUpModel.Avalue_.size();
  for (int i = 0; i != numnnz; ++i) {
    cleanedUpModel.Aindex_[i] = rIndex[cleanedUpModel.Aindex_[i]];
    assert(cleanedUpModel.Aindex_[i] < numreducedrow);
    assert(cleanedUpModel.Aindex_[i] >= 0);
  }

  cleanedUpModel.numCol_ = numreducedcol;
  cleanedUpModel.numRow_ = numreducedrow;

  mipsolver.mipdata_->rowMatrixSet = false;
  mipsolver.mipdata_->cutpool = HighsCutPool(
      cleanedUpModel.numCol_, mipsolver.options_mip_->mip_pool_age_limit,
      mipsolver.options_mip_->mip_pool_soft_limit);
  mipsolver.mipdata_->domain = HighsDomain(mipsolver);
  mipsolver.mipdata_->domain.addCutpool(mipsolver.mipdata_->cutpool);
  mipsolver.mipdata_->pseudocost = HighsPseudocost(cleanedUpModel.numCol_);
  mipsolver.mipdata_->cliquetable.rebuild(cleanedUpModel.numCol_, cIndex,
                                          rIndex);
  mipsolver.mipdata_->implications.rebuild(cleanedUpModel.numCol_, cIndex,
                                           rIndex);

  reportPresolveReductions(mipsolver.options_mip_->log_options, *origmodel,
                           cleanedUpModel);
}

void HighsMipSolverData::ModelCleanup::recoverSolution(
    const std::vector<double>& reducedSol) {
  if (int(reducedSol.size()) != cleanedUpModel.numCol_) return;
  for (int i = 0; i != origmodel->numCol_; ++i) {
    int col = cIndex[i];
    if (col < 0) continue;
    origsol[i] = reducedSol[col];
  }

  for (int k = substitutionStack.size() - 1; k >= 0; --k) {
    origsol[substitutionStack[k].substcol] =
        substitutionStack[k].offset +
        substitutionStack[k].scale * origsol[substitutionStack[k].staycol];
  }
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
  feastol = mipsolver.options_mip_->mip_feasibility_tolerance;
  epsilon = mipsolver.options_mip_->mip_epsilon;
  heuristic_effort = mipsolver.options_mip_->mip_heuristic_effort;

  firstlpsolobj = -HIGHS_CONST_INF;
  rootlpsolobj = -HIGHS_CONST_INF;

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
  tryProbing = true;
  lower_bound = -HIGHS_CONST_INF;
  upper_bound = HIGHS_CONST_INF;

  if (mipsolver.submip) pseudocost.setMinReliable(0);

  if (mipsolver.options_mip_->mip_report_level == 0)
    dispfreq = 0;
  else if (mipsolver.options_mip_->mip_report_level == 1)
    dispfreq = 50;
  else
    dispfreq = 1;
}

void HighsMipSolverData::runSetup() {
  const HighsLp& model = *mipsolver.model_;

  upper_limit = mipsolver.options_mip_->dual_objective_value_upper_bound -
                mipsolver.model_->offset_;

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

  if (mipsolver.submip) pseudocost.setMinReliable(0);

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

  if (model.numCol_ == 0) {
    mipsolver.modelstatus_ = HighsModelStatus::OPTIMAL;
    return;
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

  if (checkLimits()) return;
  // extract cliques if they have not been extracted before
  if (!cliquesExtracted) {
    cliquesExtracted = true;
    cliquetable.extractCliques(mipsolver);
    if (!domain.infeasible() && upper_limit != HIGHS_CONST_INF)
      cliquetable.extractObjCliques(mipsolver);
    if (domain.infeasible()) {
      mipsolver.modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
      lower_bound = HIGHS_CONST_INF;
      pruned_treeweight = 1.0;
      return;
    }
  }

  if (tryProbing) {
    runProbing();
    if (domain.infeasible()) {
      mipsolver.modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
      lower_bound = HIGHS_CONST_INF;
      pruned_treeweight = 1.0;
      return;
    }
    if (cliquetable.getNumFixings() == 0 &&
        cliquetable.getSubstitutions().empty() &&
        implications.substitutions.empty())
      tryProbing = false;
  }

  if (!modelcleanup) {
    int nfixed = cliquetable.getNumFixings();
    if (nfixed == 0) {
      for (int i = 0; i != model.numCol_; ++i)
        if (domain.isFixed(i)) ++nfixed;
    }

    if (nfixed != 0 || !cliquetable.getDeletedRows().empty() ||
        !cliquetable.getSubstitutions().empty() ||
        !implications.substitutions.empty() ||
        !cliquetable.getCliqueExtensions().empty()) {
      modelcleanup = decltype(modelcleanup)(new ModelCleanup(mipsolver));
      runSetup();
      return;
    }
  }

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

  continuous_cols.clear();
  integer_cols.clear();
  implint_cols.clear();
  integral_cols.clear();
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

  debugSolution.activate();
}

void HighsMipSolverData::basisTransfer() {
  // if a root basis is given, construct a basis for the root LP from
  // in the reduced problem space after presolving
  if (mipsolver.rootbasis) {
    if (mipsolver.presolve_.data_.presolve_.empty() && !modelcleanup)
      firstrootbasis = *mipsolver.rootbasis;
    else {
      int numColOrig;
      int numRowOrig;

      const HighsLp& model = *mipsolver.model_;

      if (!mipsolver.presolve_.data_.presolve_.empty()) {
        numColOrig = mipsolver.presolve_.data_.presolve_[0].numColOriginal;
        numRowOrig = mipsolver.presolve_.data_.presolve_[0].numRowOriginal;
      } else {
        numColOrig = modelcleanup->origmodel->numCol_;
        numRowOrig = modelcleanup->origmodel->numRow_;
      }

      const auto rIndex = [&](int row) {
        if (!mipsolver.presolve_.data_.presolve_.empty())
          row = mipsolver.presolve_.data_.presolve_[0].rIndex[row];

        if (modelcleanup && row != -1) row = modelcleanup->rIndex[row];

        return row;
      };

      const auto cIndex = [&](int col) {
        if (!mipsolver.presolve_.data_.presolve_.empty())
          col = mipsolver.presolve_.data_.presolve_[0].cIndex[col];

        if (modelcleanup && col != -1) col = modelcleanup->cIndex[col];

        return col;
      };

      firstrootbasis.valid_ = true;
      firstrootbasis.col_status.resize(mipsolver.numCol(),
                                       HighsBasisStatus::NONBASIC);
      firstrootbasis.row_status.resize(mipsolver.numRow(),
                                       HighsBasisStatus::NONBASIC);

      int missingbasic = mipsolver.model_->numRow_;

      for (int i = 0; i != numRowOrig; ++i) {
        int r = rIndex(i);
        if (r < 0) continue;

        HighsBasisStatus rowstatus = mipsolver.rootbasis->row_status[i];

        if (rowstatus == HighsBasisStatus::BASIC) {
          if (missingbasic != 0)
            --missingbasic;
          else
            rowstatus = HighsBasisStatus::NONBASIC;
        }

        firstrootbasis.row_status[r] = rowstatus;
      }

      for (int i = 0; i != numColOrig; ++i) {
        int c = cIndex(i);
        if (c < 0) continue;

        HighsBasisStatus colstatus = mipsolver.rootbasis->col_status[i];

        if (colstatus == HighsBasisStatus::BASIC) {
          if (missingbasic != 0)
            --missingbasic;
          else
            colstatus = HighsBasisStatus::NONBASIC;
        }

        firstrootbasis.col_status[c] = colstatus;
      }

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
        nonbasiccols.resize(
            std::min(nonbasiccols.size(), size_t(missingbasic)));
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
            if (firstrootbasis.row_status[i] == HighsBasisStatus::BASIC)
              continue;

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

    lp.getLpSolver().setBasis(firstrootbasis);
  }
}

const std::vector<double>& HighsMipSolverData::getSolution() const {
  if (modelcleanup) return modelcleanup->origsol;

  return incumbent;
}

void HighsMipSolverData::addIncumbent(const std::vector<double>& sol,
                                      double solobj, char source) {
  if (solobj < upper_bound) {
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
      debugSolution.newIncumbentFound();
      upper_limit = new_upper_limit;
      redcostfixing.propagateRootRedcost(mipsolver);
      cliquetable.extractObjCliques(mipsolver);
      pruned_treeweight += nodequeue.performBounding(upper_limit);
      printDisplayLine(source);
    }
  }
}

void HighsMipSolverData::printDisplayLine(char first) {
  double offset = mipsolver.model_->offset_;
  if (num_disp_lines % 20 == 0) {
    highsLogDev(
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

    highsLogDev(
        mipsolver.options_mip_->log_options, HighsLogType::INFO,
        " %c %6.1fs | %10lu | %10lu | %10lu | %10lu | %-14.9g | %-14.9g | "
        "%7d | %7d | %7.2f%% | %7.2f%%\n",
        first, mipsolver.timer_.read(mipsolver.timer_.solve_clock),
        nodequeue.numNodes(), num_nodes, num_leaves, total_lp_iterations, lb,
        ub, mipsolver.mipdata_->cutpool.getNumCuts(), lpcuts, gap,
        100 * double(pruned_treeweight));
  } else {
    highsLogDev(
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
  // solve the first root lp
  highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::INFO,
              "\nsolving root node LP relaxation\n");
  lp.loadModel();
  lp.getLpSolver().setHighsOptionValue("presolve", "on");

  //  lp.getLpSolver().setHighsOptionValue("log_dev_level", LOG_DEV_LEVEL_INFO);
  //  lp.getLpSolver().setHighsOptionValue("log_file",
  //  mipsolver.options_mip_->log_file);
  HighsLpRelaxation::Status status = lp.resolveLp();

  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  maxrootlpiters = lp.getNumLpIterations();
  firstrootlpiters = maxrootlpiters;

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
    if (mipsolver.submip && nseparounds == 5) break;
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

    heuristics.RENS(rootlpsol);
    heuristics.flushStatistics();

    if (upper_limit == HIGHS_CONST_INF && !mipsolver.submip) {
      heuristics.centralRounding();
      heuristics.flushStatistics();
      if (upper_limit == HIGHS_CONST_INF) {
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

  if (lower_bound <= upper_limit) {
    // add the root node to the nodequeue to initialize the search
    nodequeue.emplaceNode(std::vector<HighsDomainChange>(), lower_bound,
                          lp.getObjective(), lp.getObjective(), 1);
  }

  removeFixedIndices();
  lp.removeObsoleteRows();
  rootlpsolobj = lp.getObjective();
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

void HighsMipSolverData::runProbing() {
  const HighsLp& model = *mipsolver.model_;
  // store binary variables in vector with their number of implications on
  // other binaries
  std::vector<std::tuple<int, int, int>> binaries;
  binaries.reserve(model.numCol_);
  HighsRandom random;
  for (int i = 0; i != model.numCol_; ++i) {
    if (domain.isBinary(i))
      binaries.emplace_back(-std::min(100, cliquetable.getNumImplications(i)),
                            random.integer(), i);
  }
  if (!binaries.empty()) {
    // sort variables with many implications on other binaries first
    std::sort(binaries.begin(), binaries.end());

    int nfixed = 0;
    int contingent = 1000;
    int nprobed = 0;
    for (std::tuple<int, int, int> binvar : binaries) {
      int i = std::get<2>(binvar);

      if (cliquetable.getSubstitution(i) != nullptr) continue;

      if (domain.isBinary(i)) {
        if (nprobed % 16 == 1 && checkLimits()) return;
        --contingent;
        if (contingent < 0) break;
        ++nprobed;
        bool fixed = implications.runProbing(i, contingent);

        if (domain.infeasible()) {
          mipsolver.modelstatus_ = HighsModelStatus::PRIMAL_INFEASIBLE;
          lower_bound = HIGHS_CONST_INF;
          pruned_treeweight = 1.0;
          return;
        }
        if (fixed) {
          ++nfixed;
          contingent += nfixed;
        }
      } else
        ++nfixed;
    }

    highsLogDev(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                "%d probing evaluations: %d fixed binary variables, %d "
                "bound changes\n",
                nprobed, nfixed, int(domain.getChangedCols().size()) - nfixed);

    cliquetable.cleanupFixed(domain);
    if (!mipsolver.mipdata_->modelcleanup) cliquetable.runCliqueMerging(domain);
  }
}
