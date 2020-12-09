#include "mip/HighsMipSolverData.h"

static int64_t gcd(int64_t a, int64_t b) {
  int h;
  if (a < 0) a = -a;
  if (b < 0) b = -b;

  if (a == 0) return b;
  if (b == 0) return a;

  do {
    h = a % b;
    a = b;
    b = h;
  } while (b != 0);

  return a;
}
void HighsMipSolverData::setup() {
  const HighsLp& model = *mipsolver.model_;

#ifdef HIGHS_DEBUGSOL
  if (!mipsolver.options_mip_->mip_debug_solution_file.empty()) {
    HighsPrintMessage(mipsolver.options_mip_->output,
                      mipsolver.options_mip_->message_level, ML_MINIMAL,
                      "reading debug solution file %s\n",
                      mipsolver.options_mip_->mip_debug_solution_file.c_str());
    std::ifstream file(mipsolver.options_mip_->mip_debug_solution_file);
    if (file) {
      std::string varname;
      double varval;
      std::map<std::string, int> nametoidx;

      for (int i = 0; i != mipsolver.model_->numCol_; ++i)
        nametoidx[mipsolver.model_->col_names_[i]] = i;

      highsDebugSolution.resize(mipsolver.model_->numCol_, 0.0);
      while (!file.eof()) {
        file >> varname;
        auto it = nametoidx.find(varname);
        if (it != nametoidx.end()) {
          file >> varval;
          HighsPrintMessage(mipsolver.options_mip_->output,
                            mipsolver.options_mip_->message_level, ML_MINIMAL,
                            "%s = %g\n", varname.c_str(), varval);
          highsDebugSolution[it->second] = varval;
        }

        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    } else {
      throw std::runtime_error("debug solution: could not open file\n");
    }
  }
#endif

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
  num_disp_lines = 0;
  last_displeave = 0;
  lower_bound = -HIGHS_CONST_INF;
  upper_bound = HIGHS_CONST_INF;
  upper_limit = mipsolver.options_mip_->dual_objective_value_upper_bound -
                mipsolver.model_->offset_;

  if (mipsolver.submip) pseudocost.setMinReliable(0);

  if (mipsolver.options_mip_->mip_report_level == 0)
    dispfreq = 0;
  else if (mipsolver.options_mip_->mip_report_level == 1)
    dispfreq = 10;
  else
    dispfreq = 1;

  // create row-wise copy of model
  highsSparseTranspose(model.numRow_, model.numCol_, model.Astart_,
                       model.Aindex_, model.Avalue_, ARstart_, ARindex_,
                       ARvalue_);

  // if a root basis is given, construct a basis for the root LP from
  // in the reduced problem space after presolving
  if (mipsolver.rootbasis) {
    if (mipsolver.presolve_.data_.presolve_.empty())
      firstrootbasis = *mipsolver.rootbasis;
    else {
      int numColOrig = mipsolver.presolve_.data_.presolve_[0].numColOriginal;
      int numRowOrig = mipsolver.presolve_.data_.presolve_[0].numRowOriginal;

      const auto& rIndex = mipsolver.presolve_.data_.presolve_[0].rIndex;
      const auto& cIndex = mipsolver.presolve_.data_.presolve_[0].cIndex;
      firstrootbasis.valid_ = true;
      firstrootbasis.col_status.resize(mipsolver.numCol(),
                                       HighsBasisStatus::NONBASIC);
      firstrootbasis.row_status.resize(mipsolver.numRow(),
                                       HighsBasisStatus::NONBASIC);

      int missingbasic = mipsolver.model_->numRow_;

      for (int i = 0; i != numRowOrig; ++i) {
        if (rIndex[i] < 0) continue;

        HighsBasisStatus rowstatus = mipsolver.rootbasis->row_status[i];

        if (rowstatus == HighsBasisStatus::BASIC) {
          if (missingbasic != 0)
            --missingbasic;
          else
            rowstatus = HighsBasisStatus::NONBASIC;
        }

        firstrootbasis.row_status[rIndex[i]] = rowstatus;
      }

      for (int i = 0; i != numColOrig; ++i) {
        if (cIndex[i] < 0) continue;

        HighsBasisStatus colstatus = mipsolver.rootbasis->col_status[i];

        if (colstatus == HighsBasisStatus::BASIC) {
          if (missingbasic != 0)
            --missingbasic;
          else
            colstatus = HighsBasisStatus::NONBASIC;
        }

        firstrootbasis.col_status[cIndex[i]] = colstatus;
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
      currgcd = gcd(intval, currgcd);
      if (currgcd == 1) break;
    }

    if (currgcd != 0) objintscale /= currgcd;

    HighsPrintMessage(mipsolver.options_mip_->output,
                      mipsolver.options_mip_->message_level, ML_MINIMAL,
                      "objective is always integral with scale %g\n",
                      objintscale);
  }

  // compute row activities and propagate all rows once
  domain.computeRowActivities();
  domain.propagate();

  // extract cliques
  cliquetable.extractCliques(mipsolver);
  if (domain.infeasible()) {
    lower_bound = HIGHS_CONST_INF;
    pruned_treeweight = 1.0;
    return;
  }
  // store binary variables in vector with their number of implications on other
  // binaries
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

      if (domain.isBinary(i)) {
        --contingent;
        if (contingent < 0) break;
        ++nprobed;
        bool fixed = implications.runProbing(i, contingent);

        if (domain.infeasible()) {
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

    HighsPrintMessage(mipsolver.options_mip_->output,
                      mipsolver.options_mip_->message_level, ML_MINIMAL,
                      "%d probing evaluations: %lu bound changes and %d fixed "
                      "binary variables\n",
                      nprobed, domain.getChangedCols().size(), nfixed);

    cliquetable.cleanupFixed(domain);
  }
  if (!domain.getChangedCols().empty()) lp.flushDomain(domain, true);

  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  lp.getLpSolver().setHighsOptionValue("dual_simplex_cleanup_strategy", 2);
  // lp.getLpSolver().setHighsOptionValue("parallel", "on");
  lp.getLpSolver().setHighsOptionValue("simplex_initial_condition_check",
                                       false);
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
      upper_limit = new_upper_limit;
      pruned_treeweight += nodequeue.performBounding(upper_limit);
      printDisplayLine(source);
    }
  }
}

void HighsMipSolverData::printDisplayLine(char first) {
  double offset = mipsolver.model_->offset_;
  if (num_disp_lines % 20 == 0) {
    HighsPrintMessage(
        mipsolver.options_mip_->output, mipsolver.options_mip_->message_level,
        ML_MINIMAL,
        "   %7s | %10s | %10s | %10s | %10s | %-14s | %-14s | %7s | %7s "
        "| %8s | %8s\n",
        "time", "open nodes", "nodes", "leaves", "lpiters", "dual bound",
        "primal bound", "cutpool", "lpcuts", "gap", "progress");
  }

  ++num_disp_lines;
  last_displeave = num_leaves;

  double lb = mipsolver.mipdata_->lower_bound + offset;
  double ub = HIGHS_CONST_INF;
  double gap = HIGHS_CONST_INF;
  int lpcuts =
      mipsolver.mipdata_->lp.getNumLpRows() - mipsolver.model_->numRow_;

  if (upper_bound != HIGHS_CONST_INF) {
    ub = upper_bound + offset;
    lb = std::min(ub, lb);
    gap = 100 * (ub - lb) / std::max(1.0, std::abs(ub));

    HighsPrintMessage(
        mipsolver.options_mip_->output, mipsolver.options_mip_->message_level,
        ML_MINIMAL,
        " %c %6.1fs | %10lu | %10lu | %10lu | %10lu | %-14.9g | %-14.9g | "
        "%7d | %7d | %7.2f%% | %7.2f%%\n",
        first, mipsolver.timer_.read(mipsolver.timer_.solve_clock),
        nodequeue.numNodes(), num_nodes, num_leaves, total_lp_iterations, lb,
        ub, mipsolver.mipdata_->cutpool.getNumCuts(), lpcuts, gap,
        100 * double(pruned_treeweight));
  } else {
    HighsPrintMessage(
        mipsolver.options_mip_->output, mipsolver.options_mip_->message_level,
        ML_MINIMAL,
        " %c %6.1fs | %10lu | %10lu | %10lu | %10lu | %-14.9g | %-14.9g | "
        "%7d | %7d | %8.2f | %7.2f%%\n",
        first, mipsolver.timer_.read(mipsolver.timer_.solve_clock),
        nodequeue.numNodes(), num_nodes, num_leaves, total_lp_iterations, lb,
        ub, mipsolver.mipdata_->cutpool.getNumCuts(), lpcuts, gap,
        100 * double(pruned_treeweight));
  }
}

void HighsMipSolverData::evaluateRootNode() {
  // solve the first root lp

  HighsPrintMessage(mipsolver.options_mip_->output,
                    mipsolver.options_mip_->message_level, ML_MINIMAL,
                    "\nsolving root node LP relaxation\n");
  lp.getLpSolver().setHighsOptionValue("presolve", "on");
  lp.getLpSolver().setHighsOptionValue("message_level",
                                       mipsolver.options_mip_->message_level);
  lp.getLpSolver().setHighsLogfile(mipsolver.options_mip_->logfile);
  lp.resolveLp();
  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  maxrootlpiters = lp.getNumLpIterations();

  lp.setIterationLimit(std::max(1000, int(50 * maxrootlpiters)));
  lp.getLpSolver().setHighsLogfile(NULL);
  lp.getLpSolver().setHighsOptionValue("message_level", 0);
  lp.getLpSolver().setHighsOptionValue("parallel", "off");

  firstlpsol = lp.getLpSolver().getSolution().col_value;
  firstlpsolobj = lp.getObjective();
  firstrootbasis = lp.getLpSolver().getBasis();
  rootlpsolobj = firstlpsolobj;

  if (lp.unscaledDualFeasible(lp.getStatus()))
    mipsolver.mipdata_->lower_bound = lp.getObjective();

  // begin separation
  std::vector<double> avgdirection;
  std::vector<double> curdirection;
  avgdirection.resize(mipsolver.numCol());
  curdirection.resize(mipsolver.numCol());

  int stall = 0;
  double smoothprogress = 0.0;
  int nseparounds = 0;

  HighsLpRelaxation::Status status = lp.getStatus();

  HighsSeparation sepa;
  sepa.setLpRelaxation(&lp);

  total_lp_iterations = lp.getNumLpIterations();

  while (lp.scaledOptimal(status) && !lp.getFractionalIntegers().empty() &&
         stall < 3) {
    ++nseparounds;
    printDisplayLine();
    size_t tmpilpiters = lp.getNumLpIterations();
    maxrootlpiters =
        std::max(maxrootlpiters, lp.getNumLpIterations() - tmpilpiters);

    int ncuts = sepa.separationRound(domain, status);
    if (status == HighsLpRelaxation::Status::Infeasible) {
      pruned_treeweight = 1.0;
      lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
      total_lp_iterations = lp.getNumLpIterations();
      num_nodes = 1;
      num_leaves = 1;
      return;
    }

    const std::vector<double>& solvals =
        lp.getLpSolver().getSolution().col_value;

    HighsCDouble sqrnorm = 0.0;
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

    if (lp.unscaledDualFeasible(status))
      mipsolver.mipdata_->lower_bound = lp.getObjective();
    total_lp_iterations = lp.getNumLpIterations();

    lp.setIterationLimit(int(50 * maxrootlpiters));

    if (mipsolver.mipdata_->lower_bound > mipsolver.mipdata_->upper_limit) {
      lower_bound = std::min(HIGHS_CONST_INF, upper_bound);
      total_lp_iterations = lp.getNumLpIterations();
      pruned_treeweight = 1.0;
      num_nodes = 1;
      num_leaves = 1;
      return;
    }

    if (ncuts == 0) break;
    if (mipsolver.submip && nseparounds == 5) break;
  }

  // remove inactive cuts and resolve final root LP without iteration limit
  lp.setIterationLimit();
  cutpool.removeObsoleteRows(lp);

  status = lp.resolveLp();
  total_lp_iterations = lp.getNumLpIterations();
  if (status == HighsLpRelaxation::Status::Optimal &&
      lp.getFractionalIntegers().empty()) {
    pruned_treeweight = 1.0;
    num_nodes = 1;
    num_leaves = 1;
    addIncumbent(lp.getLpSolver().getSolution().col_value, lp.getObjective(),
                 'T');
  } else {
    rootlpsol = lp.getLpSolver().getSolution().col_value;
    rootlpsolobj = lp.getObjective();
    lp.setIterationLimit(int(50 * maxrootlpiters));

    // add the root node to the nodequeue to initialize the search
    nodequeue.emplaceNode(std::vector<HighsDomainChange>(), lp.getObjective(),
                          lp.getObjective(), 1);
  }
}