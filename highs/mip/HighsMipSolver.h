/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsCallback.h"
#include "lp_data/HighsOptions.h"
#include "mip/HighsMipAnalysis.h"

struct HighsMipSolverData;
class HighsCutPool;
struct HighsPseudocostInitialization;
class HighsCliqueTable;
class HighsImplications;

const HighsInt kMipRaceNoSolution = -1;

struct MipRaceIncumbent {
  HighsInt start_write_incumbent = kMipRaceNoSolution;
  HighsInt finish_write_incumbent = kMipRaceNoSolution;
  double objective = -kHighsInf;
  std::vector<double> solution;
  void clear();
  void initialise(const HighsInt num_col);
  void update(const double objective, const std::vector<double>& solution);
  HighsInt read(const HighsInt last_incumbent_read, double& objective_,
                std::vector<double>& solution_) const;
};

struct MipRaceRecord {
  std::vector<MipRaceIncumbent> incumbent;
  void clear();
  void initialise(const HighsInt mip_race_concurrency, const HighsInt num_col);
  HighsInt concurrency() const;
  void update(const HighsInt instance, const double objective,
              const std::vector<double>& solution);
  void report(const HighsLogOptions log_options) const;
};

struct MipRace {
  HighsInt my_instance;
  MipRaceRecord* record = nullptr;
  HighsLogOptions log_options;
  std::vector<HighsInt> last_incumbent_read;
  void clear();
  void initialise(const HighsInt mip_race_concurrency,
                  const HighsInt my_instance_, MipRaceRecord* record_,
                  const HighsLogOptions log_options_);
  HighsInt concurrency() const;
  void update(const double objective, const std::vector<double>& solution);
  bool newSolution(const HighsInt instance, double objective,
                   std::vector<double>& solution);
  void report() const;
};

struct HighsTerminator {
  HighsInt num_instance;
  HighsInt my_instance;
  HighsModelStatus* record;
  void clear();
  void initialise(HighsInt num_instance_, HighsInt my_instance_,
                  HighsModelStatus* record_);
  void terminate();
  HighsModelStatus terminated() const;
  bool notTerminated() const;
  void report(const HighsLogOptions log_options) const;
};

class HighsMipSolver {
 public:
  HighsCallback* callback_;
  const HighsOptions* options_mip_;
  const HighsLp* model_;
  const HighsLp* orig_model_;
  HighsModelStatus modelstatus_;
  std::vector<double> solution_;
  double solution_objective_;
  double bound_violation_;
  double integrality_violation_;
  double row_violation_;
  // The following are only to return data to HiGHS, and are set in
  // HighsMipSolver::cleanupSolve
  double dual_bound_;
  double primal_bound_;
  double gap_;
  int64_t node_count_;
  int64_t total_lp_iterations_;
  double primal_dual_integral_;

  FILE* improving_solution_file_;
  std::vector<HighsObjectiveSolution> saved_objective_and_solution_;

  bool submip;
  HighsInt submip_level;
  HighsInt max_submip_level;
  const HighsBasis* rootbasis;
  const HighsPseudocostInitialization* pscostinit;
  const HighsCliqueTable* clqtableinit;
  const HighsImplications* implicinit;

  std::unique_ptr<HighsMipSolverData> mipdata_;

  HighsMipAnalysis analysis_;

  MipRace mip_race_;

  HighsModelStatus termination_status_;
  HighsTerminator terminator_;

  void run();

  HighsInt numCol() const { return model_->num_col_; }

  HighsInt numRow() const { return model_->num_row_; }

  HighsInt numNonzero() const { return model_->a_matrix_.numNz(); }

  const double* colCost() const { return model_->col_cost_.data(); }

  double colCost(HighsInt col) const { return model_->col_cost_[col]; }

  const double* rowLower() const { return model_->row_lower_.data(); }

  double rowLower(HighsInt col) const { return model_->row_lower_[col]; }

  const double* rowUpper() const { return model_->row_upper_.data(); }

  double rowUpper(HighsInt col) const { return model_->row_upper_[col]; }

  const HighsVarType* variableType() const {
    return model_->integrality_.data();
  }

  HighsVarType variableType(HighsInt col) const {
    return model_->integrality_[col];
  }

  HighsMipSolver(HighsCallback& callback, const HighsOptions& options,
                 const HighsLp& lp, const HighsSolution& solution,
                 bool submip = false, HighsInt submip_level = 0);

  ~HighsMipSolver();

  void setModel(const HighsLp& model) {
    model_ = &model;
    solution_objective_ = kHighsInf;
  }

  mutable HighsTimer timer_;
  void cleanupSolve();

  void runPresolve(const HighsInt presolve_reduction_limit);
  const HighsLp& getPresolvedModel() const;
  HighsPresolveStatus getPresolveStatus() const;
  presolve::HighsPostsolveStack getPostsolveStack() const;

  void callbackGetCutPool() const;
  bool solutionFeasible(const HighsLp* lp, const std::vector<double>& col_value,
                        const std::vector<double>* pass_row_value,
                        double& bound_violation, double& row_violation,
                        double& integrality_violation, HighsCDouble& obj) const;
  std::vector<HighsModelStatus> initialiseRecord(HighsInt num_instance) const;
  void initialiseTerminator(HighsInt num_instance_ = 0,
                            HighsInt my_instance_ = kNoThreadInstance,
                            HighsModelStatus* record_ = nullptr);
  bool terminate() const {
    return this->termination_status_ != HighsModelStatus::kNotset;
  }
  HighsModelStatus terminationStatus() const {
    return this->termination_status_;
  }
};

#endif
