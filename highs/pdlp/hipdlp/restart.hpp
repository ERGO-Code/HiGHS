/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-09 14:54:26
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-11 15:19:54
 * @FilePath: /cupdlp-CPP/include/restart.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置
 * 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// restart.hpp
#ifndef RESTART_HPP
#define RESTART_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "scaling.hpp"
#include "solver_results.hpp"
#include "defs.hpp"

// Struct to communicate restart decisions
struct RestartInfo {
  bool should_restart = false;
  bool restart_to_average =
      false;  // If true, restart to average; otherwise, to current
};

class RestartScheme {
 public:
  RestartScheme() = default;

  void Initialize(const PrimalDualParams& params, const SolverResults& results);

  // Checks if a restart should be performed based on the chosen strategy
  RestartInfo Check(int current_iter, const SolverResults& current_results,
                    const SolverResults& average_results);

  int GetLastRestartIter() const { return last_restart_iter_; }

 private:
  // Computes a merit score for a given set of residuals
  double ComputeRestartScore(const SolverResults& results);

  RestartStrategy strategy_ = RestartStrategy::NO_RESTART;
  int fixed_restart_interval_ = 100;
  int last_restart_iter_ = 0;
  double beta_;

  // State for adaptive restart
  double last_restart_score_ = 1;
  double last_candidate_score_ = 1;
  double sufficient_decay_factor_ = 0.2;
  double necessary_decay_factor_ = 0.8;
};

#endif  // RESTART_HPP
