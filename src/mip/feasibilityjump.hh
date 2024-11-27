#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <cassert>
#include <algorithm>

#define FJ_LOG_PREFIX "Feasibility Jump: "

namespace external_feasibilityjump {

    enum RowType {
      Equal,
      Lte,
      Gte,
    };

    enum VarType { Continuous, Integer };

    enum CallbackControlFlow {
      Terminate,
      Continue,
    };

    struct FJStatus {
      int totalEffort;
      int effortSinceLastImprovement;
      int numVars;
      double solutionObjectiveValue;
      double* solution;
    };

    const double violationTolerance = 1.0e-5;
    const double equalityTolerance = 1.0e-5;

    // Measures if two doubles are equal within a tolerance of 1.0e-5.
    bool eq(double a, double b) { return fabs(a - b) < equalityTolerance; }

    struct IdxCoeff {
      uint32_t idx;
      double coeff;

      IdxCoeff(uint32_t idx, double coeff) : idx(idx), coeff(coeff) {}
    };

    struct Var {
      VarType vartype;
      double lb;
      double ub;
      double objectiveCoeff;
      std::vector<IdxCoeff> coeffs;
    };

    struct Constraint {
      RowType sense;
      double rhs;
      std::vector<IdxCoeff> coeffs;
      double weight;
      double incumbentLhs;
      int32_t violatedIdx;

      // Computes the constraint's contribution to the feasibility score:
      // If the constraint is satisfied by the given LHS value, returns 0.
      // If the constraint is violated by the given LHS value, returns -|lhs-rhs|.
      double score(double lhs) {
        if (sense == RowType::Equal)
          return -fabs(lhs - rhs);
        else if (sense == RowType::Lte)
          return -(std::max(0., lhs - rhs));
        else
          return -(std::max(0., rhs - lhs));
      }
    };

    // A potential new value for a varaiable, including its score.
    struct Move {
      double value;
      double score;

      static Move undef() {
        Move move;
        move.value = NAN;
        move.score = -std::numeric_limits<double>::infinity();
        return move;
      }
    };

    // Represents a modification of the LHS in a constraint, for a specific
    // variable/constraint combination.The `modifyMove` function below is used to
    // update the score of a `Move` to reflect the LHS modification.
    struct LhsModification {
      uint32_t varIdx;
      uint32_t constraintIdx;
      double coeff;
      double oldLhs;
      double newLhs;
    };

    // Stores the MIP problem, an incumbent assignment, and the set of constraints
    // that are violated in the current incumbent assignment. This set is maintained
    // when changes are given to the incumbent assignment using `setValue`.
    struct Problem {
      std::vector<Var> vars;
      std::vector<Constraint> constraints;
      std::vector<double> incumbentAssignment;
      std::vector<uint32_t> violatedConstraints;
      bool usedRelaxContinuous = false;

      size_t nNonzeros;
      double incumbentObjective = NAN;

      int addVar(VarType vartype, double lb, double ub, double objCoeff) {
        auto idx = vars.size();
        Var var;
        var.vartype = vartype;
        var.lb = lb;
        var.ub = ub;
        var.objectiveCoeff = objCoeff;
        vars.push_back(var);
        incumbentAssignment.push_back(lb);
        return idx;
      }

      int addConstraint(RowType sense, double rhs, int numCoeffs, int* rowVarIdxs,
                        double* rowCoeffs, int relax_continuous) {
        if (relax_continuous) usedRelaxContinuous = true;

        // If we are relaxing continuous variables, an equality needs to be split
        // into Gte and Lte.
        if (relax_continuous > 0 && sense == RowType::Equal)
          if (std::any_of(rowVarIdxs, rowVarIdxs + numCoeffs, [&](double varIdx) {
                return vars[varIdx].vartype == VarType::Continuous;
              })) {
            addConstraint(RowType::Gte, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                          relax_continuous);
            addConstraint(RowType::Lte, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                          relax_continuous);
            return INT_MAX;
          }

        std::vector<IdxCoeff> coeffs;
        for (int i = 0; i < numCoeffs; i += 1) {
          if (relax_continuous > 0 &&
              vars[rowVarIdxs[i]].vartype == VarType::Continuous) {
            if (sense == RowType::Lte) {
              if (rowCoeffs[i] >= 0.)
                rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
              else
                rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
            } else if (sense == RowType::Gte) {
              if (rowCoeffs[i] >= 0.)
                rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
              else
                rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
            } else
              return INT_MIN;
          } else
            coeffs.emplace_back(rowVarIdxs[i], rowCoeffs[i]);
        }

        if (coeffs.empty()) {
          bool ok;
          if (sense == RowType::Lte)
            ok = 0 <= rhs + equalityTolerance;
          else if (sense == RowType::Gte)
            ok = 0 + equalityTolerance >= rhs;
          else
            ok = eq(0, rhs);

          return ok ? INT_MAX : INT_MIN;
        }

        int newConstraintIdx = constraints.size();
        for (auto& c : coeffs) {
          vars[c.idx].coeffs.emplace_back(newConstraintIdx, c.coeff);
        }

        nNonzeros += coeffs.size();
        Constraint newConstraint;
        newConstraint.coeffs = coeffs;
        newConstraint.incumbentLhs = NAN;
        newConstraint.violatedIdx = -1;
        newConstraint.rhs = rhs;
        newConstraint.sense = sense;
        newConstraint.weight = 1.0;

        constraints.push_back(newConstraint);
        return newConstraintIdx;
      }

      void resetIncumbent(double* initialValues) {
        // Set the initial values, if given.

        if (initialValues)
          for (size_t i = 0; i < vars.size(); i += 1)
            incumbentAssignment[i] = initialValues[i];
        // std::copy(initialValues, initialValues + vars.size(),
        // incumbentAssignment);

        // Reset the incumbent objective.
        incumbentObjective = 0;
        for (size_t i = 0; i < vars.size(); i += 1)
          incumbentObjective += vars[i].objectiveCoeff * incumbentAssignment[i];

        // Reset the constraint LHSs and the violatedConstraints list.
        violatedConstraints.clear();
        for (size_t cIdx = 0; cIdx < constraints.size(); cIdx += 1) {
          Constraint& cstr = constraints[cIdx];

          cstr.incumbentLhs = 0.0;
          for (auto& vc : cstr.coeffs)
            cstr.incumbentLhs += vc.coeff * incumbentAssignment[vc.idx];

          if (cstr.score(cstr.incumbentLhs) < -violationTolerance) {
            cstr.violatedIdx = violatedConstraints.size();
            violatedConstraints.push_back(cIdx);
          } else
            cstr.violatedIdx = -1;
        }
      }

      // Updates a variable assignment for `varIdx` to `newValue`.
      // Takes a function parameter f that receives a LhsModification
      // for every variable/constraint combination (except for `varIdx` itself)
      // where the LHS of the constraint has changed.
      template <typename F>
      size_t setValue(uint32_t varIdx, double newValue, F f) {
        size_t dt = 0;
        double oldValue = incumbentAssignment[varIdx];
        double delta = (newValue - oldValue);
        incumbentAssignment[varIdx] = newValue;
        incumbentObjective += vars[varIdx].objectiveCoeff * delta;
        // printf("Setting v%d to from %g to value %g\n", varIdx, oldValue,
        // newValue);

        // Update the LHSs of all involved constraints.
        for (auto& cstrCoeff : vars[varIdx].coeffs) {
          double oldLhs = constraints[cstrCoeff.idx].incumbentLhs;
          double newLhs = oldLhs + cstrCoeff.coeff * delta;
          constraints[cstrCoeff.idx].incumbentLhs = newLhs;
          double newCost = constraints[cstrCoeff.idx].score(newLhs);

          // Add/remove from the violatedConstraints list.
          if (newCost < -violationTolerance &&
              constraints[cstrCoeff.idx].violatedIdx == -1) {
            // Became violated.
            constraints[cstrCoeff.idx].violatedIdx = violatedConstraints.size();
            violatedConstraints.push_back(cstrCoeff.idx);
          }
          if (newCost >= -violationTolerance &&
              constraints[cstrCoeff.idx].violatedIdx != -1) {
            // Became satisfied.
            auto lastViolatedIdx = violatedConstraints.size() - 1;
            auto lastConstraintIdx = violatedConstraints[lastViolatedIdx];
            auto thisViolatedIdx = constraints[cstrCoeff.idx].violatedIdx;
            std::swap(violatedConstraints[thisViolatedIdx],
                      violatedConstraints[lastViolatedIdx]);
            constraints[lastConstraintIdx].violatedIdx = thisViolatedIdx;
            constraints[cstrCoeff.idx].violatedIdx = -1;
            violatedConstraints.pop_back();
          }

          // Now, report the changes in LHS for other variables.
          dt += constraints[cstrCoeff.idx].coeffs.size();
          for (auto& varCoeff : constraints[cstrCoeff.idx].coeffs) {
            if (varCoeff.idx != varIdx) {
              LhsModification m;
              m.varIdx = varCoeff.idx;
              m.constraintIdx = cstrCoeff.idx;
              m.coeff = varCoeff.coeff;
              m.oldLhs = oldLhs;
              m.newLhs = newLhs;
              f(m);
            }
          }
        }

        return dt;
      }
    };

    void modifyMove(LhsModification mod, Problem& problem, Move& move) {
      Constraint& c = problem.constraints[mod.constraintIdx];
      auto incumbent = problem.incumbentAssignment[mod.varIdx];
      double oldModifiedLhs = mod.oldLhs + mod.coeff * (move.value - incumbent);
      double oldScoreTerm =
          c.weight * (c.score(oldModifiedLhs) - c.score(mod.oldLhs));
      double newModifiedLhs = mod.newLhs + mod.coeff * (move.value - incumbent);
      double newScoreTerm =
          c.weight * (c.score(newModifiedLhs) - c.score(mod.newLhs));
      move.score += newScoreTerm - oldScoreTerm;
    }

    // Stores current moves and computes updated jump values for
    // the "Jump" move type.
    class JumpMove {
      std::vector<Move> moves;
      std::vector<std::pair<double, double>> bestShiftBuffer;

     public:
      void init(Problem& problem) { moves.resize(problem.vars.size()); }

      template <typename F>
      void forEachVarMove(int32_t varIdx, F f) {
        f(moves[varIdx]);
      }

      void updateValue(Problem& problem, uint32_t varIdx) {
        bestShiftBuffer.clear();
        auto varIncumbentValue = problem.incumbentAssignment[varIdx];
        double currentValue = problem.vars[varIdx].lb;
        double currentScore = 0.0;
        double currentSlope = 0.0;

        // printf(" updatevalue lb %g ub %g numcells %d\n",
        //        problem.vars[varIdx].lb,
        //        problem.vars[varIdx].ub, problem.vars[varIdx].coeffs.size());

        for (auto& cell : problem.vars[varIdx].coeffs) {
          auto& constraint = problem.constraints[cell.idx];

          std::vector<std::pair<double, double>> constraintBounds;
          if (constraint.sense == RowType::Lte)
            constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(),
                                          constraint.rhs);
          else if (constraint.sense == RowType::Gte)
            constraintBounds.emplace_back(constraint.rhs,
                                          std::numeric_limits<double>::infinity());
          else {
            constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(),
                                          constraint.rhs);
            constraintBounds.emplace_back(constraint.rhs, constraint.rhs);
            constraintBounds.emplace_back(constraint.rhs,
                                          std::numeric_limits<double>::infinity());
          }

          for (auto& bound : constraintBounds) {
            double residualIncumbent =
                constraint.incumbentLhs - cell.coeff * varIncumbentValue;

            std::pair<double, double> validRange = {
                ((1.0 / cell.coeff) * (bound.first - residualIncumbent)),
                ((1.0 / cell.coeff) * (bound.second - residualIncumbent)),
            };

            if (problem.vars[varIdx].vartype == VarType::Integer)
              validRange = {
                  std::ceil(validRange.first - equalityTolerance),
                  std::floor(validRange.second + equalityTolerance),
              };

            if (validRange.first > validRange.second) continue;

            if (validRange.first > currentValue) {
              currentScore += constraint.weight * (validRange.first - currentValue);
              currentSlope -= constraint.weight;
              if (validRange.first < problem.vars[varIdx].ub)
                bestShiftBuffer.emplace_back(validRange.first, constraint.weight);
            }

            if (validRange.second <= currentValue) {
              currentScore +=
                  constraint.weight * (validRange.second - currentValue);
              currentSlope += constraint.weight;
            } else if (validRange.second < problem.vars[varIdx].ub)
              bestShiftBuffer.emplace_back(validRange.second, constraint.weight);
          }
        }

        bestShiftBuffer.emplace_back(problem.vars[varIdx].lb, 0);
        bestShiftBuffer.emplace_back(problem.vars[varIdx].ub, 0);
        std::sort(bestShiftBuffer.begin(), bestShiftBuffer.end());

        double bestScore = currentScore;
        double bestValue = currentValue;
        // printf("evaluating best shift buffer size %d \n",
        // bestShiftBuffer.size());
        for (auto& item : bestShiftBuffer) {
          currentScore += (item.first - currentValue) * currentSlope;
          currentSlope += item.second;
          currentValue = item.first;

          // printf("bestshift cscore %g cslope %g cval %g bestval %g bestscore
          // %g\n", currentScore,currentSlope, currentValue, bestScore, bestValue
          // );

          if (eq(bestValue, problem.incumbentAssignment[varIdx]) ||
              (!eq(currentValue, problem.incumbentAssignment[varIdx]) &&
               currentScore < bestScore)) {
            bestScore = currentScore;
            bestValue = currentValue;
          }

          // Slope is always increasing, so if we have a valid value, we can quit
          // as soon as the slope turns nonnegative, since we must already have
          // visited the minimum.
          if (!eq(bestValue, problem.incumbentAssignment[varIdx]) &&
              currentSlope >= 0.)
            break;
        }

        // printf("Setting jump for %d to from %g to %g\n", varIdx,
        // problem.incumbentAssignment[varIdx], moves[varIdx].value);
        moves[varIdx].value = bestValue;
      }
    };

    class FeasibilityJumpSolver {
      int verbosity;
      Problem problem;
      JumpMove jumpMove;

      std::vector<uint32_t> goodVarsSet;
      std::vector<int32_t> goodVarsSetIdx;

      std::mt19937 rng;

      double bestObjective = std::numeric_limits<double>::infinity();
      double objectiveWeight = 0.0;
      size_t bestViolationScore = SIZE_MAX;
      size_t effortAtLastCallback = 0;
      size_t effortAtLastImprovement = 0;
      size_t totalEffort = 0;

      double weightUpdateDecay;
      double weightUpdateIncrement = 1.0;

      size_t nBumps;

      // The probability of choosing a random positive-score variable.
      const double randomVarProbability = 0.001;

      // The probability of choosing a variable using a random constraint's
      // non-zero coefficient after updating weights.

      const double randomCellProbability = 0.01;

      // The number of moves to evaluate, if there are many positive-score
      // variables available.
      const size_t maxMovesToEvaluate = 25;

     public:
      FeasibilityJumpSolver(int seed = 0, int _verbosity = 0,
                            double _weightUpdateDecay = 1.0) {
        verbosity = _verbosity;
        weightUpdateDecay = _weightUpdateDecay;
        rng = std::mt19937(seed);
      }

      int addVar(VarType vartype, double lb, double ub, double objCoeff) {
        goodVarsSetIdx.push_back(-1);
        return problem.addVar(vartype, lb, ub, objCoeff);
      }

      int addConstraint(RowType sense, double rhs, int numCoeffs, int* rowVarIdxs,
                        double* rowCoeffs, int relax_continuous) {
        return problem.addConstraint(sense, rhs, numCoeffs, rowVarIdxs, rowCoeffs,
                                     relax_continuous);
      }

      int solve(double* initialValues,
                std::function<CallbackControlFlow(FJStatus)> callback) {
        assert(callback);
        if (verbosity >= 1)
          printf(FJ_LOG_PREFIX
                 "starting solve. weightUpdateDecay=%g, relaxContinuous=%d  \n",
                 weightUpdateDecay, problem.usedRelaxContinuous);

        init(initialValues);

        for (int step = 0; step < INT_MAX; step += 1) {
          if (user_terminate(callback, nullptr)) break;

          if (step % 100000 == 0) {
            if (verbosity >= 1)
              printf(FJ_LOG_PREFIX "step %d viol %zd good %zd bumps %zd\n", step,
                     problem.violatedConstraints.size(), goodVarsSet.size(),
                     nBumps);
          }

          if (problem.violatedConstraints.size() < bestViolationScore) {
            effortAtLastImprovement = totalEffort;
            bestViolationScore = problem.violatedConstraints.size();
          }

          if (problem.violatedConstraints.empty() &&
              problem.incumbentObjective < bestObjective) {
            effortAtLastImprovement = totalEffort;
            bestObjective = problem.incumbentObjective;
            if (user_terminate(callback, problem.incumbentAssignment.data())) break;
          }

          if (problem.vars.size() == 0) break;

          uint32_t var = selectVariable();
          doVariableMove(var);
        }

        return 0;
      }

     private:
      void init(double* initialValues) {
        problem.resetIncumbent(initialValues);
        jumpMove.init(problem);
        totalEffort += problem.nNonzeros;

        // Reset the variable scores.
        goodVarsSet.clear();
        for (size_t i = 0; i < problem.vars.size(); i += 1) resetMoves(i);
      }

      uint32_t selectVariable() {
        if (!goodVarsSet.empty()) {
          if (std::uniform_real_distribution<double>(0., 1.)(rng) <
              randomVarProbability)
            return goodVarsSet[rng() % goodVarsSet.size()];

          auto sampleSize = std::min(maxMovesToEvaluate, goodVarsSet.size());
          totalEffort += sampleSize;

          double bestScore = -std::numeric_limits<double>::infinity();
          uint32_t bestVar = UINT_MAX;
          for (size_t i = 0; i < sampleSize; i++) {
            auto setidx = rng() % goodVarsSet.size();
            auto varIdx = goodVarsSet[setidx];
            // assert(goodVarsSetIdx[varIdx] >= 0 && goodVarsSetIdx[varIdx] ==
            // setidx);
            Move move = bestMove(varIdx);
            // assert(move.score > equalityTolerance);
            if (move.score > bestScore) {
              bestScore = move.score;
              bestVar = varIdx;
            }
          }
          assert(bestVar != UINT_MAX);
          return bestVar;
        }

        // Local minimum, update weights.
        updateWeights();

        if (!problem.violatedConstraints.empty()) {
          size_t cstrIdx =
              problem
                  .violatedConstraints[rng() % problem.violatedConstraints.size()];
          auto& constraint = problem.constraints[cstrIdx];

          if (std::uniform_real_distribution<double>(0., 1.)(rng) <
              randomCellProbability)
            return constraint.coeffs[rng() % constraint.coeffs.size()].idx;

          double bestScore = -std::numeric_limits<double>::infinity();
          uint32_t bestVarIdx = UINT_MAX;

          for (auto& cell : constraint.coeffs) {
            Move move = bestMove(cell.idx);
            if (move.score > bestScore) {
              bestScore = move.score;
              bestVarIdx = cell.idx;
            }
          }
          return bestVarIdx;
        }

        // Fallback to random choice.
        return rng() % problem.vars.size();
      }

      void updateWeights() {
        if (verbosity >= 2) printf(FJ_LOG_PREFIX "Reached a local minimum.\n");
        nBumps += 1;
        bool rescaleAllWeights = false;
        size_t dt = 0;

        if (problem.violatedConstraints.empty()) {
          objectiveWeight += weightUpdateIncrement;
          if (objectiveWeight > 1.0e20) rescaleAllWeights = true;

          dt += problem.vars.size();
          for (size_t varIdx = 0; varIdx < problem.vars.size(); varIdx += 1)
            forEachMove(varIdx, [&](Move& move) {
              move.score += weightUpdateIncrement *
                            problem.vars[varIdx].objectiveCoeff *
                            (move.value - problem.incumbentAssignment[varIdx]);
            });
        } else {
          for (auto& cIdx : problem.violatedConstraints) {
            auto& constraint = problem.constraints[cIdx];
            constraint.weight += weightUpdateIncrement;
            if (constraint.weight > 1.0e20) rescaleAllWeights = true;

            dt += constraint.coeffs.size();
            for (auto& cell : constraint.coeffs) {
              forEachMove(cell.idx, [&](Move& move) {
                double candidateLhs =
                    constraint.incumbentLhs +
                    cell.coeff *
                        (move.value - problem.incumbentAssignment[cell.idx]);
                double diff = weightUpdateIncrement *
                              (constraint.score(candidateLhs) -
                               constraint.score(constraint.incumbentLhs));
                move.score += diff;
              });

              updateGoodMoves(cell.idx);
            }
          }
        }

        weightUpdateIncrement /= weightUpdateDecay;
        if (rescaleAllWeights) {
          weightUpdateIncrement *= 1.0e-20;
          objectiveWeight *= 1.0e-20;

          for (auto& c : problem.constraints) c.weight *= 1.0e-20;
          dt += problem.constraints.size();

          for (size_t i = 0; i < problem.vars.size(); i += 1) resetMoves(i);
        }

        totalEffort += dt;
      }

      Move bestMove(uint32_t varIdx) {
        Move best = Move::undef();
        forEachMove(varIdx, [&](Move& move) {
          if (move.score > best.score) best = move;
        });

        return best;
      }

      void doVariableMove(uint32_t varIdx) {
        // First, we get the best move for the variable;
        auto m = bestMove(varIdx);
        auto newValue = m.value;
        // assert(!isnan(newValue));

        // Update the incumbent solution.
        // printf("Setting var %d from %g to %g for a score of %g\n", varIdx,
        // oldValue, newValue, m.score);

        totalEffort += problem.setValue(varIdx, newValue, [&](LhsModification mod) {
          forEachMove(mod.varIdx, [&](Move& m) { modifyMove(mod, problem, m); });
          updateGoodMoves(mod.varIdx);
        });

        resetMoves(varIdx);
      }

      void updateGoodMoves(int32_t varIdx) {
        bool anyGoodMoves = bestMove(varIdx).score > 0.;
        if (anyGoodMoves && goodVarsSetIdx[varIdx] == -1) {
          // Became good, add to good set.
          goodVarsSetIdx[varIdx] = goodVarsSet.size();
          goodVarsSet.push_back(varIdx);
        } else if (!anyGoodMoves && goodVarsSetIdx[varIdx] != -1) {
          // Became bad, remove from good set.
          auto lastSetIdx = goodVarsSet.size() - 1;
          auto lastVarIdx = goodVarsSet[lastSetIdx];
          auto thisSetIdx = goodVarsSetIdx[varIdx];
          std::swap(goodVarsSet[thisSetIdx], goodVarsSet[lastSetIdx]);
          goodVarsSetIdx[lastVarIdx] = thisSetIdx;
          goodVarsSetIdx[varIdx] = -1;
          goodVarsSet.pop_back();
        }
      }

      template <typename F>
      void forEachMove(int32_t varIdx, F f) {
        jumpMove.forEachVarMove(varIdx, f);

        // TODO: here, we can add more move types.
        // upDownMove.forEachVarMove(varIdx, f);
      }

      void resetMoves(uint32_t varIdx) {
        totalEffort += problem.vars[varIdx].coeffs.size();
        jumpMove.updateValue(problem, varIdx);

        forEachMove(varIdx, [&](Move& move) {
          move.score = 0.0;
          move.score += objectiveWeight * problem.vars[varIdx].objectiveCoeff *
                        (move.value - problem.incumbentAssignment[varIdx]);

          for (auto& cell : problem.vars[varIdx].coeffs) {
            auto& constraint = problem.constraints[cell.idx];
            auto candidateLhs =
                constraint.incumbentLhs +
                cell.coeff * (move.value - problem.incumbentAssignment[varIdx]);
            move.score +=
                constraint.weight * (constraint.score(candidateLhs) -
                                     constraint.score(constraint.incumbentLhs));
          }
        });

        updateGoodMoves(varIdx);
      }

      bool user_terminate(std::function<CallbackControlFlow(FJStatus)> callback,
                          double* solution) {
        const int CALLBACK_EFFORT = 500000;
        if (solution != nullptr ||
            totalEffort - effortAtLastCallback > CALLBACK_EFFORT) {
          if (verbosity >= 2) printf(FJ_LOG_PREFIX "calling user termination.\n");
          effortAtLastCallback = totalEffort;

          FJStatus status;
          status.totalEffort = totalEffort;
          status.effortSinceLastImprovement = totalEffort - effortAtLastImprovement;

          status.solution = solution;
          status.numVars = problem.vars.size();
          status.solutionObjectiveValue = problem.incumbentObjective;

          auto result = callback(status);
          if (result == CallbackControlFlow::Terminate) {
            if (verbosity >= 2) printf(FJ_LOG_PREFIX "quitting.\n");
            return true;
          }
        }
        return false;
      }
    };

}  // namespace external_feasiblityjump;
