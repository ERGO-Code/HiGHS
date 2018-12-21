// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_CROSSOVER_H_
#define IPX_CROSSOVER_H_

#include <vector>
#include "basis.h"
#include "control.h"

namespace ipx {

class BasicSolution {
public:
    // Constructor initializes the object to the all-zero solution with slack
    // variables being basic and structural variables being superbasic. This is
    // a valid state which can be passed to Model::Postsolve*().
    BasicSolution(const Model& model);

    const Model& model() const { return model_; }
    Vector& x() { return x_; }
    Vector& y() { return y_; }
    Vector& z() { return z_; }
    std::vector<Int>& basic_statuses() { return basic_statuses_; }

    // Returns the maximum violation of AI*x=b.
    double presidual() const;

    // Returns the maximum violation of AI'y+z=c.
    double dresidual() const;

    // Returns the maximum violation of lb<=x<=ub.
    double pinfeas() const;

    // Returns the maximum sign violation of z. z is allowed to be
    // - positive only if x[j]==lb[j],
    // - negative only if x[j]==ub[j].
    // By this definition dual feasibility implies complementarity.
    double dinfeas() const;

    // Returns true if z[j]==0 || x[j]==lb[j] || x[j]==ub[j] for all j.
    bool complementary() const;

    // Calls Model::EvaluateBasicSolution().
    void EvaluatePostsolved(Info* info) const;

private:
    const Model& model_;
    Vector x_, y_, z_;
    std::vector<Int> basic_statuses_;
};

class Crossover {
public:
    explicit Crossover(const Control& control);

    // A crossover run consists of
    // (1) PushPhase(),
    // (2) computing the vertex solution defined by the basis, and
    // (3) setting the basic statuses.
    //
    // On return info->status_crossover is one of
    // IPX_STATUS_optimal           if (1) terminated successfully,
    // IPX_STATUS_time_limit        if (1) was interrupted,
    // IPX_STATUS_failed            if (1) failed.
    // In the latter case info->errflag is set. If (1) did not terminate
    // successfully, then (2), (3) are not executed and the basic statuses are
    // unchanged.
    void Run(const double* weights, BasicSolution* solution, Basis* basis,
             Info* info);

    // On entry solution->x, solution->y and solution->z must be initialized
    // such that solution->pinfeas() == 0.0 && solution->dinfeas() == 0.0.
    // (Otherwise an assertion will fail). On return solution->x, solution->y
    // and solution->z have been modified under the contition that
    // solution->pinfeas() == 0.0 && solution->dinfeas() == 0.0 still holds.
    // In the textbook method, solution->presidual() and solution->dresidual()
    // would have been unchanged.
    //
    // On successful return
    // - z[j]==0.0 for all basic variables j,
    // - x[j]==lb[j] || x[j]==ub[j] for all nonbasic variables j that have a
    //   finite bound.
    //
    // @weights if not NULL, then it must be an array of size # variables.
    //          Primal pushes are performed in decreasing order, and dual pushes
    //          in increasing order of the weights.
    // @solution x, y, z are set on input/output as specified above.
    //           basic_statuses is not accessed.
    // @basis must be a Basis object corresponding to the same Model object as
    //        solution. basis is updated on return.
    // @info errflag is set on return.
    void PushPhase(const double* weights, BasicSolution* solution, Basis* basis,
                   Info* info);

    // # primal pushes in last call to PushPhase().
    Int ppushes() const { return ppushes_; }

    // # dual pushes in last call to PushPhase().
    Int dpushes() const { return dpushes_; }

    // # pivots in last call to PushPhase().
    Int pivots() const { return pivots_; }

    // Runtime in last call to Run() or PushPhase().
    double time() const { return time_; }

private:
    static constexpr double kPivotZeroTol = 1e-5;

    // Bound is used by DualPushPhase() to indicates the status of x[j]:
    // NONE    if lb[j] < x[j] < ub[j]
    // LOWER   if x[j] = lb[j] < ub[j]
    // UPPER   if x[j] = ub[j] > lb[j]
    // EQ      if x[j] = lb[j] = ub[j]
    enum class Bound { NONE, LOWER, UPPER, EQ };

    // The following two methods are called from PushPhase() in that order.
    // We must execute the dual push phase first because it guarantees to
    // maintain complementarity (see the "wrong" pivots issue in the code).
    // After the dual push phase z[basic] = 0, so that an update to x[basic] in
    // the primal push phase cannot violate complementarity.
    void DualPushPhase(Info* info);
    void PrimalPushPhase(Info* info);

    // Computes the vertex solution defined by the basis and sets basic
    // statuses.
    void BuildBasicSolution();

    // Ratios tests called from DualPushPhase() and PrimalPushPhase().
    // They implement the two-pass ratio test that allows infeasibilities up to
    // feastol in order to choose a larger pivot.
    Int DualRatioTest(const Vector& z, const IndexedVector& row,
                      const std::vector<Bound>& atbound, double step,
                      double feastol);
    Int PrimalRatioTest(const Vector& xbasic, const IndexedVector& ftran,
                        const Vector& lbbasic, const Vector& ubbasic,
                        double step, double feastol, bool* block_at_lb);

    // Returns dual and primal superbasic variables in reverse order in which
    // they should be processed; i.e. the next variable to be pushed is at the
    // end of the vector.
    std::vector<Int> DualSuperbasics() const;
    std::vector<Int> PrimalSuperbasics() const;

    // Returns the # dual and primal superbasic variables.
    Int dpush_remain() const;
    Int ppush_remain() const;

    const Control& control_;
    BasicSolution* solution_{nullptr};
    Basis* basis_{nullptr};
    const double* weights_{nullptr};
    Int ppushes_{0};            // count # primal pushes
    Int dpushes_{0};            // count # dual pushes
    Int pivots_{0};             // count # pivots
    double time_{0.0};
};

}  // namespace ipx

#endif  // IPX_CROSSOVER_H_
