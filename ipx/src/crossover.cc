// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "crossover.h"
#include <algorithm>
#include <cassert>
#include "indexed_vector.h"
#include "time.h"
#include "utils.h"

namespace ipx {

BasicSolution::BasicSolution(const Model& model) : model_(model) {
        Int m = model.rows();
        Int n = model.cols();
        x_.resize(n+m);
        y_.resize(m);
        z_.resize(n+m);
        basic_statuses_.resize(n+m);
        for (Int j = 0; j < n; j++)
            basic_statuses_[j] = IPX_superbasic;
        for (Int i = 0; i < m; i++)
            basic_statuses_[n+i] = IPX_basic;
}

double BasicSolution::presidual() const {
    const Int m = model_.rows();
    const SparseMatrix& AIt = model_.AIt();
    const Vector& b = model_.b();

    double res = 0.0;
    for (Int i = 0; i < m; i++) {
        double r = b[i] - DotColumn(AIt, i, x_);
        res = std::max(res, std::abs(r));
    }
    return res;
}

double BasicSolution::dresidual() const {
    const SparseMatrix& AI = model_.AI();
    const Vector& c = model_.c();

    double res = 0.0;
    for (Int j = 0; j < c.size(); j++) {
        double r = c[j] - z_[j] - DotColumn(AI, j, y_);
        res = std::max(res, std::abs(r));
    }
    return res;
}

double BasicSolution::pinfeas() const {
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();

    double infeas = 0.0;
    for (Int j = 0; j < x_.size(); j++) {
        infeas = std::max(infeas, lb[j]-x_[j]);
        infeas = std::max(infeas, x_[j]-ub[j]);
    }
    return infeas;
}

double BasicSolution::dinfeas() const {
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();

    double infeas = 0.0;
    for (Int j = 0; j < x_.size(); j++) {
        if (x_[j] > lb[j])
            infeas = std::max(infeas, z_[j]);
        if (x_[j] < ub[j])
            infeas = std::max(infeas, -z_[j]);
    }
    return infeas;
}

bool BasicSolution::complementary() const {
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();

    for (Int j = 0; j < x_.size(); j++) {
        if (x_[j] != lb[j] && x_[j] != ub[j] && z_[j] != 0.0)
            return false;
    }
    return true;
}

void BasicSolution::EvaluatePostsolved(Info* info) const {
    model_.EvaluateBasicSolution(x_, y_, z_, basic_statuses_, info);
}

Crossover::Crossover(const Control& control) : control_(control) {}

void Crossover::Run(const double* weights, BasicSolution* solution,
                    Basis* basis, Info* info) {
    Timer timer;
    PushPhase(weights, solution, basis, info);
    if (info->errflag == IPX_ERROR_interrupt_time) {
        info->errflag = 0;
        info->status_crossover = IPX_STATUS_time_limit;
        return;
    }
    if (info->errflag) {
        info->status_crossover = IPX_STATUS_failed;
        return;
    }
    info->status_crossover = IPX_STATUS_optimal;
    BuildBasicSolution();
    control_.Debug()
        << Textline("Bound violation of basic solution:")
        << sci2(solution_->pinfeas()) << '\n'
        << Textline("Dual sign violation of basic solution:")
        << sci2(solution_->dinfeas()) << '\n';
    control_.Debug()
        << Textline("Minimum singular value of basis matrix:")
        << sci2(basis_->MinSingularValue()) << '\n';
    time_ = timer.Elapsed();
}

void Crossover::PushPhase(const double* weights, BasicSolution* solution,
                          Basis* basis, Info* info) {
    Timer timer;
    solution_ = solution;
    basis_ = basis;
    weights_ = weights;
    ppushes_ = 0;
    dpushes_ = 0;
    pivots_ = 0;
    assert(solution_->pinfeas() == 0.0);
    assert(solution_->dinfeas() == 0.0);
    info->errflag = 0;

    control_.Log()
        << Textline("Primal residual before push phase:")
        << sci2(solution_->presidual()) << '\n'
        << Textline("Dual residual before push phase:")
        << sci2(solution_->dresidual()) << '\n'
        << Textline("Number of primal pushes required:")
        << ppush_remain() << '\n'
        << Textline("Number of dual pushes required:")
        << dpush_remain() << '\n';

    // Should we do this here or in caller?
    basis_->UnfixVariables();
    basis_->UnfreeVariables();

    DualPushPhase(info);
    assert(solution_->pinfeas() == 0.0);
    assert(solution_->dinfeas() == 0.0);
    if (info->errflag != 0) {
        control_.Debug()
            << Textline("Minimum singular value of basis matrix:")
            << sci2(basis_->MinSingularValue()) << '\n';
        time_ = timer.Elapsed();
        return;
    }

    PrimalPushPhase(info);
    assert(solution_->pinfeas() == 0.0);
    assert(solution_->dinfeas() == 0.0);
    if (info->errflag != 0) {
        control_.Debug()
            << Textline("Minimum singular value of basis matrix:")
            << sci2(basis_->MinSingularValue()) << '\n';
        time_ = timer.Elapsed();
        return;
    }

    control_.Debug()
        << Textline("Primal residual after push phase:")
        << sci2(solution_->presidual()) << '\n'
        << Textline("Dual residual after push phase:")
        << sci2(solution_->dresidual()) << '\n';
    time_ = timer.Elapsed();
}

void Crossover::DualPushPhase(Info* info) {
    const Model& model = solution_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    Vector& x = solution_->x();
    Vector& y = solution_->y();
    Vector& z = solution_->z();
    IndexedVector btran(m), row(n+m);
    const double feastol = model.dualized() ?
        control_.pfeasibility_tol() : control_.dfeasibility_tol();
    std::vector<Int> superbasics = DualSuperbasics();

    std::vector<Bound> atbound(n+m);
    for (Int j = 0; j < n+m; j++) {
        if (lb[j] == ub[j])
            atbound[j] = Bound::EQ;
        else if (x[j] == lb[j])
            atbound[j] = Bound::LOWER;
        else if (x[j] == ub[j])
            atbound[j] = Bound::UPPER;
        else
            atbound[j] = Bound::NONE;
    }

    // A "wrong" pivot occurs at the dual push of jb if for a primal superbasic
    // variable jn the dual is moved away from zero. In this case we must pivot
    // immediately (without step) on jn, since otherwise complementarity would
    // be lost. The pivot is called "wrong" because it would not occur in exact
    // arithmetic and if the initial basis would minimize the number of
    // superbasics regarding any strictly complementary LP solution.
    Int num_wrong_pivots = 0;
    control_.ResetPrintInterval();

    while (!superbasics.empty()) {
        const Int jb = superbasics.back();
        assert(basis_->IsBasic(jb));
        assert(z[jb] != 0.0);
        assert(atbound[jb] == Bound::LOWER || atbound[jb] == Bound::UPPER);
        if ((info->errflag = control_.InterruptCheck()) != 0)
            break;

        basis_->TableauRow(jb, btran, row);
        double step = z[jb];
        Int jn = DualRatioTest(z, row, atbound, step, feastol);

        // If step was blocked, update basis and compute step size.
        if (jn >= 0) {
            assert(basis_->IsNonbasic(jn));
            double pivot = row[jn];
            assert(pivot);
            if (std::abs(pivot) < 1e-4)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot)) << '\n';
            bool exchanged;
            info->errflag = basis_->ExchangeIfStable(jb, jn, pivot, 1,
                                                     &exchanged);
            if (info->errflag)
                break;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            pivots_++;
            if (atbound[jn] == Bound::NONE) {
                step = 0.0;
                num_wrong_pivots++;
            } else {
                step = z[jn]/row[jn];
                if (atbound[jb] == Bound::UPPER)
                    assert(step <= 0.0);
                else
                    assert(step >= 0.0);
            }
        }
        // Update solution.
        if (step != 0.0) {
            auto update_y = [&](Int i, double x) {
                y[i] += step*x;
            };
            for_each_nonzero(btran, update_y);
            auto update_z = [&](Int j, double pivot) {
                z[j] -= step * pivot;
                if (atbound[j] == Bound::LOWER)
                    z[j] = std::max(z[j], 0.0);
                if (atbound[j] == Bound::UPPER)
                    z[j] = std::min(z[j], 0.0);
                if (atbound[j] == Bound::NONE)
                    z[j] = 0.0;
            };
            for_each_nonzero(row, update_z);
            z[jb] -= step;
        }
        if (jn >= 0)
            z[jn] = 0.0; // make clean
        else
            assert(z[jb] == 0.0);

        superbasics.pop_back();
        dpushes_++;
        control_.IntervalLog()
            << " " << Format(static_cast<Int>(superbasics.size()), 8)
            << " dual pushes remaining"
            << " (" << Format(pivots_, 7) << " pivots)\n";
    }
    control_.Debug()
        << Textline("Number of wrong pivots in dual push phase:")
        << num_wrong_pivots << '\n';
}

void Crossover::PrimalPushPhase(Info* info) {
    const Model& model = solution_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    Vector& x = solution_->x();
    Basis& basis = *basis_;
    IndexedVector ftran(m);
    const double feastol = model.dualized() ?
        control_.dfeasibility_tol() : control_.pfeasibility_tol();
    std::vector<Int> superbasics = PrimalSuperbasics();

    // Maintain a copy of primal basic variables and their bounds for faster
    // ratio test.
    Vector xbasic  = CopyBasic(x,  basis);
    Vector lbbasic = CopyBasic(lb, basis);
    Vector ubbasic = CopyBasic(ub, basis);

    control_.ResetPrintInterval();

    while (!superbasics.empty()) {
        const Int jn = superbasics.back();
        assert(basis.IsNonbasic(jn));
        assert(lb[jn] < x[jn] && x[jn] < ub[jn]);
        assert(std::isfinite(lb[jn]) || std::isfinite(ub[jn]));
        if ((info->errflag = control_.InterruptCheck()) != 0)
            break;

        // If the variable has two finite bounds, move to the nearer.
        bool move_to_lb;
        if (std::isfinite(lb[jn]) && std::isfinite(ub[jn]))
            move_to_lb = x[jn]-lb[jn] <= ub[jn]-x[jn];
        else
            move_to_lb = std::isfinite(lb[jn]);

        // A full step is such that x[jn]-step is at its bound.
        double step = move_to_lb ? x[jn]-lb[jn] : x[jn]-ub[jn];

        basis.SolveForUpdate(jn, ftran);
        bool block_at_lb;
        Int pblock = PrimalRatioTest(xbasic, ftran, lbbasic, ubbasic, step,
                                     feastol, &block_at_lb);
        Int jb = pblock >= 0 ? basis[pblock] : -1;

        // If step was blocked, update basis and compute step size.
        if (pblock >= 0) {
            double pivot = ftran[pblock];
            assert(pivot != 0.0);
            if (std::abs(pivot) < 1e-4)
                control_.Debug(3)
                    << " |pivot| = " << sci2(std::abs(pivot)) << '\n';
            bool exchanged;
            info->errflag = basis.ExchangeIfStable(jb, jn, pivot, -1,
                                                   &exchanged);
            if (info->errflag)
                break;
            if (!exchanged)     // factorization was unstable, try again
                continue;
            pivots_++;
            if (block_at_lb)
                step = (lbbasic[pblock]-xbasic[pblock]) / ftran[pblock];
            else
                step = (ubbasic[pblock]-xbasic[pblock]) / ftran[pblock];
        }
        // Update solution.
        if (step != 0.0) {
            auto update = [&](Int p, double pivot) {
                xbasic[p] += step * pivot;
                xbasic[p] = std::max(xbasic[p], lbbasic[p]);
                xbasic[p] = std::min(xbasic[p], ubbasic[p]);
            };
            for_each_nonzero(ftran, update);
            x[jn] -= step;
        }
        if (pblock >= 0) {
            xbasic[pblock] = x[jn];
            lbbasic[pblock] = lb[jn];
            ubbasic[pblock] = ub[jn];
            x[jb] = block_at_lb ? lb[jb] : ub[jb];
            assert(std::isfinite(x[jb]));
        } else {
            x[jn] = move_to_lb ? lb[jn] : ub[jn]; // make clean
            assert(std::isfinite(x[jn]));
        }

        superbasics.pop_back();
        ppushes_++;
        control_.IntervalLog()
            << " " << Format(static_cast<Int>(superbasics.size()), 8)
            << " primal pushes remaining"
            << " (" << Format(pivots_, 7) << " pivots)\n";
    }
    for (Int p = 0; p < m; p++)
        x[basis[p]] = xbasic[p];
}

void Crossover::BuildBasicSolution() {
    const Model& model = solution_->model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& b = model.b();
    const Vector& c = model.c();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    const SparseMatrix& AI = model.AI();
    Vector& x = solution_->x();
    Vector& y = solution_->y();
    Vector& z = solution_->z();
    std::vector<Int>& basic_statuses = solution_->basic_statuses();
    const Basis& basis = *basis_;

    // Compute x[basic] so that AI*x=b.
    Vector work = b;
    for (Int j = 0; j < n+m; j++)
        if (basis.IsNonbasic(j))
            ScatterColumn(AI, j, -x[j], work);
    basis.SolveDense(work, work, 'N');
    for (Int p = 0; p < m; p++)
        x[basis[p]] = work[p];

    // Compute y and z so that AI'y+z=c.
    for (Int p = 0; p < m; p++)
        work[p] = c[basis[p]];
    basis.SolveDense(work, y, 'T');
    for (Int j = 0; j < n+m; j++) {
        if (basis.IsNonbasic(j))
            z[j] = c[j] - DotColumn(AI, j, y);
        else
            z[j] = 0.0;
    }

    for (Int j = 0; j < n+m; j++) {
        if (basis_->IsBasic(j)) {
            basic_statuses[j] = IPX_basic;
        } else {
            if (lb[j] == ub[j])
                basic_statuses[j] = z[j] >= 0.0 ?
                    IPX_nonbasic_lb : IPX_nonbasic_ub;
            else if (x[j] == lb[j])
                basic_statuses[j] = IPX_nonbasic_lb;
            else if (x[j] == ub[j])
                basic_statuses[j] = IPX_nonbasic_ub;
            else
                basic_statuses[j] = IPX_superbasic;
        }
    }
}

Int Crossover::DualRatioTest(const Vector& z, const IndexedVector& row,
                             const std::vector<Bound>& atbound, double step,
                             double feastol) {
    Int jblock = -1;            // return value

    // First pass: determine maximum step size exploiting feasibility tol.
    // If wrong pivots occur, choose the maximum one.
    double max_wrong_pivot = 0.0;
    auto update_step = [&](Int j, double pivot) {
        if (std::abs(pivot) > kPivotZeroTol) {
            if (atbound[j] == Bound::NONE && std::abs(pivot) > max_wrong_pivot){
                assert(z[j] == 0.0);
                jblock = j;
                max_wrong_pivot = std::abs(pivot);
                step = 0.0;
            }
            else if (atbound[j] == Bound::LOWER && z[j]-step*pivot < -feastol) {
                step = (z[j]+feastol) / pivot;
                jblock = j;
                assert(z[j] >= 0.0);
                assert(step*pivot > 0.0);
            }
            else if (atbound[j] == Bound::UPPER && z[j]-step*pivot > feastol) {
                step = (z[j]-feastol) / pivot;
                jblock = j;
                assert(z[j] <= 0.0);
                assert(step*pivot < 0.0);
            }
        }
    };
    for_each_nonzero(row, update_step);

    // If step was not block or a wrong pivot was found, we are done.
    if (jblock < 0 || atbound[jblock] == Bound::NONE)
        return jblock;

    // Second pass: choose maximum pivot among all that block within step.
    jblock = -1;
    double max_pivot = kPivotZeroTol;
    auto update_max = [&](Int j, double pivot) {
        if (std::abs(pivot) > max_pivot &&
            std::abs(z[j]/pivot) <= std::abs(step)) {
            if (atbound[j] == Bound::LOWER && step*pivot > 0.0) {
                jblock = j;
                max_pivot = std::abs(pivot);
            }
            if (atbound[j] == Bound::UPPER && step*pivot < 0.0) {
                jblock = j;
                max_pivot = std::abs(pivot);
            }
        }
    };
    for_each_nonzero(row, update_max);
    assert(jblock >= 0);
    return jblock;
}

Int Crossover::PrimalRatioTest(const Vector& xbasic, const IndexedVector& ftran,
                               const Vector& lbbasic, const Vector& ubbasic,
                               double step, double feastol, bool* block_at_lb) {
    Int pblock = -1;            // return value
    *block_at_lb = true;

    // First pass: determine maximum step size exploiting feasibility tol.
    auto update_step = [&](Int p, double pivot) {
        if (std::abs(pivot) > kPivotZeroTol) {
            // test block at lower bound
            if (xbasic[p] + step*pivot < lbbasic[p]-feastol) {
                step = (lbbasic[p]-xbasic[p]-feastol) / pivot;
                pblock = p;
                *block_at_lb = true;
            }
            // test block at upper bound
            if (xbasic[p] + step*pivot > ubbasic[p]+feastol) {
                step = (ubbasic[p]-xbasic[p]+feastol) / pivot;
                pblock = p;
                *block_at_lb = false;
            }
        }
    };
    for_each_nonzero(ftran, update_step);

    // If the step was not blocked, we are done.
    if (pblock < 0)
        return pblock;

    // Second pass: choose maximum pivot among all that block within step.
    pblock = -1;
    double max_pivot = kPivotZeroTol;
    auto update_max = [&](Int p, double pivot) {
        if (std::abs(pivot) > max_pivot) {
            // test block at lower bound
            if (step*pivot < 0.0) {
                double step_p = (lbbasic[p]-xbasic[p]) / pivot;
                if (std::abs(step_p) <= std::abs(step)) {
                    pblock = p;
                    *block_at_lb = true;
                    max_pivot = std::abs(pivot);
                }
            }
            // test block at upper bound
            if (step*pivot > 0.0) {
                double step_p = (ubbasic[p]-xbasic[p]) / pivot;
                if (std::abs(step_p) <= std::abs(step)) {
                    pblock = p;
                    *block_at_lb = false;
                    max_pivot = std::abs(pivot);
                }
            }
        }
    };
    for_each_nonzero(ftran, update_max);
    assert(pblock >= 0);
    return pblock;
}

std::vector<Int> Crossover::DualSuperbasics() const {
    const Vector& z = solution_->z();

    std::vector<Int> perm = Sortperm(z.size(), weights_, true);
    std::vector<Int> superbasics;
    for (Int j : perm) {
        if (basis_->IsBasic(j) && z[j] != 0.0)
            superbasics.push_back(j);
    }
    return superbasics;
}

std::vector<Int> Crossover::PrimalSuperbasics() const {
    const Model& model = solution_->model();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    const Vector& x = solution_->x();

    std::vector<Int> perm = Sortperm(x.size(), weights_, false);
    std::vector<Int> superbasics;
    for (Int j : perm) {
        if (basis_->IsNonbasic(j)) {
            if ((std::isfinite(lb[j]) || std::isfinite(ub[j])) &&
                (x[j] != lb[j] && x[j] != ub[j])) {
                superbasics.push_back(j);
            }
        }
    }
    return superbasics;
}

Int Crossover::dpush_remain() const {
    const Vector& z = solution_->z();

    Int cnt = 0;
    for (Int j = 0; j < z.size(); j++) {
        if (basis_->IsBasic(j) && z[j] != 0.0)
            cnt++;
    }
    return cnt;
}

Int Crossover::ppush_remain() const {
    const Model& model = solution_->model();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    const Vector& x = solution_->x();

    Int cnt = 0;
    for (Int j = 0; j < x.size(); j++) {
        if (basis_->IsNonbasic(j))
            if ((std::isfinite(lb[j]) || std::isfinite(ub[j])) &&
                (x[j] != lb[j] && x[j] != ub[j]))
                cnt++;
    }
    return cnt;
}

}  // namespace ipx
