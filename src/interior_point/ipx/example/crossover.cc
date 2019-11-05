// Copyright (c) 2019 ERGO-Code. See license.txt for license.
//
// Example for constructing crash basis and running crossover.

#include <cmath>
#include <iostream>
#include "crossover.h"

using Int = ipxint;
using Vector = ipx::Vector;

// Consider the linear constraints A*x{<,>,=}b, lb<=x<=ub given by
//
//   [1  0  2  0  0  0  0]  * x <= [ 4]
//   [0  1  2  0  0  0  0]      >= [ 3]
//   [0  0  0  2  2  2  2]      == [ 8]
//   [0  0  0 -2 -2 -2 -2]      == [-8]
//
//   [-inf]         [inf]
//   [-inf]         [inf]
//   [-inf]         [inf]
//   [   0] <= x <= [inf]
//   [-inf]         [  1]
//   [   0]         [  1]
//   [   1]         [  1]
//
// and the feasible point x = [1 1 1 1 1 1 1]. The task is to find a vertex of
// the feasible region and an associated basis.

constexpr Int num_var = 7;      // number of structural variables
constexpr Int num_constr = 4;   // number of constraints

const double obj[num_var] = { 0 };
const double ub[] = { +INFINITY, +INFINITY, +INFINITY, +INFINITY, 1, 1, 1 };
const double lb[] = { -INFINITY, -INFINITY, -INFINITY, 0, -INFINITY, 0, 1 };
const double rhs[] = { 4, 3, 8, -8 };
const char constr_type[] = { '<', '>', '=', '=' };

// Constraint matrix in CSC format with 0-based indexing.
const Int Ap[] =    { 0, 1, 2,    4,     6,     8,    10,    12 };
const Int Ai[] =    { 0, 1, 0, 1, 2,  3, 2,  3, 2,  3, 2,  3 };
const double Ax[] = { 1, 1, 2, 2, 2, -2, 2, -2, 2, -2, 2, -2 };

void PrintVbasis(const std::vector<Int>& v) {
    std::cout << "{ ";
    for (auto vi : v) {
        char c;
        switch (vi) {
        case IPX_basic:       c = 'B'; break;
        case IPX_nonbasic_lb: c = 'L'; break;
        case IPX_nonbasic_ub: c = 'U'; break;
        case IPX_superbasic:  c = 'S'; break;
        default: c = 'O';       // v does not hold a vbasis
        }
        std::cout << c << ' ';
    }
    std::cout << "}\n";
}

void PrintCbasis(const std::vector<Int>& v) {
    std::cout << "{ ";
    for (auto vi : v) {
        char c;
        switch (vi) {
        case IPX_basic:    c = 'B'; break;
        case IPX_nonbasic: c = 'N'; break;
        default: c = 'O';       // v does not hold a cbasis
        }
        std::cout << c << ' ';
    }
    std::cout << "}\n";
}

std::ostream& operator<<(std::ostream& os, const Vector& x) {
    os << "{ ";
    for (auto xi : x) os << xi << ' ';
    os << '}';
    return os;
}

std::ostream& operator<<(std::ostream& os, const ipx::Basis& basis) {
    const Int m = basis.model().rows();
    const Int n = basis.model().cols();
    os << "{ ";
    for (Int i = 0; i < n+m; i++) os << (basis.IsBasic(i) ? 'B' : 'N') << ' ';
    os << '}';
    return os;
}

int main() {
    ipx::Control control;
    ipx::Model model;
    ipx::Info info;

    ipx::Parameters parameters;
    parameters.debug = 2;       // more output
    control.parameters(parameters);

    // Form an LP model from the problem data. Because preprocessing can
    // dualize, in general model.rows() and model.cols() are different from
    // num_var and num_constr.
    model.Load(control, num_constr, num_var, Ap, Ai, Ax, rhs, constr_type, obj,
               lb, ub, &info);
    if (info.errflag) {
        std::cerr << "model error: " << info.errflag << '\n';
        return 1;
    }
    const Int m = model.rows();
    const Int n = model.cols();

    // In our case the starting point is feasible. This is not required,
    // however. The only requirement is that lb <= x <= ub and that
    //   slack[i] >= 0 for an '<' constraint,
    //   slack[i] <= 0 for an '>' constraint and
    //   slack[i] == 0 for an '=' constraint.
    // The residual A*x+slack-b can be nonzero and (in exact arithmetic) remains
    // unchanged by the crossover. The dual point must satisfy
    //   y[i] <= 0 for an '<' constraint,
    //   y[i] >= 0 for an '>' constraint,
    //   z[j] <= 0 if x[j] > lb[j],
    //   z[j] >= 0 if x[j] < ub[j]
    // and additionally we must have y[i]==0 if slack[i]!=0. (Note that the dual
    // feasibility condition for z already implies complementarity with x.) The
    // dual residual A'y+z-c can be nonzero and (in exact arithmetic) remains
    // unchanged by the crossover.
    //
    // The preprocessed model combines the structural and slack variables into
    // one vector x of size n+m with dual variables z. The constraints are all
    // equality constraints with dual variables y.
    //
    Vector x_user(num_var);
    Vector slack_user(num_constr);
    Vector y_user(num_constr);  // can also pass nullptr instead of zero vectors
    Vector z_user(num_var);     // to the routines below
    x_user = 1.0;
    slack_user[0] = 1.0;
    // slack_user[1] = -0.1;
    // y_user[0] = 0.0;
    // y_user[1] = 2.0;
    // y_user[2] = 3.0;
    Vector x(n+m), y(m), z(n+m);
    model.PresolveStartingPoint(&x_user[0], &slack_user[0], &y_user[0],
                                &z_user[0], x, y, z);

    // Construct starting basis. Each variable is given a weight in [0,Inf] that
    // specifies the priority with which the variable is put into the basis (see
    // basis.h). Instead of giving each variable weight 1.0, we would better use
    // something like the scaling factors in the IPM, hoping for fewer basis
    // updates in crossover.
    ipx::Basis basis(control, model);
    Vector weights(n+m);
    for (Int j = 0; j < n+m; j++) {
        double lbj = model.lb(j);
        double ubj = model.ub(j);
        if (lbj == ubj)
            weights[j] = 0.0;
        else if (std::isinf(lbj) && std::isinf(ubj))
            weights[j] = INFINITY;
        else
            weights[j] = 1.0;
    }
    std::cout << "Starting basis\n";
    basis.ConstructBasisFromWeights(&weights[0], &info);
    if (info.errflag) {
        std::cerr << "starting basis error: " << info.errflag << '\n';
        return 2;
    }

    ipx::Crossover crossover(control);
    std::cout << "Crossover\n";
    crossover.PushAll(&basis, x, y, z, nullptr, &info);
    if (info.errflag) {
        std::cerr << "crossover error: " << info.errflag << '\n';
        return 3;
    }
    std::cout << "basis after crossover: " << basis << '\n';
    // The basis only says "basic" or "nonbasic". The basic status of each
    // variable (required for postsolve) also depends on the vertex solution.
    std::vector<Int> basic_statuses(n+m);
    for (Int j = 0; j < n+m; j++) {
        if (basis.IsBasic(j)) {
            basic_statuses[j] = IPX_basic;
        }
        else {
            if (model.lb(j) == model.ub(j))
                basic_statuses[j] = z[j] >= 0.0 ?
                    IPX_nonbasic_lb : IPX_nonbasic_ub;
            else if (x[j] == model.lb(j))
                basic_statuses[j] = IPX_nonbasic_lb;
            else if (x[j] == model.ub(j))
                basic_statuses[j] = IPX_nonbasic_ub;
            else
                basic_statuses[j] = IPX_superbasic;
        }
    }

    // See doc/reference.pdf for how these quantities are defined.
    std::vector<Int> cbasis(num_constr), vbasis(num_var);
    model.PostsolveBasis(basic_statuses, cbasis.data(), vbasis.data());
    model.PostsolveBasicSolution(x, y, z, basic_statuses, &x_user[0],
                                 &slack_user[0], &y_user[0], &z_user[0]);

    std::cout << "x =      " << x_user << '\n';
    std::cout << "vbasis = "; PrintVbasis(vbasis);
    std::cout << "slack =  " << slack_user << '\n';
    std::cout << "cbasis = "; PrintCbasis(cbasis);

    return 0;
}
