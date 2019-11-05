// Copyright (c) 2019 ERGO-Code. See license.txt for license.
//
// Example for using IPX from its C interface. The program solves the Netlib
// problem afiro.

#include <math.h>
#include <stdio.h>
#include "ipx_c.h"

typedef ipxint Int;

#define NUM_VAR 12
#define NUM_CONSTR 9
const double obj[] = { -0.2194, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.32,
                       -0.5564, 0.6, -0.48 };
const double lb[NUM_VAR] = { 0.0 };
const double ub[] = { 80.0, 283.303, 283.303, 312.813, 349.187, INFINITY,
                      INFINITY, INFINITY, 57.201, 500.0, 500.501, 357.501};
// Constraint matrix in CSC format with 0-based indexing.
const Int Ap[] = { 0, 2, 6, 10, 14, 18, 20, 22, 24, 26, 28, 30, 32 };
const Int Ai[] = { 0, 5,
                   1, 6, 7, 8,
                   2, 6, 7, 8,
                   3, 6, 7, 8,
                   4, 6, 7, 8,
                   1, 2,
                   2, 3,
                   2, 4,
                   0, 6,
                   0, 5,
                   2, 5,
                   5, 7 };
const double Ax[] = { -1.0, 0.301,
                      1.0, -1.0, 0.301, 1.06,
                      1.0, -1.0, 0.313, 1.06,
                      1.0, -1.0, 0.313, 0.96,
                      1.0, -1.0, 0.326, 0.86,
                      - 1.0, 0.99078,
                      1.00922, -1.0,
                      1.01802, -1.0,
                      1.4, 1.0,
                      0.109, -1.0,
                      -0.419111, 1.0,
                      1.4, -1.0 };
const double rhs[] = { 0.0, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 300.0 };
const char constr_type[] = { '<', '<', '=', '<', '<', '=', '<', '<', '<' };

int main() {
    void *lps = NULL;

    // Create new solver instance. This allocates a tiny amount of memory.
    ipx_new(&lps);
    if (!lps) return 1;

    // To change parameters from their defaults, initialize an ipx_parameters
    // object, make any changes and pass to the LP solver. See the reference
    // documentation for available parameters.
    struct ipx_parameters parameters = ipx_default_parameters();
    // parameters.crossover = 0;   // turns off crossover
    // parameters.debug = 1;       // sets first debugging level (more output)
    ipx_set_parameters(lps, parameters);

    // Solve the LP.
    Int status = ipx_solve(lps, NUM_VAR, obj, lb, ub, NUM_CONSTR, Ap, Ai, Ax,
                           rhs, constr_type);
    if (status != IPX_STATUS_solved) {
        // no solution
        // (invalid input, time/iter limit, numerical failure, out of memory)
        struct ipx_info info = ipx_get_info(lps);
        printf(" status: %ld, errflag: %ld\n", (long) status,
               (long) info.errflag);
        return 2;
    }

    // Get solver and solution information.
    struct ipx_info info = ipx_get_info(lps);

    // Get the interior solution (available if IPM was started).
    double x[NUM_VAR], xl[NUM_VAR], xu[NUM_VAR], slack[NUM_CONSTR];
    double y[NUM_CONSTR], zl[NUM_VAR], zu[NUM_VAR];
    ipx_get_interior_solution(lps, x, xl, xu, slack, y, zl, zu);

    // Get the basic solution (available if crossover terminated without error).
    double xbasic[NUM_VAR], sbasic[NUM_CONSTR];
    double ybasic[NUM_CONSTR], zbasic[NUM_VAR];
    Int cbasis[NUM_CONSTR], vbasis[NUM_VAR];
    ipx_get_basic_solution(lps, xbasic, sbasic, ybasic, zbasic, cbasis, vbasis);

    // Must call ipx_free() to deallocate memory in solver instance.
    ipx_free(&lps);
    return 0;
}
