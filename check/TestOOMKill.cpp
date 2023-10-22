#include "HCheckConfig.h"
#include "Highs.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "catch.hpp"
#include "lp_data/HConst.h"

TEST_CASE("linprog_oom", "[highs_solver]") {
    const HighsInt n_ctr = 500000;
    const HighsInt n_var = 500;

    // Seed random number generator
    std::srand(0);

    // Create c_ vector of ones
    std::vector<double> c_(n_var, 1.0);

    // Create A_ub matrix with random values between 0 and 1
    std::vector<double> Avalue(n_ctr * n_var);
    for (int i = 0; i < n_ctr * n_var; ++i) {
        Avalue[i] = (double)rand() / RAND_MAX;
    }

    // Create b_ub vector of zeros and set bounds
    std::vector<double> rowUpper(n_ctr, 0.0);
    std::vector<double> rowLower(n_ctr, -kHighsInf);

    // Assuming A_ub is in column-wise format, similar to the provided example
    std::vector<HighsInt> Astart(n_var + 1);
    std::vector<HighsInt> Aindex(n_ctr * n_var);
    for (int col = 0; col < n_var; ++col) {
        Astart[col] = col * n_ctr;
    }
    for (int idx = 0; idx < n_ctr * n_var; ++idx) {
        Aindex[idx] = idx % n_ctr;
    }
    Astart[n_var] = n_ctr * n_var;

    // Column bounds are (0, infinity)
    std::vector<double> colUpper(n_var, kHighsInf);
    std::vector<double> colLower(n_var, 0);

    Highs highs;

    // Set up the LP externally
    HighsLp lp;
    lp.num_col_ = n_var;
    lp.num_row_ = n_ctr;
    lp.col_cost_ = c_;
    lp.col_lower_ = colLower;
    lp.col_upper_ = colUpper;
    lp.row_lower_ = rowLower;
    lp.row_upper_ = rowUpper;
    lp.a_matrix_.start_ = Astart;
    lp.a_matrix_.index_ = Aindex;
    lp.a_matrix_.value_ = Avalue;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    highs.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
    // highs.setOptionValue("threads", 1);
    REQUIRE(highs.setOptionValue("presolve", "on") == HighsStatus::kOk);
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
}


TEST_CASE("linprog_oom_no_presolve", "[highs_solver]") {
    const HighsInt n_ctr = 500000;
    const HighsInt n_var = 500;

    // Seed random number generator
    std::srand(0);

    // Create c_ vector of ones
    std::vector<double> c_(n_var, 1.0);

    // Create A_ub matrix with random values between 0 and 1
    std::vector<double> Avalue(n_ctr * n_var);
    for (int i = 0; i < n_ctr * n_var; ++i) {
        Avalue[i] = (double)rand() / RAND_MAX;
    }

    // Create b_ub vector of zeros and set bounds
    std::vector<double> rowUpper(n_ctr, 0.0);
    std::vector<double> rowLower(n_ctr, -kHighsInf);

    // Assuming A_ub is in column-wise format, similar to the provided example
    std::vector<HighsInt> Astart(n_var + 1);
    std::vector<HighsInt> Aindex(n_ctr * n_var);
    for (int col = 0; col < n_var; ++col) {
        Astart[col] = col * n_ctr;
    }
    for (int idx = 0; idx < n_ctr * n_var; ++idx) {
        Aindex[idx] = idx % n_ctr;
    }
    Astart[n_var] = n_ctr * n_var;

    // Column bounds are (0, infinity)
    std::vector<double> colUpper(n_var, kHighsInf);
    std::vector<double> colLower(n_var, 0);

    Highs highs;

    // Set up the LP externally
    HighsLp lp;
    lp.num_col_ = n_var;
    lp.num_row_ = n_ctr;
    lp.col_cost_ = c_;
    lp.col_lower_ = colLower;
    lp.col_upper_ = colUpper;
    lp.row_lower_ = rowLower;
    lp.row_upper_ = rowUpper;
    lp.a_matrix_.start_ = Astart;
    lp.a_matrix_.index_ = Aindex;
    lp.a_matrix_.value_ = Avalue;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    highs.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
    // highs.setOptionValue("threads", 1);
    REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
}
