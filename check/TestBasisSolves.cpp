#include "Highs.h"

#include "catch.hpp"

#ifdef __linux__
#include <unistd.h>
#elif _WIN32
#define NOGDI
#include <windows.h>
#else

#endif

std::string GetBasisSolvesCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];

  #ifdef __linux__ 
    auto result = getcwd(buff, FILENAME_MAX);
    if (result) {
    std::string current_working_dir(buff);
    return current_working_dir;
  }
  #elif _WIN32
    GetModuleFileName( NULL, buff, FILENAME_MAX );
    string::size_type pos = string( buff ).find_last_of( "\\/" );
    return string( buff ).substr( 0, pos);
  #else

  #endif

  return "";
}

double GetBasisSolvesCheckSolution(HighsLp& lp, int* basic_variables, double* rhs, double* solution, const bool transpose=false) {
  double residual_norm = 0;
  if (transpose) {
    for (int k=0; k<lp.numRow_; k++) {
      int var = basic_variables[k];
      double residual;
      if (var < 0) {
	int row = -(1+var);
	residual = fabs(rhs[k] - solution[row]);
	if (residual > 1e-8) printf("Row |[B^Tx-b]_{%2d}| = %11.4g\n", k,  residual);
      } else {
	int col = var;
	residual = 0;
	for (int el=lp.Astart_[col]; el<lp.Astart_[col+1]; el++) {
	  int row = lp.Aindex_[col];
	  residual += lp.Avalue_[col]*solution[row];
	}
	residual = fabs(rhs[k] - residual);
	if (residual > 1e-8) printf("Col |[B^Tx-b]_{%2d}| = %11.4g\n", k,  residual);
      }
      residual_norm += residual;
    }
  } else {
    vector<double> basis_matrix_times_solution;
    basis_matrix_times_solution.assign(lp.numRow_, 0);
    for (int k=0; k<lp.numRow_; k++) {
      int var = basic_variables[k];
      if (var < 0) {
	int row = -(1+var);
	basis_matrix_times_solution[row] = solution[k];
      } else {
	int col = var;
	for (int el=lp.Astart_[col]; el<lp.Astart_[col+1]; el++) {
	  int row = lp.Aindex_[col];
	  basis_matrix_times_solution[row] += lp.Avalue_[col]*solution[k];
	}
      }
    }
    for (int k=0; k<lp.numRow_; k++) {
      double residual = fabs(rhs[k] - basis_matrix_times_solution[k]);
      if (residual > 1e-8) printf("|[B^Tx-b]_{%2d}| = %11.4g\n", k,  residual); 
      residual_norm += residual;
    }
  }
  return residual_norm;
}

// No commas in test case name.
TEST_CASE("Basis-solves", "[highs_basis_solves]") {

  std::string dir = GetBasisSolvesCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  std::string filename;
  filename = dir + "/../../check/instances/avgas.mps";
  //  filename = dir + "/../../check/instances/adlittle.mps";
  //  filename = dir + "/check/instances/adlittle.mps";

  Highs highs;

  int* basic_variables;
  double* rhs;
  double* solution;

  HighsStatus highs_status;

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisInverseRow(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisInverseCol(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisTransposeSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getReducedColumn(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.initializeFromFile(filename);
  REQUIRE(highs_status==HighsStatus::OK);
 
  HighsLp lp = highs.getLp();
  REQUIRE(highs_status == HighsStatus::OK);

  int numRow = lp.numRow_;
  basic_variables = (int*)malloc(sizeof(int) * numRow);
  solution = (double*)malloc(sizeof(double) * numRow);
  rhs = (double*)malloc(sizeof(double) * numRow);

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisInverseRow(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisInverseCol(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getBasisTransposeSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.getReducedColumn(0, solution);
  REQUIRE(highs_status==HighsStatus::Error);

  highs_status = highs.run();
  REQUIRE(highs_status == HighsStatus::OK);

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status==HighsStatus::OK);

  highs_status = highs.getBasisInverseRow(0, solution);
  REQUIRE(highs_status==HighsStatus::OK);

  highs_status = highs.getBasisInverseCol(0, solution);
  REQUIRE(highs_status==HighsStatus::OK);

  highs_status = highs.getBasisSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::OK);

  highs_status = highs.getBasisTransposeSolve(rhs, solution);
  REQUIRE(highs_status==HighsStatus::OK);

  highs_status = highs.getReducedColumn(0, solution);
  REQUIRE(highs_status==HighsStatus::OK);

  for (int row=0; row < numRow; row++) {
    printf("Basic variable %3d is ", row);
    int var = basic_variables[row];
    if (var<0) {
      printf("row %d\n", -(1+var));
    } else {
      printf("col %d\n", var);
    }
  }

}
    
