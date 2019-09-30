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

TEST_CASE("Basis-solves", "[highs_basis_solves]") {
  std::string dir = GetBasisSolvesCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  std::string filename;
  filename = dir + "/../../check/instances/adlittle.mps";
  //  filename = dir + "/check/instances/adlittle.mps";

  Highs highs;

  HighsStatus highs_status;
  highs_status = highs.initializeFromFile(filename);
  REQUIRE(highs_status == HighsStatus::OK);

  highs_status = highs.run();
  REQUIRE(highs_status == HighsStatus::OK);

  HighsLp lp = highs.getLp();
  REQUIRE(highs_status == HighsStatus::OK);

  int* basic = (int*)malloc(sizeof(int) * lp.numRow_);
  double* solution = (double*)malloc(sizeof(double) * lp.numRow_);
  double* rhs = (double*)malloc(sizeof(double) * lp.numRow_);
  
  highs_status = highs.getBasisInverseRow(0, solution);
  REQUIRE(highs_status == HighsStatus::OK);

  /*
  HighsBasis basis = getBasis();
  int basic_var_num = 0;
  for (int col=0; col < lp.numCol_; col++) {
    if (basis.col_status_[col] == HighsBasisStatus::BASIC) {
      if (basic_var_num >= lp.numRow_) break;
      printf("Basic var %d is col %d\n", basic_var_num, col);
      basic[basic_var_num++] = col;
    }
  }
  assert(basic_var_num < lp.numRow_);
  for (int row=0; row < lp.numRow_; row++) {
    if (basis.row_status_[row] == HighsBasisStatus::BASIC) {
      if (basic_var_num >= lp.numRow_) break;
      basic[basic_var_num++] = row;
    }
  }
  assert(basic_var_num == lp.numRow_);
   for (int row=0; row < lp.numRow_; row++) {
     printf(
  */
}
    
