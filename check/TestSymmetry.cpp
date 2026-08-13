#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLp.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "parallel/HighsTaskExecutor.h"
#include "presolve/HighsSymmetry.h"

const bool dev_run = false;

TEST_CASE("symmetry-orbitope-detection", "[highs_test_symmetry]") {
  // 3 tasks assigned to 4 identical machines.
  // Variables x[i][j] = 1 if task i uses machine j.
  // Constraints:
  //   x[i][1] + x[i][2] + x[i][3] + x[i][4] >= 1  (coverage, i=1,2,3)
  //   x[1][j] + x[2][j] + x[3][j] <= 1            (packing , j=1,2,3,4)
  //
  // The 4 machines are interchangeable => full orbitope with 3 rows
  // and 4 columns.
  HighsLp lp;
  lp.num_col_ = 12;
  lp.num_row_ = 7;
  // Different costs per task to break task-symmetry, keeping only machine
  // symmetry
  lp.col_cost_ = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  lp.col_lower_.assign(12, 0.0);
  lp.col_upper_.assign(12, 1.0);
  lp.integrality_.assign(12, HighsVarType::kInteger);
  lp.row_lower_ = {1, 1, 1, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_.assign(7, 1.0);
  lp.a_matrix_.start_ = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
  lp.a_matrix_.index_ = {0, 3, 0, 4, 0, 5, 0, 6, 1, 3, 1, 4,
                         1, 5, 1, 6, 2, 3, 2, 4, 2, 5, 2, 6};
  lp.a_matrix_.value_.assign(24, 1.0);

  HighsTaskExecutor::initialize(1);

  HighsSymmetryDetection detection;
  detection.loadModelAsGraph(lp, 1e-8);

  bool initialized = detection.initializeDetection();
  REQUIRE(initialized);

  // run symmetry detection
  HighsSymmetries symmetries;
  detection.run(symmetries);

  if (dev_run) {
    printf("numGenerators = %d\n", static_cast<int>(symmetries.numGenerators));
    printf("numPerms = %d\n", static_cast<int>(symmetries.numPerms));
    printf("orbitopes = %d\n", static_cast<int>(symmetries.orbitopes.size()));
  }

  // Should find generators (4 machines need 3 generators)
  REQUIRE(symmetries.numGenerators > 0);

  // Should find exactly one full orbitope (3 rows x 4 columns)
  REQUIRE(symmetries.orbitopes.size() == 1);
  REQUIRE(symmetries.orbitopes[0].numRows == 3);
  REQUIRE(symmetries.orbitopes[0].rowLength == 4);

  // Exercise determineOrbitopeType / detectSetPackingRows.
  // Register assignment cliques: for each task i, x[i][1..4] form a clique.
  HighsCliqueTable cliquetable(12);
  HighsCliqueTable::CliqueVar clique0[] = {{0, 1}, {1, 1}, {2, 1}, {3, 1}};
  HighsCliqueTable::CliqueVar clique1[] = {{4, 1}, {5, 1}, {6, 1}, {7, 1}};
  HighsCliqueTable::CliqueVar clique2[] = {{8, 1}, {9, 1}, {10, 1}, {11, 1}};
  cliquetable.doAddClique(clique0, 4, true);
  cliquetable.doAddClique(clique1, 4, true);
  cliquetable.doAddClique(clique2, 4, true);

  symmetries.orbitopes[0].determineOrbitopeType(cliquetable);

  // All 3 rows should be detected as set packing
  REQUIRE(symmetries.orbitopes[0].numSetPackingRows == 3);

  HighsTaskExecutor::shutdown();
}

TEST_CASE("symmetry-orbital-fixing", "[highs_test_symmetry]") {
  // Unit test for orbitalFixing (which dispatches to
  // orbitalFixingForPackingOrbitope when all rows are set-packing).
  // Use a Highs object to set up the MIP solver infrastructure, then
  // call orbital fixing directly on the detected orbitope.
  HighsLp lp;
  lp.num_col_ = 12;
  lp.num_row_ = 7;
  lp.col_cost_ = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  lp.col_lower_.assign(12, 0.0);
  lp.col_upper_.assign(12, 1.0);
  lp.integrality_.assign(12, HighsVarType::kInteger);
  lp.row_lower_ = {1, 1, 1, -kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_.assign(7, 1.0);
  lp.a_matrix_.start_ = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
  lp.a_matrix_.index_ = {0, 3, 0, 4, 0, 5, 0, 6, 1, 3, 1, 4,
                         1, 5, 1, 6, 2, 3, 2, 4, 2, 5, 2, 6};
  lp.a_matrix_.value_.assign(24, 1.0);

  // Use Highs to get MIP solver infrastructure
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->feastol = 1e-6;
  mipsolver.mipdata_->setupDomainPropagation();

  HighsTaskExecutor::initialize(1);

  // Detect symmetry
  HighsSymmetryDetection detection;
  detection.loadModelAsGraph(lp, 1e-8);
  REQUIRE(detection.initializeDetection());

  HighsSymmetries symmetries;
  detection.run(symmetries);
  REQUIRE(symmetries.orbitopes.size() == 1);

  // Set up orbitope as packing
  HighsCliqueTable::CliqueVar clique0[] = {{0, 1}, {1, 1}, {2, 1}, {3, 1}};
  HighsCliqueTable::CliqueVar clique1[] = {{4, 1}, {5, 1}, {6, 1}, {7, 1}};
  HighsCliqueTable::CliqueVar clique2[] = {{8, 1}, {9, 1}, {10, 1}, {11, 1}};
  mipsolver.mipdata_->cliquetable.doAddClique(clique0, 4, true);
  mipsolver.mipdata_->cliquetable.doAddClique(clique1, 4, true);
  mipsolver.mipdata_->cliquetable.doAddClique(clique2, 4, true);
  symmetries.orbitopes[0].determineOrbitopeType(
      mipsolver.mipdata_->cliquetable);
  REQUIRE(symmetries.orbitopes[0].numSetPackingRows == 3);

  const auto& orb = symmetries.orbitopes[0];

  // Branch on columns in two different rows of the orbitope.
  // Row 0 at column position 0, row 1 at column position 1.
  // This creates a scenario where orbital fixing can deduce cross-row fixings.
  HighsDomain& domain = mipsolver.mipdata_->getDomain();
  domain.changeBound(HighsBoundType::kLower, orb(0, 0), 1.0,
                     HighsDomain::Reason::branching());
  domain.changeBound(HighsBoundType::kLower, orb(1, 1), 1.0,
                     HighsDomain::Reason::branching());

  // Call the public orbitalFixing method
  HighsInt numFixed = orb.orbitalFixing(domain);

  if (dev_run) {
    printf("numFixed = %d\n", static_cast<int>(numFixed));
  }

  // Orbital fixing should deduce additional fixings from the orbitope structure
  REQUIRE(numFixed > 0);

  HighsTaskExecutor::shutdown();
}
