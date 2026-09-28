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

static HighsLp buildAssignmentModel() {
  // 3 tasks assigned to 4 identical machines.
  // Variables x[i][j] = 1 if task i uses machine j.
  // Constraints:
  //   x[i][1] + x[i][2] + x[i][3] + x[i][4] >= 1  (coverage, i=1,2,3)
  //   x[1][j] + x[2][j] + x[3][j] <= 1            (packing , j=1,2,3,4)
  //
  // Different costs per task break task-symmetry, keeping only machine
  // symmetry. The 4 machines are interchangeable => full orbitope with 3 rows
  // and 4 columns.
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
  return lp;
}

TEST_CASE("symmetry-orbitope-detection", "[highs_test_symmetry]") {
  HighsLp lp = buildAssignmentModel();

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

TEST_CASE("symmetry-get-branching-column", "[highs_test_symmetry]") {
  // Exercise getBranchingColumn: if we ask to branch on a column that is not
  // the leftmost unfixed in its orbitope row, it should return the leftmost.
  HighsLp lp = buildAssignmentModel();

  HighsTaskExecutor::initialize(1);

  HighsSymmetryDetection detection;
  detection.loadModelAsGraph(lp, 1e-8);
  REQUIRE(detection.initializeDetection());

  HighsSymmetries symmetries;
  detection.run(symmetries);
  REQUIRE(symmetries.orbitopes.size() == 1);

  HighsCliqueTable cliquetable(12);
  HighsCliqueTable::CliqueVar clique0[] = {{0, 1}, {1, 1}, {2, 1}, {3, 1}};
  HighsCliqueTable::CliqueVar clique1[] = {{4, 1}, {5, 1}, {6, 1}, {7, 1}};
  HighsCliqueTable::CliqueVar clique2[] = {{8, 1}, {9, 1}, {10, 1}, {11, 1}};
  cliquetable.doAddClique(clique0, 4, true);
  cliquetable.doAddClique(clique1, 4, true);
  cliquetable.doAddClique(clique2, 4, true);
  symmetries.orbitopes[0].determineOrbitopeType(cliquetable);

  const auto& orb = symmetries.orbitopes[0];
  std::vector<double> colLower(12, 0.0);
  std::vector<double> colUpper(12, 1.0);

  // With all columns unfixed, branching on orb(0,1) should redirect to orb(0,0)
  HighsInt redirected = orb.getBranchingColumn(colLower, colUpper, orb(0, 1));
  REQUIRE(redirected == orb(0, 0));

  // If orb(0,0) is already fixed, branching on orb(0,1) should stay at orb(0,1)
  colLower[orb(0, 0)] = 1.0;
  HighsInt notRedirected =
      orb.getBranchingColumn(colLower, colUpper, orb(0, 1));
  REQUIRE(notRedirected == orb(0, 1));

  HighsTaskExecutor::shutdown();
}

TEST_CASE("symmetry-orbital-fixing", "[highs_test_symmetry]") {
  // Unit test for orbitalFixing (which dispatches to
  // orbitalFixingForPackingOrbitope when all rows are set-packing).
  HighsLp lp = buildAssignmentModel();

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

TEST_CASE("symmetry-orbital-fixing-full-orbitope", "[highs_test_symmetry]") {
  // Unit test for orbitalFixingForFullOrbitope.
  // Without registering cliques, all rows are kRowNotPacking and the
  // full orbitope path is taken.
  HighsLp lp = buildAssignmentModel();

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

  // Call determineOrbitopeType with an empty clique table so columnToRow
  // is populated but no rows are detected as set-packing.
  HighsCliqueTable cliquetable(12);
  symmetries.orbitopes[0].determineOrbitopeType(cliquetable);
  REQUIRE(symmetries.orbitopes[0].numSetPackingRows == 0);

  const auto& orb = symmetries.orbitopes[0];

  // Branch by fixing entries to zero — this leaves free entries for the
  // full orbitope algorithm to resolve via lexicographic reasoning.
  HighsDomain& domain = mipsolver.mipdata_->getDomain();
  domain.changeBound(HighsBoundType::kUpper, orb(0, 0), 0.0,
                     HighsDomain::Reason::branching());
  domain.changeBound(HighsBoundType::kUpper, orb(1, 0), 0.0,
                     HighsDomain::Reason::branching());

  // Call orbitalFixing — should dispatch to orbitalFixingForFullOrbitope
  HighsInt numFixed = orb.orbitalFixing(domain);

  if (dev_run) {
    printf("numFixed (full orbitope) = %d\n", static_cast<int>(numFixed));
  }

  // Full orbitope fixing should deduce fixings from the lexicographic structure
  REQUIRE(numFixed > 0);

  HighsTaskExecutor::shutdown();
}

TEST_CASE("symmetry-propagate-orbitopes", "[highs_test_symmetry]") {
  // Unit test for HighsSymmetries::propagateOrbitopes which looks up
  // branched columns in columnToOrbitope and dispatches to orbitalFixing.
  HighsLp lp = buildAssignmentModel();

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

  HighsSymmetryDetection detection;
  detection.loadModelAsGraph(lp, 1e-8);
  REQUIRE(detection.initializeDetection());

  HighsSymmetries symmetries;
  detection.run(symmetries);
  REQUIRE(symmetries.orbitopes.size() == 1);

  // Set up as packing orbitope
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

  // Branch on two rows
  HighsDomain& domain = mipsolver.mipdata_->getDomain();
  domain.changeBound(HighsBoundType::kLower, orb(0, 0), 1.0,
                     HighsDomain::Reason::branching());
  domain.changeBound(HighsBoundType::kLower, orb(1, 1), 1.0,
                     HighsDomain::Reason::branching());

  // Call propagateOrbitopes on HighsSymmetries
  HighsInt numFixed = symmetries.propagateOrbitopes(domain);

  if (dev_run) {
    printf("numFixed (propagateOrbitopes) = %d\n", static_cast<int>(numFixed));
  }

  REQUIRE(numFixed > 0);

  HighsTaskExecutor::shutdown();
}

TEST_CASE("issue-3166", "[highs_test_symmetry]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/issue-3166.mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);

  const HighsLp& lp = highs.getLp();

  HighsTaskExecutor::initialize(1);

  HighsSymmetryDetection detection;
  detection.loadModelAsGraph(lp, 1e-8);

  bool initialized = detection.initializeDetection();
  REQUIRE(initialized);

  HighsSymmetries symmetries;
  detection.run(symmetries);

  REQUIRE(symmetries.numGenerators >= 0);

  HighsTaskExecutor::shutdown();
}
