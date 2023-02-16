#include <sstream>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;
TEST_CASE("highs-names", "[highs_names]") {
  const std::string model = "avgas";
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);
  const HighsLp& lp = highs.getLp();
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    std::stringstream ss;
    ss.str(std::string());
    ss << model << iCol << "\0";
    const std::string name = ss.str();
    printf("Col %d name is to be %s\n", int(iCol), name.c_str());
    //    REQUIRE(highs.passColName(iCol, name) == HighsStatus::kOk);
  }
  highs.run();
  REQUIRE(highs.writeModel("") == HighsStatus::kOk);
  highs.writeSolution("", 1);
}
