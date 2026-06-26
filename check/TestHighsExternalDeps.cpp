#include "Highs.h"
#include "HighsExternalApi.h"
#include "catch.hpp"

// Verifies whether all features are available. If so, checks provider metadata
// for BLAS, otherwise checks for a non-empty missing-features summary
TEST_CASE("HighsExternalDeps", "[highs_external_deps]") {
  const bool all_available = HighsExternalApi::isAvailable<HighsExtras::all>();
  const std::string status = HighsExternalApi::getLoadStatus();
  REQUIRE(!status.empty());

  if (all_available) {
    REQUIRE(!std::string(HighsExtras::blas::getInfo()->provider).empty());
  } else {
    REQUIRE(!HighsExternalApi::getMissingFeatures<HighsExtras::all>().empty());
  }
}

// Verifies that the notice header points to the full notice document
TEST_CASE("HighsExternalDeps-thirdPartyNoticeHeader", "[highs_external_deps]") {
  const std::string header = HighsExternalApi::thirdPartyNoticeHeader();
  REQUIRE(header.find("THIRD_PARTY_NOTICES.md") != std::string::npos);
}

// Verifies table-style summary for all compiled features
TEST_CASE("HighsExternalDeps-getThirdPartyNoticeAll", "[highs_external_deps]") {
  const std::string notice =
      HighsExternalApi::getThirdPartyNotice<HighsExtras::all>();

  REQUIRE(notice.find("Third-party components") != std::string::npos);
  REQUIRE(notice.find("key") != std::string::npos);
  REQUIRE(notice.find("license") != std::string::npos);

  if (HighsExternalApi::isAvailable<HighsExtras::pdqsort>()) {
    REQUIRE(notice.find("pdqsort") != std::string::npos);
  } else {
    REQUIRE(notice.find("pdqsort") == std::string::npos);
  }
}

// Verifies table-style summary for the HiPO feature set
TEST_CASE("HighsExternalDeps-getThirdPartyNoticeHipo",
          "[highs_external_deps]") {
  const std::string notice =
      HighsExternalApi::getThirdPartyNotice<HighsExtras::hipo>();

  REQUIRE(notice.find("Third-party components") != std::string::npos);
  REQUIRE(notice.find("key") != std::string::npos);

  if (HighsExternalApi::isAvailable<HighsExtras::hipo>()) {
    REQUIRE(notice.find("amd") != std::string::npos);
    REQUIRE(notice.find("blas") != std::string::npos);
    REQUIRE(notice.find("metis") != std::string::npos);
    REQUIRE(notice.find("rcm") != std::string::npos);
  } else {
    REQUIRE(notice.find("amd") == std::string::npos);
    REQUIRE(notice.find("blas") == std::string::npos);
    REQUIRE(notice.find("metis") == std::string::npos);
    REQUIRE(notice.find("rcm") == std::string::npos);
  }
}

// Verifies missing-feature empty for HiPO (if available)
TEST_CASE("HighsExternalDeps-getMissingFeatures", "[highs_external_deps]") {
  const std::string missing =
      HighsExternalApi::getMissingFeatures<HighsExtras::hipo>();

  if (HighsExternalApi::isAvailable<HighsExtras::hipo>()) {
    REQUIRE(missing.empty());
  } else {
    REQUIRE(!missing.empty());
    REQUIRE((missing.find("amd") != std::string::npos ||
             missing.find("blas") != std::string::npos ||
             missing.find("metis") != std::string::npos ||
             missing.find("rcm") != std::string::npos));
  }
}