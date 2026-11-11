#include "Highs.h"
#include "HighsExternalApi.h"
#include "catch.hpp"
#include "util/stringutil.h"

#if defined(HIGHS_SHARED_EXTRAS_LIBRARY) && defined(HIGHS_TEST_EXTRAS_DIR)
TEST_CASE("HighsExternalDeps-reload", "[highs_external_deps]") {
  HighsExternalApi::unload();
  HighsExternalApi::reset();

  // check fail (no empty search)
  REQUIRE(!HighsExternalApi::tryLoad("", false));
  const std::string failure = HighsExternalApi::getLoadStatus();

  // check fail (already initialized)
  REQUIRE(!HighsExternalApi::tryLoad(HIGHS_TEST_EXTRAS_DIR, false));
  REQUIRE(HighsExternalApi::getLoadStatus() == failure);

  // check unload
  HighsExternalApi::unload();
  REQUIRE(HighsExternalApi::instance().extras_feature_info_ == nullptr);
  REQUIRE(HighsExternalApi::getLoadStatus().empty());

  // check fail (still initialized)
  REQUIRE(!HighsExternalApi::tryLoad(HIGHS_TEST_EXTRAS_DIR, false));
  REQUIRE(HighsExternalApi::getLoadStatus().empty());

  // check success (correct path, initialization reset)
  HighsExternalApi::reset();
  REQUIRE(HighsExternalApi::tryLoad(HIGHS_TEST_EXTRAS_DIR, false));
  REQUIRE(HighsExternalApi::instance().extras_feature_info_ != nullptr);
  const std::string success = HighsExternalApi::getLoadStatus();

  // check success (no-op from previous)
  REQUIRE(HighsExternalApi::tryLoad("", false));
  REQUIRE(HighsExternalApi::getLoadStatus() == success);

  // try once more (unload and reload)
  HighsExternalApi::unload();
  REQUIRE(HighsExternalApi::instance().extras_feature_info_ == nullptr);
  REQUIRE(!HighsExternalApi::tryLoad(HIGHS_TEST_EXTRAS_DIR, false));
  REQUIRE(HighsExternalApi::instance().extras_feature_info_ == nullptr);
  HighsExternalApi::reset();
  REQUIRE(HighsExternalApi::tryLoad(HIGHS_TEST_EXTRAS_DIR, false));
  REQUIRE(HighsExternalApi::instance().extras_feature_info_ != nullptr);
}
#endif

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

// Verifies the case-insensitive substring check used to identify BLAS providers
TEST_CASE("HighsExternalDeps-isSubstringCase", "[highs_external_deps]") {
  REQUIRE(is_substring_case("Apple Accelerate", "apple"));
  REQUIRE(is_substring_case("OpenBLAS", "openblas"));
  REQUIRE(is_substring_case("libopenblas.so", "OPENBLAS"));
  REQUIRE(is_substring_case("OpenBLAS", "OpenBLAS"));
  REQUIRE(is_substring_case("OpenBLAS", ""));
  REQUIRE(!is_substring_case("Intel oneMKL", "openblas"));
  REQUIRE(!is_substring_case("BLAS", "openblas"));
  REQUIRE(!is_substring_case("", "openblas"));
  REQUIRE(!is_substring_case(nullptr, "openblas"));
  REQUIRE(!is_substring_case("OpenBLAS", nullptr));
}

// Verifies that the number of OpenBLAS threads can be set if, and only if,
// OpenBLAS is the BLAS provider. This catches builds where the BLAS metadata
// and HIPO_USES_OPENBLAS disagree
TEST_CASE("HighsExternalDeps-openblasThreads", "[highs_external_deps]") {
  if (!HighsExternalApi::isAvailable<HighsExtras::blas>()) return;

  const bool using_openblas =
      is_substring_case(HighsExtras::blas::getInfo()->provider, "openblas");
  const int initial_threads = HighsExtras::blas::openblas_set_num_threads(0);

  if (using_openblas) {
    REQUIRE(initial_threads >= 1);
    REQUIRE(HighsExtras::blas::openblas_set_num_threads(1) == 1);
    // querying does not change the number of threads
    REQUIRE(HighsExtras::blas::openblas_set_num_threads(0) == 1);
    HighsExtras::blas::openblas_set_num_threads(initial_threads);
  } else {
    REQUIRE(initial_threads == -1);
    REQUIRE(HighsExtras::blas::openblas_set_num_threads(1) == -1);
  }
}
