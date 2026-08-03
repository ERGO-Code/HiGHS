import re
import unittest

import highspy_extras


class TestHighsPyExtras(unittest.TestCase):
    def test_version(self):
        # Ensure the library version matches the release part of the package version,
        # e.g., "1.15.1.dev1" -> "1.15.1"
        package_release = re.match(r"\d+(?:\.\d+)*", highspy_extras.__version__)
        if package_release is None:
            self.fail(f"Could not parse a release version from {highspy_extras.__version__!r}")
        self.assertEqual(highspy_extras.library.version, package_release.group(0))

    def test_features(self):
        # Ensure the library has non-zero features available
        features = highspy_extras.library.features
        self.assertIsInstance(features, dict)
        self.assertGreater(len(features), 0)

    def test_blas_enabled(self):
        # Ensure BLAS feature is present and enabled in this build
        self.assertIn("blas", highspy_extras.library.features)
        blas_feature = highspy_extras.library["blas"]
        self.assertTrue(blas_feature.enabled)

    def test_getitem_matches_features_mapping(self):
        # Ensure __getitem__ and features[] are the same
        self.assertIs(
            highspy_extras.library["blas"],
            highspy_extras.library.features["blas"],
        )

    def test_feature_table(self):
        # Basic sanity check for the feature table
        table = highspy_extras.library.feature_table
        self.assertIsInstance(table, str)
        self.assertIn("key", table)
        self.assertIn("version", table)
        self.assertIn("blas", table)
