# nixpkgs (nixos-unstable) with `highs` and `highspy` built from this checkout.
# Used by .github/workflows/nix-revdeps.yml; evaluate with --impure.
let
  nixpkgs = builtins.getFlake "github:NixOS/nixpkgs/nixos-unstable";
  lib = nixpkgs.lib;

  src = lib.cleanSource ../..;

  # Version.txt starts with HIGHS_MAJOR=, HIGHS_MINOR=, HIGHS_PATCH= lines.
  # The derivation version must match it, otherwise nixpkgs' versionCheckPhase
  # (highs) and the dist-info METADATA check (highspy) fail.
  version = lib.concatStringsSep "." (map
    (line: lib.last (lib.splitString "=" line))
    (lib.take 3 (lib.splitString "\n" (builtins.readFile ../../Version.txt))));

  overlay = final: prev: {
    highs = prev.highs.overrideAttrs (old: { inherit version src; });
    pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
      (pyfinal: pyprev: {
        highspy = pyprev.highspy.overridePythonAttrs (old: { inherit version src; });
      })
    ];
  };
in
nixpkgs.legacyPackages.${builtins.currentSystem}.extend overlay
