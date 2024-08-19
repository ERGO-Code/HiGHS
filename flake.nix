{
  inputs = {
    nixpkgs = {
      url = "github:nixos/nixpkgs/nixos-unstable";
    };
    flake-utils = {
      url = "github:numtide/flake-utils";
    };
  };
  outputs = { nixpkgs, flake-utils, ... }: flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = import nixpkgs {
        inherit system;
      };

      version = with pkgs.lib;
        # Read the version. Note: We assume the version numbers are in
        # order in the file; i.e. Major, Minor, Patch.
        let f = builtins.readFile ./Version.txt;
        l = strings.splitString "\n" f;
        # Drop the last term; it just says if it's alpha or not.
        t = lists.take 3 l;
        # Get the numbers on the other side of the equals
        vs = lists.forEach t (v: lists.drop 1 (strings.splitString "=" v));
        # That's it!
        in concatStrings (intersperse "." (lists.flatten vs));

      # Build the binary
      highs-binary = with pkgs; stdenv.mkDerivation {
        pname = "highs";
        inherit version;
        src = pkgs.lib.cleanSource ./.;
        nativeBuildInputs = [
          clang
          cmake
        ];
      };

      # Build the python package
      highspy = pkgs.python3Packages.buildPythonPackage {
        inherit version;
        pname = "highspy";
        src = pkgs.lib.cleanSource ./.;
        format = "pyproject";
        dontUseCmakeConfigure = true;
        nativeBuildInputs = with pkgs.python3Packages; [
          numpy
          pathspec
          pybind11
          pyproject-metadata
          scikit-build-core
          pkgs.cmake
          pkgs.ninja
        ];
        buildInputs = [
          pkgs.zlib
        ];
      };

    in rec {
      defaultApp = flake-utils.lib.mkApp {
        drv = highs-binary;
      };
      defaultPackage = highs-binary;
      packages.highspy = highspy;
      devShells.highspy = pkgs.mkShell {
        buildInputs = [
          (pkgs.python3.withPackages (ps: [ highspy ]) )
        ];
      };
    }
  );
}
