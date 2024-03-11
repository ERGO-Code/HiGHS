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
      highs = (with pkgs; stdenv.mkDerivation {
          pname = "highs";
          version = "1.6.0";
          src = pkgs.lib.cleanSource ./.;
          nativeBuildInputs = [
            clang
            cmake
          ];
        }
      );
    in rec {
      defaultApp = flake-utils.lib.mkApp {
        drv = defaultPackage;
      };
      defaultPackage = highs;
      devShell = pkgs.mkShell {
        buildInputs = [
          highs
        ];
      };
    }
  );
}

