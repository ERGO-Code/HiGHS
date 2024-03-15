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
            in concatStrings (lib.intersperse "." (lists.flatten vs));
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

