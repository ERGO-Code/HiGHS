## HiGHS on nixpkgs

HiGHS now has a `flake.nix` to build the binary, allowing `nix` users to try it out

#### Python build update

Highspy with setuptools from v1.7.0 only worked on Python 3.12
For v1.7.0 we have dropped setuptools and switched to scikit-build-core
