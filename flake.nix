{
  description = "A code to compute the abundances of chemical species in the interstellar medium";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-20.03;
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system}; in
      rec {
        packages = flake-utils.lib.flattenTree {
          astrochem = pkgs.stdenv.mkDerivation {

            name = "astrochem";

            src = self;

            nativeBuildInputs = [ pkgs.autoconf pkgs.automake pkgs.libtool pkgs.ncurses ];

            buildInputs = [
              pkgs.python3Packages.python
              pkgs.hdf5
              pkgs.sundials
            ];

            propagatedBuildInputs = [
              pkgs.python3Packages.numpy
              pkgs.python3Packages.h5py
              pkgs.python3Packages.cython
              pkgs.python3Packages.matplotlib
            ];

            preConfigure = "./bootstrap";

            doCheck = true;
          };
        };
        defaultPackage = packages.astrochem;
      }                     
    );
}
