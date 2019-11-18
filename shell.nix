with import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs-channels/archive/nixos-19.03.tar.gz") {};
let
  python3Env = pkgs.python3.withPackages (ps: with ps; [ numpy h5py cython matplotlib ]);
  sundials2 = pkgs.sundials.overrideAttrs (oldAttrs: rec {
    name = "sundials-2.7.0";
    src = pkgs.fetchurl {
      url = "https://computation.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz";
      sha256 = "01513g0j7nr3rh7hqjld6mw0mcx5j9z9y87bwjc16w2x2z3wm7yk";
    };
  });
in
pkgs.mkShell rec {
  buildInputs = [
    python3Env
    pkgs.gfortran
    pkgs.hdf5
    pkgs.autoconf
    pkgs.automake
    pkgs.libtool
    pkgs.ncurses
    sundials2
    pkgs.openblas
  ];
}
