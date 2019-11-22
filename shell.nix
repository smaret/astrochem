with import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs-channels/archive/nixos-19.09.tar.gz") {};
let
  python3Env = pkgs.python3.withPackages (ps: with ps; [ numpy h5py cython matplotlib ]);
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
    pkgs.openblas
    pkgs.sundials
  ];
}
