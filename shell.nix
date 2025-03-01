with import
  (fetchTarball {
    url= "https://github.com/nixos/nixpkgs/archive/refs/heads/nixos-24.11.tar.gz";
  })
{ };
let
  python3Env = pkgs.python3.withPackages (ps: with ps; [ numpy h5py cython matplotlib ]);
in
pkgs.mkShell rec {
  buildInputs = [
    python3Env
    pkgs.hdf5
    pkgs.autoconf
    pkgs.automake
    pkgs.libtool
    pkgs.ncurses
    pkgs.sundials
  ];
}
