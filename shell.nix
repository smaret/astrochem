with import
  (fetchTarball {
    # Get the revision by choosing a version from https://github.com/nixos/nixpkgs/commits/master
    url = "https://github.com/nixos/nixpkgs/archive/22.05.tar.gz";
    # Get the hash by running `nix-prefetch-url --unpack <url>` on the above url
    sha256 = "0d643wp3l77hv2pmg2fi7vyxn4rwy0iyr8djcw1h5x72315ck9ik";
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
