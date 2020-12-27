with import
  (fetchTarball {
    # Get the revision by choosing a version from https://github.com/nixos/nixpkgs/commits/master
    url = "https://github.com/nixos/nixpkgs/archive/171144eb2850979d3841949df5a14f048da5d83e.tar.gz";
    # Get the hash by running `nix-prefetch-url --unpack <url>` on the above url
    sha256 = "0hm78p3avq9khd97x3x1aziribsz0mnrrssanq982yn5kz2jfj6s";
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
