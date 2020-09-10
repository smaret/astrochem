{
  description = "A code to compute the abundances of chemical species in the interstellar medium";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-20.03;

  outputs = { self, nixpkgs }: {

    # TODO: Support other platforms as well
    defaultPackage.x86_64-darwin =
      with import nixpkgs { system = "x86_64-darwin"; };
      stdenv.mkDerivation {

        name = "astrochem";

        src = self;

        nativeBuildInputs = [ autoconf automake libtool ncurses ];

        buildInputs = [
          python3Packages.python
          hdf5
          sundials
        ];

        propagatedBuildInputs = [
          python3Packages.numpy
          python3Packages.h5py
          python3Packages.cython
          python3Packages.matplotlib
        ];

        preConfigure = "./bootstrap";

        doCheck = true;
      };
  };
}
