{
  description = "Spatial tree library for Rust";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
    in {
      packages.x86_64-linux.shell = pkgs.stdenv.mkDerivation {
        name = "acacia";
        src = self;

        buildInputs = with pkgs; [
          cargo
          entr
          fd
          rustc
          rustfmt
        ];

        installPhase = "touch $out";
      };
      defaultPackage.x86_64-linux = self.packages.x86_64-linux.shell;
    };
}
