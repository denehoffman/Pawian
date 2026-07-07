{
  description = "PAWIAN build with standalone Minuit2 and Boost 1.70";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";

    boost170-src = {
      url = "https://archives.boost.io/release/1.70.0/source/boost_1_70_0.tar.bz2";
      flake = false;
    };

    minuit2-src = {
      url = "git+ssh://git@github.com/denehoffman/Minuit2.git";
      flake = false;
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      boost170-src,
      minuit2-src,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        boost170 = pkgs.boost.overrideAttrs (_oldAttrs: {
          version = "1.70.0";
          src = boost170-src;
          patches = [ ];
          postPatch = ''
            substituteInPlace boost/thread/pthread/thread_data.hpp \
              --replace-fail "#if PTHREAD_STACK_MIN > 0" "#ifdef PTHREAD_STACK_MIN" \
              --replace-fail "if (size<PTHREAD_STACK_MIN) size=PTHREAD_STACK_MIN;" \
                "if (size < static_cast<std::size_t>(PTHREAD_STACK_MIN)) size = static_cast<std::size_t>(PTHREAD_STACK_MIN);"
          '';
        });

        minuit2 = pkgs.stdenv.mkDerivation {
          pname = "minuit2-standalone";
          version = "unstable";
          src = minuit2-src;

          nativeBuildInputs = [
            pkgs.cmake
          ];
        };

        pawian = pkgs.stdenv.mkDerivation {
          pname = "pawian";
          version = "unstable";
          src = pkgs.lib.cleanSource self;

          nativeBuildInputs = [
            pkgs.cmake
            pkgs.pkg-config
          ];

          buildInputs = [
            boost170
            minuit2
            pkgs.gsl
            pkgs.root
            pkgs.log4cpp
          ];

          cmakeFlags = [
            "-DPAWIAN_USE_STANDALONE_MINUIT2=ON"
            "-DBOOST_ROOT=${boost170}"
            "-DBoost_NO_BOOST_CMAKE=ON"
            "-DBoost_NO_SYSTEM_PATHS=ON"
            "-DBoost_INCLUDE_DIR=${boost170.dev}/include"
            "-DBoost_LIBRARY_DIR_RELEASE=${boost170.out}/lib"
          ];

          ROOTSYS = pkgs.root;
        };
      in
      {
        packages = {
          default = pawian;
          inherit boost170 minuit2 pawian;
        };

        devShells.default = pkgs.mkShell {
          inputsFrom = [
            pawian
          ];

          packages = [
            pkgs.cmake
            pkgs.ninja
          ];

          ROOTSYS = pkgs.root;
          BOOST_ROOT = boost170;
          MINUIT2_STANDALONE = minuit2;
        };
      }
    );
}
