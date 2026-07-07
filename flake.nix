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

          cmakeFlags = [
            "-DCMAKE_BUILD_TYPE=Release"
            "-DBUILD_SHARED_LIBS=ON"
            "-DCMAKE_INSTALL_LIBDIR=lib"
          ];
        };

        pawian = pkgs.stdenv.mkDerivation {
          pname = "pawian";
          version = "unstable";
          src = pkgs.lib.cleanSource self;

          nativeBuildInputs = [
            pkgs.boost-build
            pkgs.pkg-config
          ];

          buildInputs = [
            boost170
            minuit2
            pkgs.gsl
            pkgs.root
            pkgs.log4cpp
          ];

          ROOTSYS = pkgs.root;

          BOOSTROOT = boost170;
          BOOST_ROOT = boost170;
          BOOST_INCLUDEDIR = "${boost170.dev}/include";
          BOOST_LIBRARYDIR = "${boost170.out}/lib";
          BOOSTLIBPATH = "${boost170.out}/lib";
          BOOSTINCLUDE = "${boost170.dev}/include";

          MINUIT2_STANDALONE = minuit2;

          NIX_CFLAGS_COMPILE = [
            "-I${boost170.dev}/include"
            "-I${minuit2}/include"
            "-I${pkgs.gsl}/include"
            "-I${pkgs.log4cpp}/include"
          ];

          NIX_LDFLAGS = [
            "-L${boost170.out}/lib"
            "-L${minuit2}/lib"
            "-L${pkgs.gsl}/lib"
            "-L${pkgs.log4cpp}/lib"
          ];

          enableParallelBuilding = true;
          dontConfigure = true;
      buildPhase = ''
            runHook preBuild

            ln -sfn Jamroot_andromeda Jamroot

            b2 \
              -j"$NIX_BUILD_CORES" \
              variant=release \
              include="${boost170.dev}/include" \
              library-path="${boost170.out}/lib"

            runHook postBuild
          '';

          installPhase = ''
            runHook preInstall

            mkdir -p "$out"

            b2 \
              -j"$NIX_BUILD_CORES" \
              variant=release \
              install \
              --prefix="$out" \
              include="${boost170.dev}/include" \
              library-path="${boost170.out}/lib"

            runHook postInstall
          '';
        };
      in
      {
        packages = {
          default = pawian;
          inherit pawian boost170 minuit2;
        };

        devShells.default = pkgs.mkShell {
          inputsFrom = [
            pawian
          ];

          packages = [
            pkgs.boost-build
            pkgs.pkg-config
            pkgs.cmake
            pkgs.ninja
            pkgs.gdb
          ];

          ROOTSYS = pkgs.root;

          BOOSTROOT = boost170;
          BOOST_ROOT = boost170;
          BOOST_INCLUDEDIR = "${boost170.dev}/include";
          BOOST_LIBRARYDIR = "${boost170.out}/lib";
          BOOSTLIBPATH = "${boost170.out}/lib";
          BOOSTINCLUDE = "${boost170.dev}/include";

          MINUIT2_STANDALONE = minuit2;

          shellHook = ''
            ln -sfn Jamroot_andromeda Jamroot
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

            echo "PAWIAN dev shell"
            echo "  ROOTSYS=$ROOTSYS"
            echo "  BOOSTROOT=$BOOSTROOT"
            echo "  MINUIT2_STANDALONE=$MINUIT2_STANDALONE"
            echo "  b2=$(command -v b2)"
          '';
        };
      }
    );
}
