{
  description = "PAWIAN build with standalone Minuit2 and Boost 1.70";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    nixpkgs-unstable.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";

    boost170-src = {
      url = "https://archives.boost.io/release/1.70.0/source/boost_1_70_0.tar.bz2";
      flake = false;
    };

    minuit2-src = {
      url = "git+ssh://git@github.com/denehoffman/Minuit2.git";
      flake = false;
    };

    amptools-src = {
      url = "github:mashephe/AmpTools/v0.15.4";
      flake = false;
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      nixpkgs-unstable,
      flake-utils,
      boost170-src,
      minuit2-src,
      amptools-src,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        pkgsUnstable = import nixpkgs-unstable {
          inherit system;
        };

        python = pkgsUnstable.python314;

        root = pkgs.root.override {
          inherit python;
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

        amptools = pkgs.stdenv.mkDerivation {
          pname = "amptools";
          version = "0.15.4";
          src = amptools-src;

          nativeBuildInputs = [
            pkgs.gnumake
            pkgs.pkg-config
          ];

          buildInputs = [
            root
          ];

          ROOTSYS = root;

          buildPhase = ''
            runHook preBuild

            export PATH="${root}/bin:$PATH"
            export AMPTOOLS_HOME="$PWD"

            make -C AmpTools -j"$NIX_BUILD_CORES"

            runHook postBuild
          '';

          installPhase = ''
            runHook preInstall

            mkdir -p "$out/include" "$out/lib" "$out/share/AmpTools"

            # Install usable include tree.
            cp -r AmpTools/IUAmpTools "$out/include/" 2>/dev/null || true
            cp -r AmpTools/DataIO "$out/include/" 2>/dev/null || true
            cp -r AmpTools/Amplitude "$out/include/" 2>/dev/null || true
            cp -r AmpTools/Utilities "$out/include/" 2>/dev/null || true

            # Install built libraries.
            cp -r AmpTools/lib/* "$out/lib/" 2>/dev/null || true

            # Keep full source tree for AmpTools' own Makefile-style workflows.
            cp -r AmpTools Makefile.settings "$out/share/AmpTools/"

            runHook postInstall
          '';
        };

        pawian = pkgs.stdenv.mkDerivation {
          pname = "pawian";
          version = "unstable";
          src = pkgs.lib.cleanSource self;

          nativeBuildInputs = [
            pkgs.boost-build
            pkgs.pkg-config
          ] ++ root.nativeBuildInputs;
          propagatedBuildInputs = root.propagatedBuildInputs;

          buildInputs = [
            boost170
            minuit2
            pkgs.gsl
            root
            pkgs.log4cpp
          ] ++ root.buildInputs;

          ROOTSYS = root;

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
          inherit pawian boost170 minuit2 amptools;
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
            python
            pkgs.uv
          ];

          ROOTSYS = root;

          BOOSTROOT = boost170;
          BOOST_ROOT = boost170;
          BOOST_INCLUDEDIR = "${boost170.dev}/include";
          BOOST_LIBRARYDIR = "${boost170.out}/lib";
          BOOSTLIBPATH = "${boost170.out}/lib";
          BOOSTINCLUDE = "${boost170.dev}/include";

          MINUIT2_STANDALONE = minuit2;

          AMPTOOLS_HOME = "${amptools}/share/AmpTools";
          AMPTOOLS = "${amptools}/share/AmpTools/AmpTools";

          shellHook = ''
            export TOP_DIR="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"

            ln -sfn "$TOP_DIR/Jamroot_andromeda" "$TOP_DIR/Jamroot"

            export ROOTSYS="${root}"
            export BOOST_BUILD_PATH="${pkgs.boost-build}/share/boost-build"

            export AMPTOOLS_HOME="${amptools}/share/AmpTools"
            export AMPTOOLS="$AMPTOOLS_HOME/AmpTools"

            export PATH="$ROOTSYS/bin:$TOP_DIR/bin:$PATH"

            export LD_LIBRARY_PATH="$TOP_DIR/lib:${pkgs.lib.makeLibraryPath [
              root
              boost170
              minuit2
              amptools
              pkgs.gsl
              pkgs.log4cpp
            ]}''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

            export PYTHONPATH="${root}/${python.sitePackages}:''${PYTHONPATH:-}"
            export NIX_PYTHON="${python}/bin/python"
            export UV_NO_MANAGED_PYTHON=1

            export EVT_DIR="''${EVT_DIR:-$TOP_DIR/Event}"
            export KMAT_DIR="''${KMAT_DIR:-$TOP_DIR/KMatStore}"
            mkdir -p "$EVT_DIR" "$KMAT_DIR"

            echo "PAWIAN dev shell"
            echo "  TOP_DIR=$TOP_DIR"
            echo "  ROOTSYS=$ROOTSYS"
            echo "  EVT_DIR=$EVT_DIR"
            echo "  KMAT_DIR=$KMAT_DIR"
            echo "  BOOST_BUILD_PATH=$BOOST_BUILD_PATH"
          '';
        };
      }
    );
}
