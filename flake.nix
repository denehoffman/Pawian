{
  description = "PAWIAN and AmpTools comparison development environment";

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
    flake-utils.lib.eachSystem
      [
        "aarch64-darwin"
        "x86_64-linux"
        "aarch64-linux"
      ]
      (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          pkgsUnstable = import nixpkgs-unstable { inherit system; };
          lib = pkgs.lib;
          isDarwin = pkgs.stdenv.isDarwin;
          isLinux = pkgs.stdenv.isLinux;
          sharedLibrarySuffix = if isDarwin then "dylib" else "so";

          python = pkgsUnstable.python312;

          ladduWheel =
            {
              aarch64-darwin = {
                url = "https://files.pythonhosted.org/packages/aa/b2/3028e9899661d0220c4062ae13660c28bc806f2661033d8a2e0c663cf51b/laddu-0.21.5-cp311-abi3-macosx_11_0_arm64.whl";
                hash = "sha256-OjndxhVzzmZsSnB8qUeu6hMgIdmIX1fqmSCO7+mXPH8=";
              };
              x86_64-linux = {
                url = "https://files.pythonhosted.org/packages/cb/dd/63c72a64151ba8971bcb1e9bf4af9f96a4e1fcd0902d2dbcebdbbf85b933/laddu-0.21.5-cp311-abi3-manylinux_2_17_x86_64.manylinux2014_x86_64.whl";
                hash = "sha256-cmViAwLD3tqHiVw/Iz3ssVrymTTS6MqQAYWuqhrL1Zw=";
              };
              aarch64-linux = {
                url = "https://files.pythonhosted.org/packages/ee/cb/3287c12f71903788f7bfa92e2e5d9f133a67df848867f8a984cda1763d6b/laddu-0.21.5-cp311-abi3-manylinux_2_17_aarch64.manylinux2014_aarch64.whl";
                hash = "sha256-6CbOCtdZHuAuI6e6M7Pmex0WynTO6O9vFbvidFv8764=";
              };
            }
            .${system};

          laddu = python.pkgs.buildPythonPackage rec {
            pname = "laddu";
            version = "0.21.5";
            format = "wheel";

            src = pkgs.fetchurl {
              inherit (ladduWheel) url hash;
            };

            nativeBuildInputs = lib.optionals isLinux [ pkgs.autoPatchelfHook ];
            buildInputs = lib.optionals isLinux [ pkgs.stdenv.cc.cc.lib ];
            propagatedBuildInputs = [ python.pkgs.numpy ];
            pythonImportsCheck = [ "laddu" ];
            doCheck = false;
          };

          pythonEnvironment = python.withPackages (ps: [
            laddu
            ps.matplotlib
            ps.numpy
            ps.pyarrow
            ps.uproot
          ]);

          root = pkgs.root.override { inherit python; };

          boost170Base = pkgs.callPackage "${nixpkgs}/pkgs/development/libraries/boost/generic.nix" {
            version = "1.70.0";
            src = boost170-src;
          };

          boost170 = boost170Base.overrideAttrs (oldAttrs: {
            patches = builtins.filter (
              patch:
              !(lib.hasSuffix "cmake-paths-173.patch" (toString patch))
              && !(lib.hasSuffix "77f1e20bd69c2e7a9e25e6a9818ae6105f7d070c.patch" (toString patch))
            ) oldAttrs.patches;
            postPatch = (oldAttrs.postPatch or "") + ''
              substituteInPlace boost/container_hash/hash.hpp \
                --replace-fail "#if defined(_HAS_AUTO_PTR_ETC) && !_HAS_AUTO_PTR_ETC" \
                  "#if (defined(_HAS_AUTO_PTR_ETC) && !_HAS_AUTO_PTR_ETC) || defined(BOOST_NO_CXX98_FUNCTION_BASE)"
              substituteInPlace boost/log/utility/strictest_lock.hpp \
                --replace-fail "#include <boost/mpl/integral_c.hpp>" \
                  "#include <boost/type_traits/integral_constant.hpp>" \
                --replace-fail "mpl::integral_c< lock_access_mode," \
                  "boost::integral_constant< lock_access_mode,"
              substituteInPlace boost/thread/pthread/thread_data.hpp \
                --replace-fail "#if PTHREAD_STACK_MIN > 0" "#ifdef PTHREAD_STACK_MIN" \
                --replace-fail "if (size<PTHREAD_STACK_MIN) size=PTHREAD_STACK_MIN;" \
                  "if (size < static_cast<std::size_t>(PTHREAD_STACK_MIN)) size = static_cast<std::size_t>(PTHREAD_STACK_MIN);"
            '';
          });

          minuit2 = pkgs.stdenv.mkDerivation {
            pname = "minuit2-standalone";
            version = "unstable-2026-05-03";
            src = minuit2-src;
            nativeBuildInputs = [ pkgs.cmake ];
            postPatch = ''
              substituteInPlace src/MnPrintImpl.cxx \
                --replace-fail \
                  'const char *label[4] = {"[Error]", "[Warn]", "[Info]", "[Debug]"};' \
                  'const char *label[5] = {"[Error]", "[Warn]", "[Info]", "[Debug]", "[Trace]"};'
            '';
            cmakeFlags = [
              "-DCMAKE_BUILD_TYPE=Release"
              "-DBUILD_SHARED_LIBS=ON"
              "-DCMAKE_INSTALL_LIBDIR=lib"
            ];
            doInstallCheck = true;
            installCheckPhase = ''
              runHook preInstallCheck
              "$CXX" -std=c++17 -x c++ - -o minuit2-trace-smoke \
                -I"$out/include/Minuit2" -L"$out/lib" -lMinuit2 <<'EOF'
              #include "Minuit2/MnPrint.h"

              int main() {
                ROOT::Minuit2::MnPrint print("trace-smoke", 4);
                print.Trace("trace logging works");
              }
              EOF
              ./minuit2-trace-smoke
              runHook postInstallCheck
            '';
          };

          amptools = pkgs.stdenv.mkDerivation {
            pname = "amptools";
            version = "0.15.4";
            src = amptools-src;
            nativeBuildInputs = [
              pkgs.gnumake
              pkgs.pkg-config
            ];
            buildInputs = [ root ];
            ROOTSYS = root;

            buildPhase = ''
              runHook preBuild
              export PATH="${root}/bin:$PATH"
              export AMPTOOLS_HOME="$PWD"
              make -C AmpTools -j"$NIX_BUILD_CORES" CXX="$CXX"
              runHook postBuild
            '';

            installPhase = ''
              runHook preInstall
              mkdir -p "$out/include" "$out/lib" "$out/share/AmpTools"
              cp -r AmpTools/IUAmpTools "$out/include/"
              for optional_headers in DataIO Amplitude Utilities; do
                if [[ -d "AmpTools/$optional_headers" ]]; then
                  cp -r "AmpTools/$optional_headers" "$out/include/"
                fi
              done
              cp -r AmpTools/lib/. "$out/lib/"
              cp -r AmpTools Makefile.settings "$out/share/AmpTools/"
              runHook postInstall
            '';
          };

          buildInputs = [
            boost170
            minuit2
            pkgs.gsl
            root
            pkgs.log4cpp
          ];

          libraryPath = lib.makeLibraryPath (
            buildInputs ++ lib.optionals isLinux [ pkgs.stdenv.cc.cc.lib ]
          );

          commonEnvironment = {
            ROOTSYS = root;
            BOOSTROOT = boost170;
            BOOST_ROOT = boost170;
            BOOST_INCLUDEDIR = "${boost170.dev}/include";
            BOOST_LIBRARYDIR = "${boost170.out}/lib";
            BOOSTLIBPATH = "${boost170.out}/lib";
            BOOSTINCLUDE = "${boost170.dev}/include";
            MINUIT2_STANDALONE = minuit2;
            MINUIT2_LIBRARY = "${minuit2}/lib/libMinuit2.${sharedLibrarySuffix}";
          };

          pawian = pkgs.stdenv.mkDerivation (
            commonEnvironment
            // {
              pname = "pawian";
              version = "unstable";
              src = lib.cleanSource self;
              nativeBuildInputs = [
                pkgs.boost-build
                pkgs.pkg-config
              ] ++ root.nativeBuildInputs;
              propagatedBuildInputs = root.propagatedBuildInputs;
              inherit buildInputs;

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
                Scripts/run-pawian-build b2 \
                  -j"$NIX_BUILD_CORES" \
                  variant=release \
                  gammap-install
                runHook postBuild
              '';
              installPhase = ''
                runHook preInstall
                mkdir -p "$out"
                cp -r bin lib "$out/"
                runHook postInstall
              '';
            }
          );
        in
        {
          packages = {
            default = pawian;
            inherit pawian boost170 minuit2 amptools laddu;
          };

          checks = {
            development-shell-tools = pkgs.runCommand "development-shell-tools" {
              nativeBuildInputs = [ pythonEnvironment ];
            } ''
              python -c 'import importlib.metadata; assert importlib.metadata.version("laddu") == "0.21.0"'
              python -c 'import sys; assert sys.version_info[:2] == (3, 12)'
              touch "$out"
            '';
            pawian-gammap = pawian;
            inherit amptools;
          };

          devShells.default = pkgs.mkShell (
            commonEnvironment
            // {
              packages = [
                pkgs.boost-build
                pkgs.cmake
                pkgs.gnumake
                pkgs.just
                pkgs.ninja
                pkgs.pkg-config
                pythonEnvironment
                amptools
                root
                boost170
                minuit2
                pkgs.gsl
                pkgs.log4cpp
              ]
              ++ lib.optionals isDarwin [
                pkgs.clang
                pkgs.lldb
              ]
              ++ lib.optionals isLinux [
                pkgs.gcc
                pkgs.gdb
              ];

              AMPTOOLS_SRC = amptools-src;
              AMPTOOLS_PREFIX = amptools;
              AMPTOOLS_HOME = "${amptools}/share/AmpTools";
              NIX_PYTHON = "${pythonEnvironment}/bin/python";
              UV_NO_MANAGED_PYTHON = 1;

              shellHook = ''
                export TOP_DIR="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
                export BOOST_BUILD_PATH="${pkgs.boost-build}/share/boost-build"
                export PATH="${pythonEnvironment}/bin:$ROOTSYS/bin:$TOP_DIR/bin:$PATH"
                export PYTHONPATH="${root}/${python.sitePackages}:''${PYTHONPATH:-}"
                export EVT_DIR="''${EVT_DIR:-$TOP_DIR/Event}"
                export KMAT_DIR="''${KMAT_DIR:-$TOP_DIR/KMatStore}"

                if [[ "${if isDarwin then "1" else "0"}" == 1 ]]; then
                  export DYLD_FALLBACK_LIBRARY_PATH="$TOP_DIR/lib:${libraryPath}''${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
                  export PAWIAN_DEBUGGER=lldb
                else
                  export LD_LIBRARY_PATH="$TOP_DIR/lib:${libraryPath}''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
                  export PAWIAN_DEBUGGER=gdb
                fi

                echo "PAWIAN comparison development shell"
                echo "  platform=${system}"
                echo "  build: just build-comparison"
                echo "  demo:  just demo"
              '';
            }
          );
        }
      );
}
