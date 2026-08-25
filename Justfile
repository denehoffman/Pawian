set shell := ["zsh", "-c"]

root := justfile_directory()
b2_variant := "release"
b2_jobs := env_var_or_default("JOBS", `getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`)
b2_runner := root / "Scripts/run-nix-b2"

default:
    @just --list

build target variant=b2_variant:
    b2 -j{{ b2_jobs }} variant={{ variant }} {{ target }}

build-debug target:
    just build {{ target }} debug

clean target:
    b2 --clean {{ target }}

build-gammap:
    just build Examples/gammap//install

# Build the gamma-proton application with the pinned Nix dependencies.
nix-build-gammap:
    {{ b2_runner }} -j{{ b2_jobs }} variant=release gammap-install

# Build both executables needed by the local comparison.
build-comparison: nix-build-gammap build-amptools

# Build the pinned AmpTools package explicitly, then incrementally link fitZlm.
build-amptools:
    #!/usr/bin/env bash
    set -euo pipefail
    amptools_prefix="${AMPTOOLS_PREFIX:-$(nix build .#amptools --no-link --print-out-paths)}"
    make -C {{ root }}/Examples/gammap/comparison \
      -j{{ b2_jobs }} \
      AMPTOOLS_HOME="$amptools_prefix/share/AmpTools"

# Confirm that entering nix develop does not mutate the checkout.
check-development-shell:
    {{ root }}/Scripts/check-development-shell

# Run the current comparison workflow with the pinned build environment.
demo:
    #!/usr/bin/env bash
    set -euo pipefail
    amptools_prefix="${AMPTOOLS_PREFIX:-$(nix build .#amptools --no-link --print-out-paths)}"
    amptools_home="$amptools_prefix/share/AmpTools"
    amptools_lib="$amptools_home/AmpTools/lib"
    JOBS={{ b2_jobs }} \
      AMPTOOLS_HOME="$amptools_home" \
      LD_LIBRARY_PATH="$amptools_lib:${LD_LIBRARY_PATH:-}" \
      DYLD_FALLBACK_LIBRARY_PATH="$amptools_lib:${DYLD_FALLBACK_LIBRARY_PATH:-}" \
      {{ root }}/Scripts/with-nix-jamroot \
      just --justfile {{ root }}/Examples/gammap/comparison/Justfile plots

# Delegate to the self-contained laddu/PAWIAN/AmpTools comparison.
zlm-comparison recipe="all":
    #!/usr/bin/env bash
    set -euo pipefail
    amptools_prefix="${AMPTOOLS_PREFIX:-$(nix build .#amptools --no-link --print-out-paths)}"
    amptools_home="$amptools_prefix/share/AmpTools"
    amptools_lib="$amptools_home/AmpTools/lib"
    JOBS={{ b2_jobs }} \
      AMPTOOLS_HOME="$amptools_home" \
      LD_LIBRARY_PATH="$amptools_lib:${LD_LIBRARY_PATH:-}" \
      DYLD_FALLBACK_LIBRARY_PATH="$amptools_lib:${DYLD_FALLBACK_LIBRARY_PATH:-}" \
      {{ root }}/Scripts/with-nix-jamroot \
      just --justfile {{ root }}/Examples/gammap/comparison/Justfile {{ recipe }}
