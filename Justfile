set shell := ["bash", "-eu", "-o", "pipefail", "-c"]

root := justfile_directory()
jobs := env_var_or_default("JOBS", num_cpus())
variant := env_var_or_default("VARIANT", "release")
pawian_build := root / "Scripts/run-pawian-build"
amptools_runner := root / "Examples/gammap/comparison/run-amptools-reference"
comparison_justfile := root / "Examples/gammap/comparison/Justfile"

default:
    @just --justfile {{ quote(root / "Justfile") }} --list

# Show resolved native dependencies and actionable missing packages.
[group('environment')]
doctor:
    {{ quote(pawian_build) }} doctor

# Prepare the isolated Python comparison environment.
[group('environment')]
setup:
    just --justfile {{ quote(comparison_justfile) }} setup

# Build any Boost.Build target with the portable dependency adapter.
[group('build')]
build target variant=variant:
    {{ quote(pawian_build) }} b2 -j{{ quote(jobs) }} variant={{ quote(variant) }} {{ quote(target) }}

# Build a target with debug symbols.
[group('build')]
build-debug target:
    {{ quote(pawian_build) }} b2 -j{{ quote(jobs) }} variant=debug {{ quote(target) }}

# Build the gamma-proton PAWIAN executable.
[group('build')]
build-gammap:
    {{ quote(pawian_build) }} b2 -j{{ quote(jobs) }} variant={{ quote(variant) }} gammap-install

# Build the AmpTools comparison and reference executables.
[group('build')]
build-amptools:
    {{ quote(amptools_runner) }} build

# Build the PAWIAN and AmpTools executables used by the comparison.
[group('build')]
build-comparison: build-gammap build-amptools

# Run the fast deterministic workflow tests without compiling PAWIAN.
[group('test')]
test:
    {{ quote(amptools_runner) }} test

# Run the deterministic AmpTools reference end to end.
[group('test')]
amptools-reference:
    {{ quote(amptools_runner) }} run

# Run fast tests, build PAWIAN, and run the reference (which builds AmpTools).
[group('test')]
check: test build-gammap amptools-reference

# Run a comparison recipe; defaults to the complete fit workflow.
[group('comparison')]
comparison recipe="all":
    JOBS={{ quote(jobs) }} just --justfile {{ quote(comparison_justfile) }} {{ quote(recipe) }}

# Generate, fit, and plot the full comparison.
[group('comparison')]
demo:
    JOBS={{ quote(jobs) }} just --justfile {{ quote(comparison_justfile) }} plots

# Clean one Boost.Build target.
[group('maintenance')]
clean target:
    {{ quote(pawian_build) }} b2 --clean {{ quote(target) }}

# Clean all Boost.Build and local AmpTools executable outputs.
[group('maintenance')]
clean-all:
    {{ quote(pawian_build) }} b2 --clean-all
    make -C {{ quote(root / "Examples/gammap/comparison") }} clean
