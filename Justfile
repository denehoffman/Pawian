set shell := ["zsh", "-c"]

root := justfile_directory()
b2_variant := "release"
b2_jobs := `nproc`

default:
    @just --list

build target variant=b2_variant:
    b2 -j{{b2_jobs}} variant={{variant}} {{target}}

build-debug target:
    just build {{target}} debug

clean target:
    b2 --clean {{target}}

build-gammap:
    just build Examples/gammap//install

# Delegate to the self-contained laddu/PAWIAN/AmpTools comparison.
zlm-comparison recipe="all":
    just --justfile {{root}}/Examples/gammap/amptools/Justfile {{recipe}}
