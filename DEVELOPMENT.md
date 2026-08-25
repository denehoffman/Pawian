# Reproducible development shell

The flake supports Apple-Silicon macOS (`aarch64-darwin`) and NixOS on
`x86_64-linux` or `aarch64-linux`. Determinate Nix is recommended on macOS;
upstream Nix is sufficient on NixOS. The private Minuit2 input is fetched over
SSH, so GitHub authentication must already work before Nix evaluates the
flake.

Enter the dependency-only shell from the repository root:

```bash
nix develop
```

Shell entry does not build PAWIAN or AmpTools and does not create `Jamroot`,
`KMatStore`, build products, or symlinks. The Nix-specific build definition is
`nix/Jamroot.nix`; legacy `Jamroot*` files and `SetEnv_*` scripts remain
available for manual environments.

The shell provides Python 3.12 with the released `laddu==0.20.0` wheel, ROOT,
Boost 1.70, Boost.Build, the authenticated Minuit2 fork, GSL, log4cpp, CMake,
Make, Just, plotting tools, and platform-native compiler/debugger tools. The
AmpTools 0.15.4 source is pinned as `AMPTOOLS_SRC`, but its library is built
only by an explicit build command.

## Incremental commands

```bash
just build-comparison  # PAWIAN gamma-proton app and fitZlm
just nix-build-gammap  # PAWIAN gamma-proton app only
just build-amptools    # pinned AmpTools plus fitZlm only
just demo              # build and run the current comparison workflow
```

The existing target-oriented recipes remain available for manual environments,
including `just build <target>` and `just build-gammap`. Nix recipes select
`nix/Jamroot.nix` only while their command runs and restore any pre-existing
manual `Jamroot` symlink afterward.

`just nix-build-gammap` writes normal incremental outputs under `bin/` and
`lib/`. `just build-amptools` asks Nix for the pinned AmpTools package without
creating a `result` symlink, then lets the comparison Makefile rebuild only
stale objects.

Set `JOBS` to override the portable CPU-count default. Set `AMPTOOLS_PREFIX` to
an already built AmpTools package prefix to skip `nix build .#amptools`.

## Verification

On each supported host, run:

```bash
nix flake check
just check-development-shell
just build-comparison
```

The shell check records both Git state and ignored legacy paths before and
after `nix develop`, and fails if shell entry changes the checkout.
