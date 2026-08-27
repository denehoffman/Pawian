# Mass-binned Zlm comparison

This directory is a self-contained comparison of laddu 0.21, PAWIAN, and
AmpTools for

```text
gamma p -> X p,  X -> Kshort Kshort.
```

Laddu generates one toy sample containing an `f0(1500)` in `S0+` and an
`f2(1270)` in `D2+`, with a relative magnitude of 0.70 and phase of 0.65 rad.
The fit is deliberately mass-independent: ten 100 MeV bins from 1 to 2 GeV
each fit one complex `D2+/S0+` coefficient. The resonance line shapes occur
only in generation, so the three fit programs compare angular and
normalization conventions rather than line-shape implementations.

## Common model

Only positive reflectivity and the proton-nonflip sector are active. The two
incoherent components are

```text
I = |S sqrt(1+P) Re[Y00 exp(-i Phi)]
     + D sqrt(1+P) Re[Y22 exp(-i Phi)]|^2
  + |S sqrt(1-P) Im[Y00 exp(-i Phi)]
     + D sqrt(1-P) Im[Y22 exp(-i Phi)]|^2.
```

Here `P` and the laboratory polarization angle vary event by event, and `Phi`
is the angle to the production plane calculated in `Zlm.cc`. Laddu samples
`P` uniformly from `[0.2, 0.4)` and samples the angle in radians from a sparse
histogram whose four occupied one-degree bins are centered at 0, 45, 90, and
135 degrees. The `S` coefficient is fixed to `1 + 0i`; the magnitude and phase
of `D` are reported in each mass bin after fitting its Cartesian real and
imaginary components.

- `model.py` writes the same expression directly with laddu 0.21 and reads the
  generated polarization scalar columns in both generation and fitting.
- AmpTools uses `Zlm 0 0 +1 +1`/`Zlm 2 2 +1 +1` for the real sum and
  `Zlm 0 0 -1 -1`/`Zlm 2 2 -1 -1` for the imaginary sum. `constrain` shares
  each production coefficient between the two sums. The ROOT converter stores
  `P cos(angle)` and `P sin(angle)` in the beam x and y momentum fields consumed
  by Zlm's eventwise-polarization mode.
- PAWIAN uses `productionFormalism = Reflectivity`. Only
  `Reflectivity_f0(1500)_J0_M0_R+_K0` and
  `Reflectivity_f2(1270)_J2_M2_R+_K0` remain active; every other production
  coefficient and both decay couplings are fixed. Its channel normalization
  floats so that fixing the S-wave reference does not also fix the absolute
  event yield, and every available normalization-MC event is used in each bin.
  The ASCII converter writes each event's polarization magnitude and angle
  immediately after its weight.

## Running

From this directory (or with `just --justfile` from elsewhere):

```bash
just setup
just generate
just laddu
just amptools
just pawian
just plots
```

`just all` runs all three fits. `just plots` runs them and writes the comparison
figure and its source table to `generated/fit-comparison.png` and
`generated/fit-comparison.csv`. In addition to the fitted D/S ratio, the figure
shows the S0+ and D2+ event projections from each fit and from the generation
truth. The diagonal wave yields are normalized with the common weighted MC;
because S0+ and D2+ interfere, their projected counts need not add to the data
count. Override sample sizes with `DATA_EVENTS` and
`MC_EVENTS`, the parallel worker count with `JOBS`, and the managed Python
interpreter with `COMPARISON_PYTHON`. All three programs dispatch independent
mass-bin jobs, with up to `JOBS` processes running concurrently (default: all
available CPUs). Each Laddu job fits a single bin using one JIT thread;
the results are combined into `generated/laddu-fit.json` in bin order.
Fit output is retained under `generated/timings/*.log` instead
of being streamed to the terminal. At the end of `just plots` or `just demo`, a
compact table reports each framework's wall-clock time for fitting all ten bins;
generation, conversion, configuration, and build time are excluded.
`setup` creates `generated/python` and installs the released `laddu==0.21.0`;
inside the Nix shell it validates and uses `NIX_PYTHON` without modifying it.
Developers may point `COMPARISON_PYTHON` at an environment containing the
adjacent laddu checkout.

Generated parquet, ROOT, ASCII, configuration, normalization-integral, and fit
files stay under `generated/`. Without `AMPTOOLS_HOME` or `AMPTOOLS_PREFIX`, the
AmpTools runner bootstraps its pinned source into the ignored generated cache.
PAWIAN's portable build adapter resolves the active system or Nix dependencies.

## Frozen AmpTools S-wave reference

Run the smallest compatibility reference from the repository root:

```bash
just amptools-reference
```

This command does not invoke Nix. It requires a C++ compiler, Make, Python 3,
and a ROOT installation whose `root-config` is on `PATH`. On macOS, those
system prerequisites can be installed with `brew install root python just`.
On Linux, install ROOT and source its `thisroot.sh` before running the command.

The first run creates an ignored cache under
`generated/amptools-reference-dependencies/`. It installs the pinned minimal
Python dependencies there, downloads the SHA-256-protected AmpTools v0.15.4
source archive, and builds it against the active ROOT installation. Later runs
reuse that cache. Set `AMPTOOLS_HOME` to use an existing built source tree,
`AMPTOOLS_PREFIX` to use a packaged AmpTools prefix, or
`AMPTOOLS_REFERENCE_CACHE` to move the cache.

Nix remains an optional way to provide the prerequisites:

```bash
nix develop --command just amptools-reference
```

This verifies recorded SHA-256 digests for `Zlm.cc`, `Zlm.h`, `wignerD.cc`,
and `wignerD.h`; translates the committed four-event data and six-event
normalization samples; explicitly writes the AmpTools normalization cache;
then evaluates a fixed unit-coefficient S0+ amplitude at zero beam
polarization. The validator requires all intensities to equal AmpTools'
identical-particle result `1 / (8 pi)` within floating-point tolerance.

Auditable output is retained in `generated/amptools-reference/`:

- `input-metadata.json` records event identifiers, particle ordering, counts,
  weights, the sample digest, and Pawian/AmpTools revisions.
- `normalization.ni` is prepared in a separate step before evaluation.
- `normalization-audit.json` records the exact generated and accepted event
  identifiers and weights read by AmpTools.
- `result.json` records per-event kinematic checks and intensities plus
  normalization and revision metadata.

The compatibility manifest notes that this S-wave reference has no local
barrier-factor artifact. Any later reference that introduces one must add it
to the manifest before its result can be accepted.
