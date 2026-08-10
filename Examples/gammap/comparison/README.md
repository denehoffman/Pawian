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

Here `P = 0.3519`, the laboratory polarization angle is zero, and `Phi` is the
angle to the production plane calculated in `Zlm.cc`. The `S` coefficient is
fixed to `1 + 0i`; the magnitude and phase of `D` float independently in each
mass bin.

- `model.py` writes the same expression directly with laddu 0.21.
- AmpTools uses `Zlm 0 0 +1 +1`/`Zlm 2 2 +1 +1` for the real sum and
  `Zlm 0 0 -1 -1`/`Zlm 2 2 -1 -1` for the imaginary sum. `constrain` shares
  each production coefficient between the two sums.
- PAWIAN uses `productionFormalism = Reflectivity`. Only
  `Reflectivity_f0(1500)_J0_M0_R+_K0` and
  `Reflectivity_f2(1270)_J2_M2_R+_K0` remain active; every other production
  coefficient and both decay couplings are fixed.

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
`MC_EVENTS`, the parallel worker count with `JOBS`, and the Python interpreter
with `PYTHON`. By default Laddu uses its multithreaded JIT backend and the other
two programs fit independent mass bins in parallel using all available CPUs.
`setup` installs the released `laddu==0.21.0`; developers may instead install
the adjacent checkout with `python -m pip install -e ../../../../laddu/python/laddu`.

Generated parquet, ROOT, ASCII, configuration, normalization-integral, and fit
files stay under `generated/`. AmpTools additionally requires `AMPTOOLS_HOME`;
PAWIAN uses this checkout's normal b2 environment.
