# Pawian/Laddu Kshort-Kshort Comparison

This directory is a self-contained workspace for comparing Pawian's current
`gammap` channel to a laddu model for

```text
gamma p -> X p,  X -> Kshort Kshort.
```

The Pawian ASCII files contain only final-state particles, in this order:

```text
recoil proton, Kshort1, Kshort2
```

The incident photon and target proton are part of the channel definition. The
laddu parquet files keep `beam`, `target`, `kshort1`, `kshort2`, `recoil`, plus
the fixed `pol_magnitude` and `pol_angle` auxiliary columns.

## First Pass

From the Pawian checkout, using this environment or the pinned requirements:

```bash
python Examples/gammap/laddu_compare/generate_laddu_ksks.py
python Examples/gammap/laddu_compare/translate_laddu_to_pawian.py \
  Examples/gammap/laddu_compare/generated/data.parquet \
  Examples/gammap/laddu_compare/generated/ksks_data.txt
python Examples/gammap/laddu_compare/translate_laddu_to_pawian.py \
  Examples/gammap/laddu_compare/generated/accmc.parquet \
  Examples/gammap/laddu_compare/generated/ksks_mc.txt
```

Then use `Examples/gammap/gammapKshortKshort.cfg` from the Pawian checkout. It
points at these translated ASCII files and uses the same fixed beam
polarization as the generated laddu samples.

The comparison-side laddu models keep the photon SDME contraction visible.
For parity checks, prefer direct intensity projection over fitting:

```bash
just --justfile Examples/gammap/laddu_compare/Justfile compare-intensity
```

The Justfile uses deterministic parameters from
`comparison_parameters.py`; it does not run an NLL fit for the parity check.
The fit target is kept only for later toy studies after event-by-event
intensity parity is understood.

If you do want to run the old exploratory fit:

```bash
python Examples/gammap/laddu_compare/fit_laddu_pawian_like.py \
  --data Examples/gammap/laddu_compare/generated/data.parquet \
  --accmc Examples/gammap/laddu_compare/generated/accmc.parquet \
  --genmc Examples/gammap/laddu_compare/generated/genmc.parquet \
  --out Examples/gammap/laddu_compare/generated/laddu_pawian_like_fit.pkl
python Examples/gammap/laddu_compare/compare_outputs.py \
  --data Examples/gammap/laddu_compare/generated/data.parquet \
  --genmc Examples/gammap/laddu_compare/generated/genmc.parquet \
  --fit Examples/gammap/laddu_compare/generated/laddu_pawian_like_fit.pkl \
  --out Examples/gammap/laddu_compare/generated/ksks_mass_compare.png
```

To recreate the Python environment:

```bash
python -m pip install -r Examples/gammap/laddu_compare/requirements.txt
```

The standard `zlm` model can be selected with `--model zlm` to cross-check
against laddu's reflectivity-basis amplitude. The default Justfile model is
`pawian-ls`, which uses laddu's canonical LS factors for both the `X -> Kshort
Kshort` decay and an approximate `gamma p -> X p` production node.

## Intensity Parity

Before using rejection-sampled toys, compare event-by-event intensities on the
same flat MC sample. Pawian writes `gammapIntensity<name>.csv` in
`dumpIntensity` mode:

```bash
gammapReactionApp -c Examples/gammap/gammapKshortKshort.cfg \
  --mode dumpIntensity --name _ksks
```

Dump the matching laddu intensities for the accepted MC parquet:

```bash
python Examples/gammap/laddu_compare/dump_laddu_intensity.py \
  --input Examples/gammap/laddu_compare/generated/accmc.parquet \
  --sample mc \
  --model pawian-like \
  --out Examples/gammap/laddu_compare/generated/laddu_intensity_mc.csv
```

Then compare ratios and shapes:

```bash
python Examples/gammap/laddu_compare/compare_intensity_csv.py \
  --pawian gammapIntensity_ksks.csv \
  --laddu Examples/gammap/laddu_compare/generated/laddu_intensity_mc.csv \
  --out Examples/gammap/laddu_compare/generated/intensity_ratio.png
```

For parity work, the important result is a flat Pawian/laddu ratio. The
absolute normalization can differ.

The Justfile does not use Pawian's raw dumped defaults directly. The
`pawian-params` recipe runs `prepare_pawian_params.py` after
`dumpDefaultParams` and writes `generated/ksks_startParams.dat`. The default
`minimal-j12-minus` preset keeps only:

```text
L1S1/2_J1/2P-1C0Tof0(1500)_proton
L1S3/2_J1/2P-1C0Tof2'(1525)_proton
L3S5/2_J1/2P-1C0Tof2'(1525)_proton
f0(1500) -> Kshort1 Kshort2
f2'(1525) -> Kshort1 Kshort2
```

All other Pawian production magnitudes are set to zero before the intensity
dump. The retained magnitudes/phases are assigned from
`comparison_parameters.py`, and the laddu evaluator uses the same values by
parameter name.

The `pawian-ls` model uses laddu's two-to-two production view and its public
`Production.canonical_factor(...)` API for `gamma p -> X p`. It is closer to
Pawian than the older `pawian-like` scalar-production model, but it is still a
convention/proxy model until Pawian and laddu agree exactly on the production
frame and photon-proton coupling basis. It also follows the current Pawian cfg
and leaves the `f0(1500)` and `f2'(1525)` decays at unit `WoDynamics`; add
matched dynamics to Pawian and laddu before using lineshapes in this parity
test.

## What To Compare

Start with flat generated samples and hand-fixed or dumped Pawian parameters.
The first invariant should be event handling: row counts, four-vector order,
and `m(Kshort1 Kshort2)` should agree between the laddu parquet files and the
translated Pawian ASCII files.

After that, compare model evaluations on the same accepted MC sample. The
overall scale can float, but the event-by-event intensity ratio should be
constant once the same spin terms and parameter values are active in both
models.

Only then move to generated toy data from a non-flat intensity and compare fit
results and projected wave yields.
