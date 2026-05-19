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

The comparison-side laddu fit keeps the photon SDME contraction visible:

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

The `pawian-like` model currently uses a small direct random-search likelihood
driver around laddu evaluators. The standard `zlm` model can be selected with
`--model zlm` and uses laddu's normal `NLL` machinery.

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
