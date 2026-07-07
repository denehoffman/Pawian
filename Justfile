set shell := ["zsh", "-c"]
export MPLBACKEND := "Agg"

root := justfile_directory()
python := root + "/.venv/bin/python"

laddu_dir := root + "/Examples/gammap/laddu_compare"
generated := laddu_dir + "/generated"

cfg := root + "/Examples/gammap/gammapKshortKshort.cfg"
suffix := "_ksks"
params := generated + "/ksks_startParams.dat"

pawian_param_preset := "minimal-j12-minus"
model := "pawian-ls"

data_events := "200"
accmc_events := "1000"
genmc_events := "1000"
niters := "3"

b2_variant := "release"
b2_jobs := `nproc`

default:
    @just --list

build target:
  b2 -j{{b2_jobs}} variant={{b2_variant}} {{target}}

build-debug target:
    just build {{target}} debug

clean target:
    b2 --clean {{target}}

build-gammap:
    just build Examples/gammap//install

generate:
    {{python}} {{laddu_dir}}/generate_laddu_ksks.py \
      --out-dir {{generated}} \
      --data-events {{data_events}} \
      --accmc-events {{accmc_events}} \
      --genmc-events {{genmc_events}}

translate-data:
    {{python}} {{laddu_dir}}/translate_laddu_to_pawian.py \
      {{generated}}/data.parquet \
      {{generated}}/ksks_data.txt

translate-mc:
    {{python}} {{laddu_dir}}/translate_laddu_to_pawian.py \
      {{generated}}/accmc.parquet \
      {{generated}}/ksks_mc.txt

translate: translate-data translate-mc

pawian-params: build-gammap
    {{root}}/bin/gammapReactionApp -c {{cfg}} \
      --mode dumpDefaultParams --name {{suffix}}

    {{python}} {{laddu_dir}}/prepare_pawian_params.py \
      --input {{root}}/defaultparams{{suffix}}.dat \
      --out {{params}} \
      --preset {{pawian_param_preset}}

pawian-intensity: pawian-params translate
    {{root}}/bin/gammapReactionApp -c {{cfg}} \
      --mode dumpIntensity --name {{suffix}}

laddu-intensity:
    {{python}} {{laddu_dir}}/dump_laddu_intensity.py \
      --data {{generated}}/data.parquet \
      --accmc {{generated}}/accmc.parquet \
      --out {{generated}}/laddu_intensity_{{model}}.csv

compare-intensity: pawian-intensity laddu-intensity
    {{python}} {{laddu_dir}}/compare_intensity_csv.py \
      --pawian {{root}}/gammapIntensity{{suffix}}.csv \
      --laddu {{generated}}/laddu_intensity_{{model}}.csv \
      --out {{generated}}/intensity_ratio_{{model}}.png

laddu-analysis: generate translate compare-intensity

parity: generate compare-intensity

all: parity
