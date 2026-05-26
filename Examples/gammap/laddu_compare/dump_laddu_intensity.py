#!/usr/bin/env python3
"""Dump laddu model intensities with the same schema as Pawian dumpIntensity."""

from __future__ import annotations

import argparse
from pathlib import Path

import laddu as ld
import numpy as np
from common import AUX_COLUMNS, P4_COLUMNS, RESONANCE_MASS
from comparison_parameters import parameter_vector
from pawian_like_model import pawian_ls_model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--accmc', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()

    ds_data = ld.io.read_parquet(args.data, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    ds_accmc = ld.io.read_parquet(args.accmc, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    model = pawian_ls_model()
    evaluator_data = model.load(ds_data)
    evaluator_accmc = model.load(ds_accmc)
    parameters = np.asarray(parameter_vector(evaluator_data), dtype=float)
    intensities_data = evaluator_data.evaluate(parameters).real
    intensities_accmc = evaluator_accmc.evaluate(parameters).real
    masses_data = RESONANCE_MASS.value_on(ds_data)
    masses_accmc = RESONANCE_MASS.value_on(ds_accmc)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', encoding='utf-8') as output:
        output.write('sample,event,m_ksks,intensity,weight\n')
        for event_index, (mass, intensity, weight) in enumerate(
            zip(masses_data, intensities_data, ds_data.weights, strict=True)
        ):
            output.write(f'data,{event_index},{mass:.17g},{intensity:.17g},{float(weight):.17g}\n')
        for event_index, (mass, intensity, weight) in enumerate(
            zip(masses_accmc, intensities_accmc, ds_accmc.weights, strict=True)
        ):
            output.write(f'mc,{event_index},{mass:.17g},{intensity:.17g},{float(weight):.17g}\n')
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
