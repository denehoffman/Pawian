#!/usr/bin/env python3
"""Run the mass-binned S0+ + D2+ fit with laddu 0.21."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import laddu as ld
import numpy as np
from model import MASS_EDGES, bin_model, channel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('generated/data.parquet'))
    parser.add_argument('--mc', type=Path, default=Path('generated/mc.parquet'))
    parser.add_argument('--output', type=Path, default=Path('generated/laddu-fit.json'))
    parser.add_argument('--max-steps', type=int, default=500)
    parser.add_argument(
        '--threads',
        type=int,
        default=None,
        help='JIT worker threads (default: all CPUs available to laddu)',
    )
    args = parser.parse_args()

    reaction = channel()
    mass = reaction.mass('X')
    bins = ld.Bin(np.asarray(MASS_EDGES, dtype=np.float64))
    data_bins = ld.read_parquet(args.data).bin_by(mass, bins=bins)
    mc_bins = ld.read_parquet(args.mc).bin_by(mass, bins=bins)

    execution = ld.Execution('jit', threads=args.threads, precision='f64')
    total_nll = 0.0
    parameters = {}
    fitted_bins = []
    for observed, normalization in zip(data_bins, mc_bins, strict=True):
        if len(observed.dataset) == 0 or len(normalization.dataset) == 0:
            continue
        prefix = f'bin_{observed.index:02d}'
        likelihood = ld.Likelihood(
            [ld.NLL(
                bin_model(reaction, prefix),
                data=observed.dataset,
                accepted_mc=normalization.dataset,
                name=prefix,
            )],
            execution=execution,
        )
        fits = [
            likelihood.fit(
                initial={
                    f'{prefix}_d_magnitude': magnitude,
                    f'{prefix}_d_phase': 0.35,
                },
                terminators=[ld.ganesh.MaxSteps(args.max_steps)],
                observers=[ld.ganesh.DebugObserver()],
            )
            for magnitude in (0.15, 0.5)
        ]
        fit = min(fits, key=lambda candidate: candidate.fx)
        names = fit.parameter_names or likelihood.parameter_names
        bin_parameters = dict(zip(names, np.asarray(fit.x, dtype=float), strict=True))
        parameters.update(bin_parameters)
        total_nll += fit.fx
        fitted_bins.append(
            {
                'index': observed.index,
                'low': observed.low,
                'high': observed.high,
                'data_events': len(observed.dataset),
                'mc_events': len(normalization.dataset),
                'd_magnitude': bin_parameters[f'{prefix}_d_magnitude'],
                'd_phase': bin_parameters[f'{prefix}_d_phase'],
            }
        )
        print(f'bin {observed.index}: {fit}')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({'nll': total_nll, 'parameters': parameters, 'bins': fitted_bins}, indent=2) + '\n',
        encoding='utf-8',
    )
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
