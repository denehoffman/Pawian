#!/usr/bin/env python3
"""Run the mass-binned S0+ + D2+ fit with laddu 0.21."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import json
import math
import multiprocessing
from pathlib import Path

import laddu as ld
import numpy as np
from model import MASS_EDGES, bin_model, channel


def fit_bin(index: int, args: argparse.Namespace) -> dict | None:
    """Load and fit one mass bin within an independent worker process."""
    reaction = channel()
    mass = reaction.mass('X')
    bins = ld.Bin(np.asarray(MASS_EDGES, dtype=np.float64))
    data_bins = ld.read_parquet(args.data).bin_by(mass, bins=bins)
    mc_bins = ld.read_parquet(args.mc).bin_by(mass, bins=bins)

    execution = ld.Execution(args.backend, threads=args.threads)
    for observed, normalization in zip(data_bins, mc_bins, strict=True):
        if observed.index != index:
            continue
        if len(observed.dataset) == 0 or len(normalization.dataset) == 0:
            return None
        prefix = f'bin_{observed.index:02d}'
        likelihood = ld.Likelihood(
            [
                ld.NLL(
                    bin_model(reaction, prefix),
                    data=observed.dataset,
                    accepted_mc=normalization.dataset,
                    name=prefix,
                )
            ],
            execution=execution,
        )
        fit = likelihood.fit(
            initial={
                f'{prefix}_d_real': 0.15,
                f'{prefix}_d_imaginary': 0.15,
            },
            # terminators=[ld.ganesh.MaxSteps(args.max_steps)],
        )
        names = fit.parameter_names or likelihood.parameter_names
        bin_parameters = dict(zip(names, np.asarray(fit.x, dtype=float), strict=True))
        real = bin_parameters[f'{prefix}_d_real']
        imaginary = bin_parameters[f'{prefix}_d_imaginary']
        return {
            'nll': fit.fx,
            'parameters': bin_parameters,
            'bin': {
                'index': observed.index,
                'low': observed.low,
                'high': observed.high,
                'data_events': len(observed.dataset),
                'mc_events': len(normalization.dataset),
                'd_magnitude': math.hypot(real, imaginary),
                'd_phase': math.atan2(imaginary, real),
            },
        }
    raise ValueError(f'mass bin {index} was not found')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('generated/data.parquet'))
    parser.add_argument('--mc', type=Path, default=Path('generated/mc.parquet'))
    parser.add_argument('--output', type=Path, default=Path('generated/laddu-fit.json'))
    parser.add_argument('--max-steps', type=int, default=500)
    parser.add_argument('--backend', choices=('cpu', 'jit'), default='jit')
    parser.add_argument('--jobs', type=int, default=1, help='Concurrent mass-bin processes')
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Laddu threads per mass-bin process (default: 1)',
    )
    args = parser.parse_args()
    if args.jobs < 1 or args.threads < 1:
        parser.error('--jobs and --threads must be positive')

    total_nll = 0.0
    parameters = {}
    fitted_bins = []
    # Spawn keeps native Laddu thread pools out of forked process state.
    with ProcessPoolExecutor(
        max_workers=min(args.jobs, len(MASS_EDGES) - 1),
        mp_context=multiprocessing.get_context('spawn'),
    ) as pool:
        for result in pool.map(fit_bin, range(len(MASS_EDGES) - 1), repeat(args)):
            if result is None:
                continue
            total_nll += result['nll']
            parameters.update(result['parameters'])
            fitted_bins.append(result['bin'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({'nll': total_nll, 'parameters': parameters, 'bins': fitted_bins}, indent=2) + '\n',
        encoding='utf-8',
    )
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
