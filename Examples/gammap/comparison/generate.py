#!/usr/bin/env python3
"""Generate common laddu 0.21 pseudo-data and normalization MC."""

from __future__ import annotations

import argparse
from pathlib import Path

import laddu as ld
from model import channel, polarization_sources, truth_model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('generated'))
    parser.add_argument('--data-events', type=int, default=10_000)
    parser.add_argument('--mc-events', type=int, default=100_000)
    parser.add_argument('--seed', type=int, default=20260809)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    reaction = channel()
    generator = ld.Generator(reaction, scalars=polarization_sources())
    execution = ld.Execution('cpu')

    data, data_report = generator.unweighted(
        args.data_events,
        truth_model(reaction),
        execution=execution,
        seed=args.seed,
        pilot_proposals=50_000,
        safety_factor=2.0,
        grow_envelope=True,
        max_proposals=100_000_000,
    )
    normalization, mc_report = generator.unweighted(
        args.mc_events,
        execution=execution,
        seed=args.seed + 1,
        proven_envelope=True,
    )
    data.write_to(ld.ParquetSink(args.output / 'data.parquet'))
    normalization.write_to(ld.ParquetSink(args.output / 'mc.parquet'))
    print(f'data: {data_report}')
    print(f'normalization MC: {mc_report}')


if __name__ == '__main__':
    main()
