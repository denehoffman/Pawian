#!/usr/bin/env python3
"""Generate flat gamma p -> Kshort Kshort p samples with laddu."""

from __future__ import annotations

import argparse
from pathlib import Path

import laddu as ld
from common import AUX_COLUMNS, GEN_REACTION, P4_COLUMNS


def write_sample(
    path: Path,
    *,
    n_events: int,
    seed: int,
    pol_magnitude: float,
    pol_angle: float,
) -> None:
    generator = ld.EventGenerator(
        GEN_REACTION,
        aux_generators={
            'pol_magnitude': ld.Distribution.fixed(pol_magnitude),
            'pol_angle': ld.Distribution.fixed(pol_angle),
        },
        seed=seed,
        storage=ld.GeneratedStorage.only(tuple(P4_COLUMNS)),
    )
    dataset = generator.generate_dataset(n_events)
    ld.io.write_parquet(dataset, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, default=Path(__file__).parent / 'generated')
    parser.add_argument('--data-events', type=int, default=20_000)
    parser.add_argument('--accmc-events', type=int, default=100_000)
    parser.add_argument('--genmc-events', type=int, default=100_000)
    parser.add_argument('--pol-magnitude', type=float, default=0.3519)
    parser.add_argument('--pol-angle', type=float, default=0.0)
    parser.add_argument('--seed', type=int, default=12345)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    samples = [
        ('data.parquet', args.data_events, args.seed),
        ('accmc.parquet', args.accmc_events, args.seed + 1),
        ('genmc.parquet', args.genmc_events, args.seed + 2),
    ]
    for filename, n_events, seed in samples:
        write_sample(
            args.out_dir / filename,
            n_events=n_events,
            seed=seed,
            pol_magnitude=args.pol_magnitude,
            pol_angle=args.pol_angle,
        )
        print(f'wrote {args.out_dir / filename} ({n_events} events)')
    print(f'p4 columns: {", ".join(P4_COLUMNS)}')
    print(f'aux columns: {", ".join(AUX_COLUMNS)}')


if __name__ == '__main__':
    main()
