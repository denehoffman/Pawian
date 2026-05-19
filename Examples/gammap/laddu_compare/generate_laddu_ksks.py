#!/usr/bin/env python3
"""Generate flat gamma p -> Kshort Kshort p samples with laddu."""

from __future__ import annotations

import argparse
from pathlib import Path

import laddu as ld
from common import (
    AUX_COLUMNS,
    MASS_KSHORT,
    MASS_PROTON,
    P4_COLUMNS,
)


def build_reaction(
    beam_energy: float,
    mass_min: float,
    mass_max: float,
    t_min: float,
    t_max: float,
) -> ld.GeneratedReaction:
    beam = ld.GeneratedParticle.initial(
        'beam',
        ld.InitialGenerator.beam_with_fixed_energy(0.0, beam_energy),
        ld.Reconstruction.stored(),
    )
    target = ld.GeneratedParticle.initial(
        'target',
        ld.InitialGenerator.target(MASS_PROTON),
        ld.Reconstruction.missing(),
    )
    kshort1 = ld.GeneratedParticle.stable(
        'kshort1',
        ld.StableGenerator(MASS_KSHORT),
        ld.Reconstruction.stored(),
    )
    kshort2 = ld.GeneratedParticle.stable(
        'kshort2',
        ld.StableGenerator(MASS_KSHORT),
        ld.Reconstruction.stored(),
    )
    kk = ld.GeneratedParticle.composite(
        'kk',
        ld.CompositeGenerator(mass_min, mass_max),
        (kshort1, kshort2),
        ld.Reconstruction.composite(),
    )
    recoil = ld.GeneratedParticle.stable(
        'recoil',
        ld.StableGenerator(MASS_PROTON),
        ld.Reconstruction.stored(),
    )
    t_edges = [t_min + (t_max - t_min) * i / 100 for i in range(101)]
    t_hist = ld.Histogram(t_edges, [1.0] * 100)
    return ld.GeneratedReaction.two_to_two(
        beam,
        target,
        kk,
        recoil,
        ld.MandelstamTDistribution.histogram(t_hist),  # NOTE: maybe use exponential
    )


def write_sample(
    path: Path,
    *,
    n_events: int,
    seed: int,
    beam_energy: float,
    mass_min: float,
    mass_max: float,
    pol_magnitude: float,
    pol_angle: float,
    t_min: float,
    t_max: float,
) -> None:
    reaction = build_reaction(beam_energy, mass_min, mass_max, t_min, t_max)
    generator = ld.EventGenerator(
        reaction,
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
    parser.add_argument('--beam-energy', type=float, default=8.0)
    parser.add_argument('--mass-min', type=float, default=1.0)
    parser.add_argument('--mass-max', type=float, default=2.0)
    parser.add_argument('--pol-magnitude', type=float, default=0.3519)
    parser.add_argument('--pol-angle', type=float, default=0.0)
    parser.add_argument('--t-min', type=float, default=-3.0)
    parser.add_argument('--t-max', type=float, default=-0.2)
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
            beam_energy=args.beam_energy,
            mass_min=args.mass_min,
            mass_max=args.mass_max,
            pol_magnitude=args.pol_magnitude,
            pol_angle=args.pol_angle,
            t_min=args.t_min,
            t_max=args.t_max,
        )
        print(f'wrote {args.out_dir / filename} ({n_events} events)')
    print(f'p4 columns: {", ".join(P4_COLUMNS)}')
    print(f'aux columns: {", ".join(AUX_COLUMNS)}')


if __name__ == '__main__':
    main()
