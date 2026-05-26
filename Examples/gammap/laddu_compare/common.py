#!/usr/bin/env python3
"""Shared setup for the Pawian/laddu gamma-p comparison scripts."""

from __future__ import annotations

from pathlib import Path

import laddu as ld

HERE = Path(__file__).resolve().parent
PAWIAN_ROOT = HERE.parents[2]

P4_COLUMNS = ['beam', 'target', 'kshort1', 'kshort2', 'recoil']
PAWIAN_FINAL_STATE = ['recoil', 'kshort1', 'kshort2']
AUX_COLUMNS = ['pol_magnitude', 'pol_angle']

MASS_PROTON = 0.938272
MASS_KSHORT = 0.497611


BEAM = ld.GeneratedParticle.initial(
    'beam',
    ld.InitialGenerator.beam_with_fixed_energy(0.0, 8.0),
    ld.Reconstruction.stored(),
)
TARGET = ld.GeneratedParticle.initial(
    'target',
    ld.InitialGenerator.target(MASS_PROTON),
    ld.Reconstruction.missing(),
)
KSHORT1 = ld.GeneratedParticle.stable(
    'kshort1',
    ld.StableGenerator(MASS_KSHORT),
    ld.Reconstruction.stored(),
)
KSHORT2 = ld.GeneratedParticle.stable(
    'kshort2',
    ld.StableGenerator(MASS_KSHORT),
    ld.Reconstruction.stored(),
)
RESONANCE = ld.GeneratedParticle.composite(
    'kk',
    ld.CompositeGenerator(1.0, 2.0),
    (KSHORT1, KSHORT2),
    ld.Reconstruction.composite(),
)
RECOIL = ld.GeneratedParticle.stable(
    'recoil',
    ld.StableGenerator(MASS_PROTON),
    ld.Reconstruction.stored(),
)
GEN_REACTION = ld.GeneratedReaction.two_to_two(
    BEAM,
    TARGET,
    RESONANCE,
    RECOIL,
    ld.MandelstamTDistribution.exponential(2.80),
)

REACTION = GEN_REACTION.reconstructed_reaction()
DECAY = REACTION.decay('kk')
POLARIZATION = REACTION.polarization(pol_magnitude='pol_magnitude', pol_angle='pol_angle')
RESONANCE_MASS = DECAY.parent_mass()


def open_dataset(path: Path) -> ld.Dataset:
    return ld.io.read_parquet(path, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
