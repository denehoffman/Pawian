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


def open_dataset(path: Path) -> ld.Dataset:
    return ld.io.read_parquet(path, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
