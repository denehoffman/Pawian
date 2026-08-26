#!/usr/bin/env python3
"""Translate laddu parquet files to PAWIAN ASCII and mass-binned AmpTools ROOT."""

from __future__ import annotations

import argparse
from pathlib import Path

import laddu as ld
import numpy as np
import uproot
from model import MASS_EDGES

P4_NAMES = ('beam', 'recoil', 'ks1', 'ks2')
FINAL_STATE = ('recoil', 'ks1', 'ks2')


def arrays(dataset: ld.Dataset) -> tuple[dict[str, np.ndarray], np.ndarray]:
    p4s = {}
    for name in P4_NAMES:
        p4 = ld.Vec4.event(name)
        p4s[name] = np.column_stack(
            [dataset.evaluate(component, real=True) for component in (p4.e(), p4.px(), p4.py(), p4.pz())]
        )
    return p4s, np.asarray(dataset.weights(), dtype=np.float64)


def write_ascii(path: Path, p4s: dict[str, np.ndarray], weights: np.ndarray) -> None:
    rows = np.column_stack([weights, *(p4s[name][:, [1, 2, 3, 0]] for name in FINAL_STATE)])
    np.savetxt(path, rows, fmt='%.17g')


def write_root_bins(directory: Path, stem: str, p4s: dict[str, np.ndarray], weights: np.ndarray) -> None:
    mass = np.sqrt(
        np.maximum(
            (p4s['ks1'][:, 0] + p4s['ks2'][:, 0]) ** 2
            - np.sum(
                (p4s['ks1'][:, 1:] + p4s['ks2'][:, 1:]) ** 2,
                axis=1,
            ),
            0.0,
        )
    )
    for index, (low, high) in enumerate(zip(MASS_EDGES[:-1], MASS_EDGES[1:], strict=True)):
        selected = (mass >= low) & (mass < high if index + 1 < len(MASS_EDGES) - 1 else mass <= high)
        n_events = int(np.count_nonzero(selected))
        final = [p4s[name][selected] for name in FINAL_STATE]
        branches = {
            'NumFinalState': np.full(n_events, len(FINAL_STATE), dtype=np.int32),
            'E_Beam': p4s['beam'][selected, 0].astype(np.float32),
            'Px_Beam': p4s['beam'][selected, 1].astype(np.float32),
            'Py_Beam': p4s['beam'][selected, 2].astype(np.float32),
            'Pz_Beam': p4s['beam'][selected, 3].astype(np.float32),
            'Weight': weights[selected].astype(np.float32),
            'E_FinalState': np.column_stack([item[:, 0] for item in final]).astype(np.float32),
            'Px_FinalState': np.column_stack([item[:, 1] for item in final]).astype(np.float32),
            'Py_FinalState': np.column_stack([item[:, 2] for item in final]).astype(np.float32),
            'Pz_FinalState': np.column_stack([item[:, 3] for item in final]).astype(np.float32),
        }
        with uproot.recreate(directory / f'{stem}-{index:02d}.root') as output:
            tree = output.mktree(
                'kin',
                {
                    'NumFinalState': 'int32',
                    'E_Beam': 'float32',
                    'Px_Beam': 'float32',
                    'Py_Beam': 'float32',
                    'Pz_Beam': 'float32',
                    'Weight': 'float32',
                    'E_FinalState': '3 * float32',
                    'Px_FinalState': '3 * float32',
                    'Py_FinalState': '3 * float32',
                    'Pz_FinalState': '3 * float32',
                },
            )
            tree.extend(branches)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('generated/data.parquet'))
    parser.add_argument('--mc', type=Path, default=Path('generated/mc.parquet'))
    parser.add_argument('--output', type=Path, default=Path('generated'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    for stem, source in (('data', args.data), ('mc', args.mc)):
        p4s, weights = arrays(ld.read_parquet(source))
        write_ascii(args.output / f'{stem}.txt', p4s, weights)
        write_root_bins(args.output, stem, p4s, weights)
    print(f'wrote PAWIAN and AmpTools inputs under {args.output}')


if __name__ == '__main__':
    main()
