#!/usr/bin/env python3
"""Make first-pass comparable mass plots from laddu data and fit weights."""

from __future__ import annotations

import argparse
import csv
import os
import pickle
import tempfile
from pathlib import Path

import laddu as ld
import numpy as np
from common import AUX_COLUMNS, P4_COLUMNS

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'matplotlib'))
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--genmc', type=Path, required=True)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=Path('ksks_mass_compare.png'))
    parser.add_argument('--bins', type=int, default=50)
    args = parser.parse_args()

    data = ld.io.read_parquet(args.data, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    genmc = ld.io.read_parquet(args.genmc, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    mass = ld.Mass(['kshort1', 'kshort2'])
    m_data = mass.value_on(data)
    m_genmc = mass.value_on(genmc)
    with args.fit.open('rb') as fit_file:
        fit = pickle.load(fit_file)  # noqa: S301
    weights = fit.get('genmc_weights')

    data_counts, edges = np.histogram(m_data, bins=args.bins, range=(1.0, 2.0))
    fit_counts, _ = np.histogram(
        m_genmc,
        weights=weights,
        bins=args.bins,
        range=(1.0, 2.0),
    )
    csv_path = args.out.with_suffix('.csv')
    with csv_path.open('w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['m_low', 'm_high', 'data', 'fit'])
        for low, high, data_count, fit_count in zip(
            edges[:-1],
            edges[1:],
            data_counts,
            fit_counts,
            strict=True,
        ):
            writer.writerow([low, high, data_count, fit_count])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(m_data, bins=args.bins, range=(1.0, 2.0), histtype='step', color='black', label='data')
    ax.hist(
        m_genmc,
        weights=weights,
        bins=args.bins,
        range=(1.0, 2.0),
        histtype='stepfilled',
        alpha=0.35,
        color='#0072b2',
        label='laddu fit on genmc',
    )
    ax.set_xlabel(r'$m(K_S K_S)$ [GeV]')
    ax.set_ylabel('counts')
    ax.legend()
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    print(f'wrote {args.out}')
    print(f'wrote {csv_path}')


if __name__ == '__main__':
    main()
