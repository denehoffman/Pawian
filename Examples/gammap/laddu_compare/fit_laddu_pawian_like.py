#!/usr/bin/env python3
"""Fit the generated Kshort sample with the explicit laddu comparison model."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import laddu as ld
import numpy as np
from common import AUX_COLUMNS, P4_COLUMNS
from pawian_like_model import intensity_model, zlm_baseline_model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--accmc', type=Path, required=True)
    parser.add_argument('--genmc', type=Path)
    parser.add_argument('--out', type=Path, default=Path('laddu_pawian_like_fit.pkl'))
    parser.add_argument('--model', choices=['pawian-like', 'zlm'], default='pawian-like')
    parser.add_argument('--niters', type=int, default=10)
    parser.add_argument('--seed', type=int, default=1234)
    args = parser.parse_args()

    data = ld.io.read_parquet(args.data, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    accmc = ld.io.read_parquet(args.accmc, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
    model = intensity_model() if args.model == 'pawian-like' else zlm_baseline_model()

    nll = ld.NLL(model, data, accmc)
    rng = np.random.default_rng(args.seed)
    best = None
    best_nll = np.inf
    for _i in range(args.niters):
        p0 = rng.uniform(-1.0, 1.0, len(nll.parameters.free))
        for ipar, par in enumerate(nll.parameters.free):
            if par.name.endswith('_width'):
                p0[ipar] = 0.1
        fit = nll.minimize(p0)
        if fit.fx < best_nll:
            best = fit
            best_nll = fit.fx
    if best is None:
        msg = 'all laddu fits failed'
        raise RuntimeError(msg)
    parameter_names = nll.parameters.names
    best_x = best.x

    result = {'model': args.model, 'best': best, 'parameters': parameter_names}
    if args.genmc:
        genmc = ld.io.read_parquet(args.genmc, p4s=P4_COLUMNS, aux=AUX_COLUMNS)
        gen_eval = model.load(genmc)
        result['genmc_weights'] = gen_eval.evaluate(best_x).real
    with args.out.open('wb') as output:
        pickle.dump(result, output)
    print(f'wrote {args.out}')
    print(best)


if __name__ == '__main__':
    main()
