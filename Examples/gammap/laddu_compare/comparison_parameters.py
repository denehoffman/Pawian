#!/usr/bin/env python3
"""Shared deterministic parameters for Pawian/laddu intensity comparisons."""

from __future__ import annotations

import numpy as np

PARAMETER_VALUES = {
    'L1S1/2_J1/2P-1C0Tof0(1500)_protonMag': 0.8,
    'L1S1/2_J1/2P-1C0Tof0(1500)_protonPhi': 0.1,
    'L0S0_f0(1500)ToKshort1_Kshort2Mag': 1.0,
    'L0S0_f0(1500)ToKshort1_Kshort2Phi': 0.0,
    "L1S3/2_J1/2P-1C0Tof2'(1525)_protonMag": 0.5,
    "L1S3/2_J1/2P-1C0Tof2'(1525)_protonPhi": -0.2,
    "L3S5/2_J1/2P-1C0Tof2'(1525)_protonMag": 0.3,
    "L3S5/2_J1/2P-1C0Tof2'(1525)_protonPhi": 0.4,
    "L2S0_f2'(1525)ToKshort1_Kshort2Mag": 1.0,
    "L2S0_f2'(1525)ToKshort1_Kshort2Phi": 0.0,
    'channelType4protonKshort1Kshort2Scaling': 1.0,
}


def parameter_vector(evaluator) -> np.typing.NDArray[np.float64]:
    values = []
    for parameter in evaluator.parameters.free:
        if parameter.name not in PARAMETER_VALUES:
            msg = f'no deterministic comparison value configured for laddu parameter {parameter.name!r}'
            raise KeyError(msg)
        values.append(PARAMETER_VALUES[parameter.name])
    return np.asarray(values, dtype=float)


def reflectivity_parameters() -> dict[str, float]:
    pass
