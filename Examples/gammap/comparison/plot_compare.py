#!/usr/bin/env python3
"""Compare the mass-binned D/S fit results from Laddu, AmpTools, and PAWIAN."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import laddu as ld
import matplotlib.pyplot as plt
import numpy as np
from model import MASS_EDGES, bin_model, bin_wave_models, channel, truth_model, truth_wave_models

MIN_PARAMETER_FIELDS = 2


@dataclass(frozen=True)
class Point:
    program: str
    index: int
    magnitude: float
    phase: float

    @property
    def low(self) -> float:
        return MASS_EDGES[self.index]

    @property
    def high(self) -> float:
        return MASS_EDGES[self.index + 1]

    @property
    def center(self) -> float:
        return (self.low + self.high) / 2.0


@dataclass(frozen=True)
class Projection:
    program: str
    index: int
    s0_events: float
    d2_events: float
    data_events: float

    @property
    def center(self) -> float:
        return (MASS_EDGES[self.index] + MASS_EDGES[self.index + 1]) / 2.0


def read_laddu(path: Path) -> list[Point]:
    payload = json.loads(path.read_text(encoding='utf-8'))
    return [
        Point('Laddu', int(item['index']), float(item['d_magnitude']), float(item['d_phase']))
        for item in payload['bins']
    ]


def read_amptools_file(path: Path) -> dict[str, float]:
    lines = path.read_text(encoding='utf-8').splitlines()
    header = '+++ Parameter Values and Errors +++'
    try:
        start = lines.index(header) + 1
    except ValueError as error:
        message = f'{path} has no AmpTools parameter section'
        raise ValueError(message) from error
    count = int(lines[start].strip())
    parameters = {}
    for line in lines[start + 1 : start + 1 + count]:
        name, value = line.split()[:2]
        parameters[name] = float(value)
    return parameters


def indexed_path(path: Path) -> int:
    match = re.search(r'(\d+)$', path.stem)
    if match is None:
        match = re.search(r'(\d+)$', path.parent.name)
    if match is None:
        message = f'cannot determine mass-bin index from {path}'
        raise ValueError(message)
    return int(match.group(1))


def read_amptools(directory: Path) -> list[Point]:
    points = []
    paths = sorted(directory.glob('zlm-bin-*.fit'))
    if not paths:
        message = f'no AmpTools fit files found under {directory}'
        raise FileNotFoundError(message)
    for path in paths:
        parameters = read_amptools_file(path)
        s_real_name = next(
            (name for name in parameters if name.endswith('::PosRe::S0+_re')),
            None,
        )
        s_imag_name = next(
            (name for name in parameters if name.endswith('::PosRe::S0+_im')),
            None,
        )
        real_name = next(
            (name for name in parameters if name.endswith('::PosRe::D2+_re')),
            None,
        )
        imag_name = next(
            (name for name in parameters if name.endswith('::PosRe::D2+_im')),
            None,
        )
        if None in (s_real_name, s_imag_name, real_name, imag_name):
            message = f'{path} is missing an S0+ or D2+ production coefficient'
            raise ValueError(message)
        s = complex(parameters[s_real_name], parameters[s_imag_name])
        if abs(s) == 0.0:
            message = f'{path} has a zero S0+ production coefficient'
            raise ValueError(message)
        d = complex(parameters[real_name], parameters[imag_name])
        ratio = d / s
        points.append(
            Point(
                'AmpTools',
                indexed_path(path),
                abs(ratio),
                math.atan2(ratio.imag, ratio.real),
            )
        )
    return points


def read_pawian(directory: Path) -> list[Point]:
    magnitude_name = 'Reflectivity_f2(1270)_J2_M2_R+_K0Mag'
    phase_name = 'Reflectivity_f2(1270)_J2_M2_R+_K0Phi'
    points = []
    paths = sorted(directory.glob('bin-*/finalResult_*.dat'))
    if not paths:
        message = f'no PAWIAN result files found under {directory}'
        raise FileNotFoundError(message)
    for path in paths:
        parameters = {
            fields[0]: float(fields[1])
            for line in path.read_text(encoding='utf-8').splitlines()
            if len(fields := line.split()) >= MIN_PARAMETER_FIELDS
        }
        try:
            magnitude, phase = parameters[magnitude_name], parameters[phase_name]
        except KeyError as error:
            message = f'{path} is missing the D2+ production coefficient'
            raise ValueError(message) from error
        points.append(Point('PAWIAN', indexed_path(path), magnitude, phase))
    return points


def weighted_sum(values: np.ndarray, dataset: ld.Dataset) -> float:
    return float(np.dot(values, dataset.weights()))


def project_waves(points: list[Point], data_path: Path, mc_path: Path) -> list[Projection]:
    """Project fitted and generated diagonal wave intensities onto mass."""
    reaction = channel()
    bins = ld.Bin(np.asarray(MASS_EDGES, dtype=np.float64))
    mass = reaction.mass('X')
    data_bins = ld.read_parquet(data_path).bin_by(mass, bins=bins)
    mc_bins = ld.read_parquet(mc_path).bin_by(mass, bins=bins)
    unit_s, unit_d = bin_wave_models(reaction)
    truth_s, truth_d = truth_wave_models(reaction)
    generated = truth_model(reaction)
    by_bin = {index: [] for index in range(len(MASS_EDGES) - 1)}
    for point in points:
        by_bin[point.index].append(point)

    projections = []
    for observed, normalization in zip(data_bins, mc_bins, strict=True):
        index = observed.index
        dataset = normalization.dataset
        data_events = observed.dataset.sum_weights()
        s_integral = weighted_sum(unit_s.evaluate(dataset, real=True), dataset)
        d_integral = weighted_sum(unit_d.evaluate(dataset, real=True), dataset)
        for point in by_bin[index]:
            prefix = f'projection_{point.program.lower()}_{index:02d}'
            fitted = bin_model(reaction, prefix)
            parameters = {
                f'{prefix}_d_magnitude': point.magnitude,
                f'{prefix}_d_phase': point.phase,
            }
            total = weighted_sum(fitted.evaluate(dataset, parameters=parameters, real=True), dataset)
            projections.append(
                Projection(
                    point.program,
                    index,
                    data_events * s_integral / total,
                    data_events * point.magnitude**2 * d_integral / total,
                    data_events,
                )
            )

        truth_total = weighted_sum(generated.evaluate(dataset, real=True), dataset)
        projections.append(
            Projection(
                'Truth',
                index,
                data_events * weighted_sum(truth_s.evaluate(dataset, real=True), dataset) / truth_total,
                data_events * weighted_sum(truth_d.evaluate(dataset, real=True), dataset) / truth_total,
                data_events,
            )
        )
    return projections


def write_table(path: Path, points: list[Point], projections: list[Projection]) -> None:
    projected = {(item.program, item.index): item for item in projections}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(
            (
                'program',
                'bin',
                'mass_low_gev',
                'mass_high_gev',
                'd_magnitude',
                'd_phase_rad',
                's0_projected_events',
                'd2_projected_events',
                'data_events',
            )
        )
        for point in points:
            item = projected[point.program, point.index]
            writer.writerow(
                (
                    point.program,
                    point.index,
                    point.low,
                    point.high,
                    point.magnitude,
                    point.phase,
                    item.s0_events,
                    item.d2_events,
                    item.data_events,
                )
            )
        for item in (projection for projection in projections if projection.program == 'Truth'):
            writer.writerow(
                (
                    item.program,
                    item.index,
                    MASS_EDGES[item.index],
                    MASS_EDGES[item.index + 1],
                    '',
                    '',
                    item.s0_events,
                    item.d2_events,
                    item.data_events,
                )
            )


def plot(path: Path, points: list[Point], projections: list[Projection]) -> None:
    figure, axes = plt.subplots(2, 2, sharex='col', figsize=(12, 8), constrained_layout=True)
    markers = {'Laddu': 'o', 'AmpTools': 's', 'PAWIAN': '^'}
    for program, marker in markers.items():
        selected = sorted((point for point in points if point.program == program), key=lambda point: point.index)
        mass = [point.center for point in selected]
        axes[0, 0].plot(mass, [point.magnitude for point in selected], marker=marker, label=program)
        axes[1, 0].plot(mass, [point.phase for point in selected], marker=marker, label=program)
        projected = sorted(
            (item for item in projections if item.program == program),
            key=lambda item: item.index,
        )
        axes[0, 1].plot(
            [item.center for item in projected],
            [item.s0_events for item in projected],
            marker=marker,
            label=program,
        )
        axes[1, 1].plot(
            [item.center for item in projected],
            [item.d2_events for item in projected],
            marker=marker,
            label=program,
        )

    truth = sorted((item for item in projections if item.program == 'Truth'), key=lambda item: item.index)
    data_events = [item.data_events for item in truth]
    for axis, attribute in ((axes[0, 1], 's0_events'), (axes[1, 1], 'd2_events')):
        axis.stairs(data_events, MASS_EDGES, color='0.75', linewidth=1.5, label='Data total')
        axis.plot(
            [item.center for item in truth],
            [getattr(item, attribute) for item in truth],
            color='black',
            linestyle='--',
            marker='.',
            label='Truth',
        )

    axes[0, 0].set_ylabel(r'$|D2^+/S0^+|$')
    axes[1, 0].set_ylabel(r'$\arg(D2^+/S0^+)$ [rad]')
    axes[0, 1].set_ylabel(r'Projected $S0^+$ events / bin')
    axes[1, 1].set_ylabel(r'Projected $D2^+$ events / bin')
    axes[1, 0].set_xlabel(r'$m(K_S K_S)$ [GeV]')
    axes[1, 1].set_xlabel(r'$m(K_S K_S)$ [GeV]')
    axes[1, 0].set_ylim(-math.pi, math.pi)
    axes[1, 0].set_yticks(
        (-math.pi, -math.pi / 2, 0.0, math.pi / 2, math.pi),
        (r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'),
    )
    for axis in axes.flat:
        axis.grid(alpha=0.25)
        axis.legend()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('generated'))
    parser.add_argument('--output', type=Path, default=Path('generated/fit-comparison.png'))
    args = parser.parse_args()
    points = [
        *read_laddu(args.directory / 'laddu-fit.json'),
        *read_amptools(args.directory / 'amptools'),
        *read_pawian(args.directory / 'pawian'),
    ]
    projections = project_waves(
        points,
        args.directory / 'data.parquet',
        args.directory / 'mc.parquet',
    )
    write_table(args.output.with_suffix('.csv'), points, projections)
    plot(args.output, points, projections)
    print(f'wrote {args.output} and {args.output.with_suffix(".csv")}')


if __name__ == '__main__':
    main()
