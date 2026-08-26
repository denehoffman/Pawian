#!/usr/bin/env python3
"""Prepare and protect the deterministic AmpTools S-wave reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[2]
DEFAULT_MANIFEST = HERE / 'reference/compatibility-artifacts.json'
DEFAULT_SAMPLE = HERE / 'reference/events.json'
PARTICLE_ORDER = ('beam', 'recoil', 'ks1', 'ks2')
FINAL_STATE = PARTICLE_ORDER[1:]
TARGET_PROTON_MASS_GEV = 0.9382720813
EPS = 2e-12
BIG_EPS = 2e-6


class DataReferenceError(RuntimeError):
    """An input cannot be accepted as part of the compatibility reference."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError) as error:
        msg = f'cannot read {path}: {error}'
        raise DataReferenceError(msg) from error
    if not isinstance(value, dict):
        msg = f'{path} must contain a JSON object'
        raise DataReferenceError(msg)
    return value


def verify_artifacts(manifest_path: Path, root: Path) -> int:
    manifest = load_json(manifest_path)
    artifacts = manifest.get('artifacts')
    if not isinstance(artifacts, list) or not artifacts:
        msg = f'{manifest_path} has no compatibility artifacts'
        raise DataReferenceError(msg)

    failures = []
    for entry in artifacts:
        relative = Path(str(entry.get('path', '')))
        expected = str(entry.get('sha256', ''))
        artifact = root / relative
        if not artifact.is_file():
            failures.append(f'{relative}: missing')
            continue
        actual = sha256(artifact)
        if actual != expected:
            failures.append(f'{relative}: digest mismatch (expected {expected}, got {actual})')
    if failures:
        raise DataReferenceError('compatibility artifact verification failed:\n  ' + '\n  '.join(failures))

    print(f'verified {len(artifacts)} compatibility artifacts')
    return len(artifacts)


def four_momentum_sum(events: list[list[float]]) -> list[float]:
    return [sum(momentum[index] for momentum in events) for index in range(4)]


def validate_sample(sample: dict[str, Any], path: Path) -> None:  # noqa: PLR0912
    if sample.get('format_version') != 1:
        msg = f'{path}: unsupported format_version'
        raise DataReferenceError(msg)
    if tuple(sample.get('particle_order', [])) != PARTICLE_ORDER:
        msg = f'{path}: particle_order must be {list(PARTICLE_ORDER)}'
        raise DataReferenceError(msg)
    if sample.get('polarization_fraction') != 0.0:
        msg = f'{path}: reference must have zero beam polarization'
        raise DataReferenceError(msg)

    beam_energy = float(sample.get('beam_energy_gev', math.nan))
    if not math.isfinite(beam_energy) or beam_energy <= 0.0:
        msg = f'{path}: invalid beam_energy_gev'
        raise DataReferenceError(msg)

    all_ids: set[int] = set()
    expected_total = [
        beam_energy + TARGET_PROTON_MASS_GEV,
        0.0,
        0.0,
        beam_energy,
    ]
    for collection_name in ('data', 'normalization'):
        collection = sample.get(collection_name)
        if not isinstance(collection, list) or not collection:
            msg = f'{path}: {collection_name} must be non-empty'
            raise DataReferenceError(msg)
        for event in collection:
            event_id = int(event.get('event_id'))
            if event_id in all_ids:
                msg = f'{path}: duplicate event_id {event_id}'
                raise DataReferenceError(msg)
            all_ids.add(event_id)
            weight = float(event.get('weight', math.nan))
            if not math.isfinite(weight) or weight <= 0.0:
                msg = f'{path}: event {event_id} has invalid weight'
                raise DataReferenceError(msg)
            particles = event.get('particles')
            if not isinstance(particles, dict) or tuple(particles) != PARTICLE_ORDER:
                msg = f'{path}: event {event_id} particle ordering is not explicit'
                raise DataReferenceError(msg)
            momenta = []
            for name in PARTICLE_ORDER:
                momentum = particles.get(name)
                if (
                    not isinstance(momentum, list)
                    or len(momentum) != 4  # noqa: PLR2004
                    or not all(math.isfinite(float(value)) for value in momentum)
                ):
                    msg = f'{path}: event {event_id} has invalid {name} four-momentum'
                    raise DataReferenceError(msg)
                momenta.append([float(value) for value in momentum])
            if momenta[0] != [beam_energy, 0.0, 0.0, beam_energy]:
                msg = f'{path}: event {event_id} beam is inconsistent'
                raise DataReferenceError(msg)
            total = four_momentum_sum(momenta[1:])
            if any(abs(actual - expected) > EPS for actual, expected in zip(total, expected_total)):
                msg = f'{path}: event {event_id} violates four-momentum conservation'
                raise DataReferenceError(msg)


def write_root(path: Path, events: list[dict[str, Any]]) -> None:
    try:
        import numpy as np  # noqa: PLC0415
        import uproot  # noqa: PLC0415
    except ImportError as error:
        msg = (
            'preparation requires numpy and uproot; run through '
            'run-amptools-reference or install the reference dependencies'
        )
        raise DataReferenceError(msg) from error

    four_momenta = {
        name: np.asarray([event['particles'][name] for event in events], dtype=np.float64) for name in PARTICLE_ORDER
    }
    final_state_momenta = [four_momenta[name] for name in FINAL_STATE]
    branches = {
        'EventID': np.asarray([event['event_id'] for event in events], dtype=np.int64),
        'NumFinalState': np.full(len(events), len(FINAL_STATE), dtype=np.int32),
        'E_Beam': four_momenta['beam'][:, 0].astype(np.float32),
        'Px_Beam': four_momenta['beam'][:, 1].astype(np.float32),
        'Py_Beam': four_momenta['beam'][:, 2].astype(np.float32),
        'Pz_Beam': four_momenta['beam'][:, 3].astype(np.float32),
        'Weight': np.asarray([event['weight'] for event in events], dtype=np.float32),
        'E_FinalState': np.column_stack([item[:, 0] for item in final_state_momenta]).astype(np.float32),
        'Px_FinalState': np.column_stack([item[:, 1] for item in final_state_momenta]).astype(np.float32),
        'Py_FinalState': np.column_stack([item[:, 2] for item in final_state_momenta]).astype(np.float32),
        'Pz_FinalState': np.column_stack([item[:, 3] for item in final_state_momenta]).astype(np.float32),
    }
    with uproot.recreate(path) as output:
        tree = output.mktree(
            'kin',
            {
                'EventID': 'int64',
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


def repository_revision() -> str:
    completed = subprocess.run(
        ['/usr/bin/git', 'describe', '--always', '--dirty'],
        cwd=REPOSITORY,
        check=True,
        text=True,
        capture_output=True,
    )
    return completed.stdout.strip()


def config_text(output: Path, *, use_prepared_normalization: bool) -> str:
    normalization = (output / 'normalization.ni').resolve()
    input_marker = ' input' if use_prepared_normalization else ''
    return f"""fit amptools-unpolarized-s-wave-reference

reaction reference Beam Proton KShort KShort
data reference ROOTDataReader {(output / 'data.root').resolve()}
genmc reference ROOTDataReader {(output / 'normalization.root').resolve()}
accmc reference ROOTDataReader {(output / 'normalization.root').resolve()}
normintfile reference {normalization}{input_marker}

sum reference PositiveReal
# Four arguments read the zero polarization vector from beam Px/Py.  The
# six-argument Zlm form treats a zero fraction as a histogram lookup request.
amplitude reference::PositiveReal::S0+ Zlm 0 0 1 1
initialize reference::PositiveReal::S0+ cartesian 1 0 fixed
"""


def write_configs(output: Path) -> None:
    (output / 'prepare-normalization.cfg').write_text(
        config_text(output, use_prepared_normalization=False), encoding='utf-8'
    )
    (output / 'reference.cfg').write_text(config_text(output, use_prepared_normalization=True), encoding='utf-8')


def prepare(sample_path: Path, output: Path) -> None:
    verify_artifacts(DEFAULT_MANIFEST, REPOSITORY)
    sample = load_json(sample_path)
    validate_sample(sample, sample_path)
    output.mkdir(parents=True, exist_ok=True)

    data = sample['data']
    normalization = sample['normalization']
    write_root(output / 'data.root', data)
    write_root(output / 'normalization.root', normalization)
    write_configs(output)

    manifest = load_json(DEFAULT_MANIFEST)
    metadata = {
        'format_version': 1,
        'reaction': sample['reaction'],
        'particle_order': list(PARTICLE_ORDER),
        'beam_energy_gev': sample['beam_energy_gev'],
        'polarization_fraction': sample['polarization_fraction'],
        'polarization_angle_degrees': sample['polarization_angle_degrees'],
        'sample_sha256': sha256(sample_path),
        'pawian_revision': repository_revision(),
        'amptools_revision': manifest['reference']['revision'],
        'data': {
            'event_count': len(data),
            'event_ids': [event['event_id'] for event in data],
            'weight_sum': sum(event['weight'] for event in data),
        },
        'normalization': {
            'event_count': len(normalization),
            'event_ids': [event['event_id'] for event in normalization],
            'weight_sum': sum(event['weight'] for event in normalization),
        },
        'generated_files': [
            'data.root',
            'normalization.root',
            'prepare-normalization.cfg',
            'reference.cfg',
        ],
    }
    (output / 'input-metadata.json').write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + '\n',
        encoding='utf-8',
    )
    print(f'prepared {len(data)} data and {len(normalization)} normalization events under {output}')


def validate_result(result_path: Path, sample_path: Path) -> None:  # noqa: PLR0912, PLR0915
    result = load_json(result_path)
    sample = load_json(sample_path)
    validate_sample(sample, sample_path)
    manifest = load_json(DEFAULT_MANIFEST)
    expected_revision = manifest['reference']['revision']
    expected_events = sample['data']
    expected_ids = [event['event_id'] for event in expected_events]

    if result.get('sample_sha256') != sha256(sample_path):
        msg = f'{result_path}: sample digest does not match {sample_path}'
        raise DataReferenceError(msg)

    framework = result.get('framework', {})
    if framework.get('name') != 'AmpTools' or framework.get('revision') != expected_revision:
        msg = f'{result_path}: unexpected AmpTools revision metadata'
        raise DataReferenceError(msg)
    if tuple(result.get('particle_order', [])) != PARTICLE_ORDER:
        msg = f'{result_path}: particle order was not preserved'
        raise DataReferenceError(msg)
    if result.get('polarization_fraction') != 0.0:
        msg = f'{result_path}: reference is not unpolarized'
        raise DataReferenceError(msg)

    events: dict[str, object] = result.get('events')  # ty:ignore[invalid-assignment]
    if not isinstance(events, list) or result.get('event_count') != len(events):
        msg = f'{result_path}: inconsistent event count'
        raise DataReferenceError(msg)
    actual_ids = [event.get('event_id') for event in events]
    if actual_ids != expected_ids:
        msg = f'{result_path}: event identifiers or ordering changed'
        raise DataReferenceError(msg)
    expected_data_weight = sum(event['weight'] for event in expected_events)
    if not math.isclose(
        float(result.get('weight_sum', math.nan)),
        expected_data_weight,
        rel_tol=1e-7,
        abs_tol=1e-7,
    ):
        msg = f'{result_path}: data weight sum changed'
        raise DataReferenceError(msg)

    # AmpTools' identical-KShort permutation normalization contributes one
    # half to the unit-coefficient |Y00|^2 = 1 / (4 pi) intensity.
    expected_intensity = 1.0 / (8.0 * math.pi)
    expected_total = [
        float(sample['beam_energy_gev']) + TARGET_PROTON_MASS_GEV,
        0.0,
        0.0,
        float(sample['beam_energy_gev']),
    ]
    for event, expected_event in zip(events, expected_events, strict=True):
        event_id = event['event_id']
        if not math.isclose(
            float(event.get('weight', math.nan)),
            float(expected_event['weight']),
            rel_tol=1e-7,
            abs_tol=1e-7,
        ):
            msg = f'{result_path}: event {event_id} weight changed'
            raise DataReferenceError(msg)
        particles = event.get('particles')
        if not isinstance(particles, list) or len(particles) != len(PARTICLE_ORDER):
            msg = f'{result_path}: event {event_id} lacks ordered particles'
            raise DataReferenceError(msg)
        for name, actual_momentum in zip(PARTICLE_ORDER, particles, strict=True):
            expected_momentum = expected_event['particles'][name]
            if not isinstance(actual_momentum, list) or len(actual_momentum) != len(expected_momentum):
                msg = f'{result_path}: event {event_id} has invalid {name} momentum'
                raise DataReferenceError(msg)
            if any(
                abs(float(actual) - float(expected)) > BIG_EPS  # ty:ignore[invalid-argument-type]
                for actual, expected in zip(actual_momentum, expected_momentum, strict=True)
            ):
                msg = f'{result_path}: event {event_id} {name} momentum or ordering changed'
                raise DataReferenceError(msg)
        intensity = float(event.get('intensity', math.nan))
        if not math.isfinite(intensity) or not math.isclose(
            intensity, expected_intensity, rel_tol=5e-12, abs_tol=5e-12
        ):
            msg = (
                f'{result_path}: event {event_id} intensity {intensity} does not '
                f'match the fixed S-wave value {expected_intensity}'
            )
            raise DataReferenceError(msg)
        final_state = event.get('final_state_sum')
        if not isinstance(final_state, list) or len(final_state) != 4:  # noqa: PLR2004
            msg = f'{result_path}: event {event_id} lacks kinematic checks'
            raise DataReferenceError(msg)
        if any(abs(float(actual) - expected) > BIG_EPS for actual, expected in zip(final_state, expected_total)):  # ty:ignore[invalid-argument-type]
            msg = f'{result_path}: event {event_id} failed four-momentum checks'
            raise DataReferenceError(msg)

    normalization = result.get('normalization', {})
    expected_normalization = sample['normalization']
    normalization_source = result.get('normalization_source', {})
    for stream_name in ('generated', 'accepted'):
        stream = normalization_source.get(stream_name, {})
        source_events = stream.get('events')
        if (
            stream.get('event_count') != len(expected_normalization)
            or not isinstance(source_events, list)
            or len(source_events) != len(expected_normalization)
        ):
            msg = f'{result_path}: normalization {stream_name} event count changed'
            raise DataReferenceError(msg)
        for actual, expected in zip(source_events, expected_normalization, strict=True):
            if actual.get('event_id') != expected['event_id']:  # ty:ignore[unresolved-attribute]
                msg = f'{result_path}: normalization {stream_name} event identifiers or ordering changed'
                raise DataReferenceError(msg)
            if not math.isclose(
                float(actual.get('weight', math.nan)),  # ty:ignore[unresolved-attribute]
                float(expected['weight']),
                rel_tol=1e-7,
                abs_tol=1e-7,
            ):
                msg = f'{result_path}: normalization {stream_name} event {expected["event_id"]} weight changed'
                raise DataReferenceError(msg)
    if not normalization.get('generated_before_evaluation'):
        msg = f'{result_path}: normalization was not prepared explicitly'
        raise DataReferenceError(msg)
    if normalization.get('generated_events') != len(expected_normalization):
        msg = f'{result_path}: normalization event count changed'
        raise DataReferenceError(msg)
    expected_weight = sum(event['weight'] for event in expected_normalization)
    if not math.isclose(
        float(normalization.get('accepted_weight_sum', math.nan)),
        expected_weight,
        rel_tol=1e-12,
        abs_tol=1e-12,
    ):
        msg = f'{result_path}: normalization weight sum changed'
        raise DataReferenceError(msg)
    print(f'validated AmpTools unpolarized S-wave reference: {len(events)} events, intensity={expected_intensity:.17g}')


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)

    verify = commands.add_parser('verify-artifacts')
    verify.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    verify.add_argument('--root', type=Path, default=REPOSITORY)

    prepare_command = commands.add_parser('prepare')
    prepare_command.add_argument('--sample', type=Path, default=DEFAULT_SAMPLE)
    prepare_command.add_argument('--output', type=Path, required=True)

    digest = commands.add_parser('sample-digest')
    digest.add_argument('--sample', type=Path, default=DEFAULT_SAMPLE)

    validate = commands.add_parser('validate')
    validate.add_argument('--result', type=Path, required=True)
    validate.add_argument('--sample', type=Path, default=DEFAULT_SAMPLE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == 'verify-artifacts':
            verify_artifacts(args.manifest, args.root)
        elif args.command == 'prepare':
            prepare(args.sample, args.output)
        elif args.command == 'sample-digest':
            sample = load_json(args.sample)
            validate_sample(sample, args.sample)
            print(sha256(args.sample))
        elif args.command == 'validate':
            validate_result(args.result, args.sample)
    except (DataReferenceError, subprocess.CalledProcessError) as error:
        print(f'error: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
