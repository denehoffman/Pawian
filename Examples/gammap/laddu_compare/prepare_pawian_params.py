#!/usr/bin/env python3
"""Prepare Pawian start parameters for the laddu comparison examples."""

from __future__ import annotations

import argparse
from pathlib import Path

from comparison_parameters import PARAMETER_VALUES

ACTIVE_PRODUCTION_KEYS = {
    'minimal-j12-minus': (
        'J1/2P-1C0Tof0(1500)_proton',
        "J1/2P-1C0Tof2'(1525)_proton",
    ),
}


def should_keep_parameter(name: str, active_production_keys: tuple[str, ...]) -> bool:
    if name == 'channelType4protonKshort1Kshort2Scaling':
        return True
    if 'ToKshort1_Kshort2' in name:
        return True
    if 'Tof0(1500)_proton' in name or "Tof2'(1525)_proton" in name:
        return any(key in name for key in active_production_keys)
    return True


def zero_parameter_line(line: str) -> str:
    fields = line.split()
    if len(fields) < 2:  # noqa: PLR2004
        return line
    fields[1] = '0'
    if len(fields) >= 3:  # noqa: PLR2004
        fields[2] = '0'
    return '\t'.join(fields)


def set_parameter_line(line: str, value: float) -> str:
    fields = line.split()
    if len(fields) < 2:  # noqa: PLR2004
        return line
    fields[1] = f'{value:.17g}'
    return '\t'.join(fields)


def prepare_params(input_path: Path, output_path: Path, preset: str) -> None:
    active_keys = ACTIVE_PRODUCTION_KEYS[preset]
    output_lines = []
    for raw_line in input_path.read_text(encoding='utf-8').splitlines():
        line = raw_line.strip()
        if not line or line.startswith('#'):
            output_lines.append(raw_line)
            continue
        name = line.split()[0]
        if should_keep_parameter(name, active_keys):
            output_lines.append(set_parameter_line(raw_line, PARAMETER_VALUES.get(name, 0.0)))
        else:
            output_lines.append(zero_parameter_line(raw_line))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text('\n'.join(output_lines) + '\n', encoding='utf-8')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--preset', choices=sorted(ACTIVE_PRODUCTION_KEYS), default='minimal-j12-minus')
    args = parser.parse_args()

    prepare_params(args.input, args.out, args.preset)
    print(f'wrote {args.out} using preset {args.preset}')


if __name__ == '__main__':
    main()
