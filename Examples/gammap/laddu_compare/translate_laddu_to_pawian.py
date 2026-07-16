#!/usr/bin/env python3
"""Convert laddu Kshort samples to Pawian default ASCII event files."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import FINAL_STATE, open_dataset


def write_pawian_ascii(input_path: Path, output_path: Path, *, write_weight: bool) -> None:
    dataset = open_dataset(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as output:
        for event in dataset.events_local:
            fields = []
            if write_weight:
                fields.append(f'{event.weight:.17g}')
            for name in FINAL_STATE:
                p4 = event.p4(name)
                fields.extend(
                    [
                        f'{p4.px:.17g}',
                        f'{p4.py:.17g}',
                        f'{p4.pz:.17g}',
                        f'{p4.e:.17g}',
                    ]
                )
            output.write(' '.join(fields))
            output.write('\n')
    print(f'wrote {output_path} with final-state order: {", ".join(FINAL_STATE)}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='laddu parquet input')
    parser.add_argument('output', type=Path, help='Pawian ASCII output')
    parser.add_argument(
        '--write-weight',
        action='store_true',
        help='prepend one event weight per row for useDataEventWeight=true',
    )
    args = parser.parse_args()
    write_pawian_ascii(args.input, args.output, write_weight=args.write_weight)


if __name__ == '__main__':
    main()
