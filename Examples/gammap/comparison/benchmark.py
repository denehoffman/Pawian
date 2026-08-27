#!/usr/bin/env python3
"""Run comparison fits quietly and report their wall-clock timings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def run_fit(args: argparse.Namespace) -> int:
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        raise ValueError('a fit command is required after --')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with args.log.open('w', encoding='utf-8') as log:
        result = subprocess.run(
            command,
            cwd=args.cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    elapsed = time.perf_counter() - started
    args.output.write_text(
        json.dumps(
            {
                'name': args.name,
                'seconds': elapsed,
                'status': result.returncode,
                'log': str(args.log.resolve()),
            },
            indent=2,
        )
        + '\n',
        encoding='utf-8',
    )
    if result.returncode != 0:
        sys.stderr.write(args.log.read_text(encoding='utf-8', errors='replace'))
    return result.returncode


def report(args: argparse.Namespace) -> int:
    timings = [json.loads(path.read_text(encoding='utf-8')) for path in args.timings]
    failures = [timing for timing in timings if timing['status'] != 0]
    if failures:
        for timing in failures:
            print(
                f"{timing['name']} failed; see {timing['log']}",
                file=sys.stderr,
            )
        return 1

    fastest = min(timing['seconds'] for timing in timings)
    width = max(len('Framework'), *(len(timing['name']) for timing in timings))
    print('\nTen-bin fit wall times')
    print(f"  {'Framework':<{width}}  Wall time  Relative")
    for timing in timings:
        relative = timing['seconds'] / fastest if fastest > 0.0 else 1.0
        print(
            f"  {timing['name']:<{width}}  {timing['seconds']:8.2f} s"
            f'  {relative:7.2f}x'
        )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='action', required=True)

    run = subparsers.add_parser('run')
    run.add_argument('--name', required=True)
    run.add_argument('--output', type=Path, required=True)
    run.add_argument('--log', type=Path, required=True)
    run.add_argument('--cwd', type=Path)
    run.add_argument('command', nargs=argparse.REMAINDER)
    run.set_defaults(handler=run_fit)

    summary = subparsers.add_parser('report')
    summary.add_argument('timings', nargs='+', type=Path)
    summary.set_defaults(handler=report)

    args = parser.parse_args()
    return args.handler(args)


if __name__ == '__main__':
    raise SystemExit(main())
