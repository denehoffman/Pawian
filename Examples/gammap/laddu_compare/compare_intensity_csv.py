#!/usr/bin/env python3
"""Compare Pawian and laddu event-by-event intensity CSV dumps."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
import matplotlib.pyplot as plt


def read_rows(path: Path) -> dict[tuple[str, int], dict[str, float | str]]:
    rows = {}
    with path.open(encoding="utf-8") as input_file:
        for row in csv.DictReader(input_file):
            key = (row["sample"], int(row["event"]))
            rows[key] = {
                "sample": row["sample"],
                "event": int(row["event"]),
                "m_ksks": float(row["m_ksks"]),
                "intensity": float(row["intensity"]),
                "weight": float(row["weight"]),
            }
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pawian", type=Path, required=True)
    parser.add_argument("--laddu", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("intensity_ratio.png"))
    parser.add_argument("--joined-csv", type=Path)
    args = parser.parse_args()

    pawian = read_rows(args.pawian)
    laddu = read_rows(args.laddu)
    keys = sorted(set(pawian).intersection(laddu))
    if not keys:
        msg = "no matching sample,event rows found"
        raise RuntimeError(msg)

    samples = []
    events = []
    masses = []
    pawian_i = []
    laddu_i = []
    ratios = []
    for key in keys:
        p_row = pawian[key]
        l_row = laddu[key]
        if l_row["intensity"] == 0.0:
            continue
        samples.append(key[0])
        events.append(key[1])
        masses.append(p_row["m_ksks"])
        p_intensity = p_row["intensity"]
        assert isinstance(p_intensity, float)
        l_intensity = l_row["intensity"]
        assert isinstance(l_intensity, float)
        pawian_i.append(p_intensity)
        laddu_i.append(l_intensity)
        ratios.append(p_intensity / l_intensity)

    masses_arr = np.asarray(masses, dtype=float)
    pawian_arr = np.asarray(pawian_i, dtype=float)
    laddu_arr = np.asarray(laddu_i, dtype=float)
    ratio_arr = np.asarray(ratios, dtype=float)
    finite = np.isfinite(ratio_arr)
    if not finite.any():
        msg = "no finite intensity ratios found"
        raise RuntimeError(msg)

    median = float(np.median(ratio_arr[finite]))
    rel_spread = float(np.std(ratio_arr[finite] / median)) if median != 0.0 else np.inf
    print(f"matched rows: {len(keys)}")
    print(f"finite ratios: {int(finite.sum())}")
    print(f"median Pawian/laddu ratio: {median:.8g}")
    print(f"relative ratio spread: {rel_spread:.8g}")

    joined_csv = args.joined_csv or args.out.with_suffix(".csv")
    joined_csv.parent.mkdir(parents=True, exist_ok=True)
    with joined_csv.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            ["sample", "event", "m_ksks", "pawian_intensity", "laddu_intensity", "ratio"]
        )
        for sample, event, mass, p_intensity, l_intensity, ratio in zip(
            samples, events, masses, pawian_i, laddu_i, ratios, strict=True
        ):
            writer.writerow([sample, event, mass, p_intensity, l_intensity, ratio])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(ncols=2, figsize=(12, 5))
    axes[0].scatter(masses_arr[finite], ratio_arr[finite] / median, s=8, alpha=0.7)
    axes[0].axhline(1.0, color="black", linewidth=1)
    axes[0].set_xlabel(r"$m(K_S K_S)$ [GeV]")
    axes[0].set_ylabel("ratio / median")

    axes[1].scatter(laddu_arr[finite], pawian_arr[finite], s=8, alpha=0.7)
    if laddu_arr[finite].size:
        x_min = float(np.min(laddu_arr[finite]))
        x_max = float(np.max(laddu_arr[finite]))
        axes[1].plot([x_min, x_max], [median * x_min, median * x_max], color="black", linewidth=1)
    axes[1].set_xlabel("laddu intensity")
    axes[1].set_ylabel("Pawian intensity")

    fig.tight_layout()
    fig.savefig(args.out, dpi=160)
    print(f"wrote {args.out}")
    print(f"wrote {joined_csv}")


if __name__ == "__main__":
    main()
