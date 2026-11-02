#!/usr/bin/env python3
"""Generate report-ready log-log plots and side-by-side spy plots from run_all_tests.sh output."""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def plot_scaling(root: Path, variant: str, filename: str, use_total: bool) -> None:
    path = root / "all_summary.csv"
    if not path.exists():
        return
    points = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            if row["variant"] != variant or row["status"] != "ok":
                continue
            time_field = "total_s" if use_total else "factor_s"
            try:
                x = int(row["n"])
                y = float(row[time_field])
                if x > 0 and math.isfinite(y) and y > 0:
                    points[row["method"]].append((x, y))
            except (KeyError, ValueError):
                continue
    if not points:
        return
    plt.figure(figsize=(8, 5))
    for method, data in sorted(points.items()):
        data.sort()
        plt.loglog([p[0] for p in data], [p[1] for p in data], "o-", label=method)
    plt.xlabel("Matrix dimension n")
    plt.ylabel("Total time (s)" if use_total else "Factorization time (s)")
    plt.title("Thomas scaling" if use_total else "Dense LU / Cholesky scaling")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(root / (filename + ".png"), dpi=180)
    plt.savefig(root / (filename + ".pdf"))
    plt.close()


def compare_spy(root: Path) -> None:
    for folder in root.iterdir():
        if not folder.is_dir() or not folder.name.startswith("matrix_"):
            continue
        after = folder / "matrix_rcm_spy_bonito.png"
        before = [p for p in folder.glob("*_spy_bonito.png")
                  if p.name != after.name]
        if not after.exists() or not before:
            continue
        before.sort()
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        for ax, path, title in zip(axes, (before[0], after), ("Original", "After RCM")):
            ax.imshow(plt.imread(path))
            ax.set_title(title)
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(folder / "spy_before_after.png", dpi=180)
        fig.savefig(folder / "spy_before_after.pdf")
        plt.close(fig)


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python3 scripts/plot_results.py <benchmark_dir>", file=sys.stderr)
        return 1
    root = Path(sys.argv[1])
    plot_scaling(root, "tridiag_scale", "thomas_scaling_loglog", use_total=True)
    plot_scaling(root, "dense_scale", "dense_scaling_loglog", use_total=False)
    compare_spy(root)
    print("Plots generated where matching data were available.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
