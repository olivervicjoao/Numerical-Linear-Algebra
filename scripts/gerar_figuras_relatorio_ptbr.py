#!/usr/bin/env python3
"""Regenera as cinco figuras do relatório a partir dos CSVs e matrizes.

Uso na raiz do repositório:
    python3 scripts/gerar_figuras_relatorio_ptbr.py --root . --output figuras_lista_1

Os experimentos NÃO são executados novamente. Requer: numpy, scipy, matplotlib.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread
from scipy.sparse import issparse

PLOTS = {
    "thomas": ("Thomas (implementação própria)", "LAPACKE (dgtsv)"),
    "dense": ("LU (implementação própria)", "LU (LAPACKE)",
              "Cholesky (implementação própria)", "Cholesky (LAPACKE)"),
}


def rows_from(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV não encontrado: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write(fig: plt.Figure, directory: Path, stem: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        dest = directory / (stem + "." + ext)
        fig.savefig(dest, bbox_inches="tight", dpi=240)
        print(f"Gerado: {dest}")
    plt.close(fig)


def scaling(root: Path, output: Path, flavor: str) -> None:
    config = {
        "thomas": {
            "path": "tests/thomas_scaling/all_summary.csv",
            "column": "total_s", "variant": "tridiag_scale",
            "title": "Escalabilidade do algoritmo de Thomas",
            "ylabel": "Tempo total (s)",
            "stem": "escalabilidade_thomas_cluster",
            "names": [
                ("custom_thomas", "Thomas (implementação própria)"),
                ("lapack_dgtsv", "LAPACKE (dgtsv)"),
            ],
        },
        "dense": {
            "path": "tests/dense_scaling/all_summary.csv",
            "column": "factor_s", "variant": "dense_scale",
            "title": "Escalabilidade das fatorações LU e Cholesky",
            "ylabel": "Tempo de fatoração (s)",
            "stem": "escalabilidade_lu_cholesky_cluster",
            "names": [
                ("custom_lu_dense", "LU (implementação própria)"),
                ("lapack_dgetrf", "LU (LAPACKE)"),
                ("custom_cholesky_dense", "Cholesky (implementação própria)"),
                ("lapack_dpotrf", "Cholesky (LAPACKE)"),
            ],
        },
    }[flavor]
    data = defaultdict(list)
    for row in rows_from(root / config["path"]):
        if row.get("status") != "ok" or row.get("variant") != config["variant"]:
            continue
        try:
            x, y = int(row["n"]), float(row[config["column"]])
            if x > 0 and math.isfinite(y) and y > 0:
                data[row["method"]].append((x, y))
        except (TypeError, ValueError, KeyError):
            continue
    fig, ax = plt.subplots(figsize=(8.7, 5.3))
    for method, label in config["names"]:
        if method not in data:
            raise ValueError(f"Falta {method} no CSV {config['path']}")
        points = sorted(data[method])
        ax.loglog([p[0] for p in points], [p[1] for p in points], "o-", lw=2, ms=6, label=label)
    ax.set_xlabel("Dimensão da matriz ($n$)")
    ax.set_ylabel(config["ylabel"])
    ax.set_title(config["title"], pad=13)
    ax.grid(which="both", alpha=.24)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    write(fig, output, config["stem"])


def dense_advantage(root: Path, output: Path) -> None:
    collection = defaultdict(dict)
    for row in rows_from(root / "tests/dense_scaling/all_summary.csv"):
        if row.get("variant") != "dense_scale" or row.get("status") != "ok":
            continue
        collection[int(row["n"])][row["method"]] = float(row["factor_s"])
    nx = sorted(collection)
    series = [
        ("LU densa", "custom_lu_dense", "lapack_dgetrf"),
        ("Cholesky densa", "custom_cholesky_dense", "lapack_dpotrf"),
    ]
    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    for label, own, lap in series:
        try:
            ys = [collection[n][own] / collection[n][lap] for n in nx]
        except KeyError as exc:
            raise ValueError("Dados de LU/Cholesky ausentes para o gráfico de aceleração") from exc
        ax.plot(nx, ys, "o-", lw=2, ms=6, label=label)
    ax.set_xscale("log", base=2)
    ax.set_xticks(nx)
    ax.set_xticklabels([str(n) for n in nx])
    ax.set_xlabel("Dimensão da matriz ($n$)")
    ax.set_ylabel("Razão de tempos (própria / LAPACKE)")
    ax.set_title("Vantagem da LAPACKE nas fatorações densas", pad=13)
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    write(fig, output, "desempenho_fatoracoes_densas_ptbr")


def nonzero_positions(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Matriz não encontrada: {path}")
    mat = mmread(path)
    if issparse(mat):
        mat = mat.tocoo()
        mat.sum_duplicates()
        mat.eliminate_zeros()
        return mat.shape, mat.row, mat.col
    r, c = np.nonzero(mat)
    return mat.shape, r, c


def paired_spy(root: Path, output: Path, names: list[str], stem: str) -> None:
    paths = {
        "bcsstk11": ("matrix/SPD/bcsstk11.mtx", "matrix_1_bcsstk11"),
        "bcsstk14": ("matrix/SPD/bcsstk14.mtx", "matrix_2_bcsstk14"),
        "sherman4": ("matrix/NS/sherman4.mtx", "matrix_3_sherman4"),
        "orsreg_1": ("matrix/NS/orsreg_1.mtx", "matrix_4_orsreg_1"),
    }
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 9.8))
    for i, name in enumerate(names):
        original, dir_name = paths[name]
        p0 = root / original
        p1 = root / "tests/nist_final_20" / dir_name / "matrix_rcm.mtx"
        for j, (path, prefix) in enumerate(((p0, "Original"), (p1, "Após RCM"))):
            ax = axes[i, j]
            shape, r, c = nonzero_positions(path)
            ax.scatter(c, r, s=.22, marker="s", alpha=.94, rasterized=True)
            ax.set_xlim(-.5, shape[1]-.5)
            ax.set_ylim(shape[0]-.5, -.5)
            ax.set_aspect("equal")
            ax.set_title(f"{name} — {prefix}", fontsize=12)
            ax.set_xlabel("Coluna (j)", fontsize=10)
            ax.set_ylabel("Linha (i)", fontsize=10)
            ax.tick_params(labelsize=8)
    fig.suptitle("Padrões de esparsidade: original e após RCM", fontsize=15, y=.99)
    fig.tight_layout(rect=[0, 0, 1, .967])
    write(fig, output, stem)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."), help="raiz do repositório")
    parser.add_argument("--output", type=Path, default=Path("figuras_lista_1"), help="pasta de saída")
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    scaling(root, output, "thomas")
    scaling(root, output, "dense")
    dense_advantage(root, output)
    paired_spy(root, output, ["bcsstk11", "bcsstk14"], "rcm_estruturais_ptbr")
    paired_spy(root, output, ["sherman4", "orsreg_1"], "rcm_reservatorios_ptbr")


if __name__ == "__main__":
    main()
