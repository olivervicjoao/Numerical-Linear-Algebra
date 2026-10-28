"""Gera um gráfico do padrão de esparsidade a partir de uma matriz Matrix Market.

Uso: python spy_plot_bonito.py caminho/para/matriz.mtx [-o pasta/de/saida]
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.io import mmread
from scipy.sparse import csr_matrix, issparse

BLUE = "#265DCD"
GRID = "#D7E0EA"
BG = "#F2F6FB"
TEXT = "#18283B"
MUTED = "#66788D"


def save_spy_plot(A, title: str, output_path: Path) -> None:
    """Desenha quadrículas para matrizes pequenas e spy para matrizes grandes."""

    if not issparse(A):
        A = csr_matrix(A)
    else:
        A = A.tocsr()

    A.sum_duplicates()
    A.eliminate_zeros()

    rows, cols = A.shape
    nnz = A.nnz
    density = 100.0 * nnz / (rows * cols)

    fig, ax = plt.subplots(figsize=(7.0, 7.5), dpi=170)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    fig.text(
        0.11,
        0.953,
        "Padrão de esparsidade",
        fontsize=19,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )

    fig.text(
        0.11,
        0.903,
        f"{title}   •   {rows} × {cols}   •   nnz = {nnz:,}   •   densidade = {density:.1f}%",
        fontsize=10.5,
        color=MUTED,
        ha="left",
        va="top",
    )

    if max(rows, cols) <= 30:
        # Matrizes pequenas: representação em células
        mask = (A.toarray() != 0).astype(int)

        cmap = ListedColormap([BG, BLUE])

        ax.imshow(
            mask,
            cmap=cmap,
            vmin=0,
            vmax=1,
            interpolation="nearest",
            aspect="equal",
            origin="upper",
        )

        ax.set_xticks(np.arange(cols))
        ax.set_yticks(np.arange(rows))

        ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)

        ax.grid(which="minor", color=GRID, linewidth=1.3)

        ax.tick_params(which="minor", bottom=False, left=False)

        ax.tick_params(which="major", length=0, pad=9, labelsize=11, labelcolor=MUTED)

        ax.legend(
            handles=[
                Patch(facecolor=BLUE, edgecolor="none", label="Elemento não nulo"),
                Patch(facecolor=BG, edgecolor=GRID, label="Zero"),
            ],
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.13),
            ncol=2,
            fontsize=10,
        )

    else:
        # Matrizes grandes: spy plot tradicional
        ax.spy(
            A,
            color=BLUE,
            markersize=max(0.3, min(2.0, 600 / max(rows, cols))),
            marker="s",
            aspect="equal",
        )

        ax.tick_params(length=0, pad=8, labelsize=10, labelcolor=MUTED)

        ax.grid(False)

    ax.set_xlabel("Coluna (j)", color=TEXT, fontsize=12, labelpad=10)

    ax.set_ylabel("Linha (i)", color=TEXT, fontsize=12, labelpad=10)

    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.93, bottom=0.16, top=0.86)

    # Cria o diretório de saída, se necessário
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")

    plt.close(fig)

    print(f"Figura salva em: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Gera o spy plot de uma matriz Matrix Market."
    )

    parser.add_argument("matrix", type=Path, help="Caminho para o arquivo .mtx")

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Pasta onde salvar a figura (padrão: mesma pasta do .mtx)",
    )

    args = parser.parse_args()

    filename = args.matrix

    A = mmread(filename, spmatrix=False)

    output_dir = args.output_dir if args.output_dir is not None else filename.parent

    output = output_dir / f"{filename.stem}_spy_bonito.png"

    save_spy_plot(A, filename.stem, output)


if __name__ == "__main__":
    main()
