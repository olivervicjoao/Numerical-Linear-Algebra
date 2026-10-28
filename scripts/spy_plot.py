import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from scipy.io import mmread
from scipy.sparse import csr_array, issparse

BLUE = "#2563EB"
BACKGROUND = "#F8FAFC"
GRID = "#E2E8F0"
TEXT = "#172B4D"
MUTED = "#64748B"


def generate_spy_plot(
    mtx_path: Path,
    output_dir: Path,
    dpi: int = 300,
) -> None:

    # ---------------------------------------
    # 1. Leitura da matriz
    # ---------------------------------------
    A = mmread(mtx_path, spmatrix=False)

    A = A.tocsr() if issparse(A) else csr_array(A)

    A.sum_duplicates()
    A.eliminate_zeros()

    rows, cols = A.shape
    nnz = A.nnz
    density = 100.0 * nnz / (rows * cols)

    # ---------------------------------------
    # 2. Preparação da figura
    # ---------------------------------------
    fig, ax = plt.subplots(
        figsize=(8.3, 8.9),
        layout="constrained",
    )

    fig.patch.set_facecolor("white")
    ax.set_facecolor(BACKGROUND)

    fig.suptitle(
        "Padrão de esparsidade",
        fontsize=20,
        fontweight="bold",
        color=TEXT,
    )

    ax.set_title(
        f"{mtx_path.stem}  |  {rows} × {cols}  |  "
        f"nnz = {nnz:,}  |  densidade = {density:.3f}%",
        fontsize=10.5,
        color=MUTED,
        pad=14,
    )

    # ---------------------------------------
    # 3. Matrizes pequenas
    # ---------------------------------------
    if max(rows, cols) <= 30:
        pattern = (A.toarray() != 0).astype(int)

        ax.imshow(
            pattern,
            cmap=ListedColormap([BACKGROUND, BLUE]),
            vmin=0,
            vmax=1,
            interpolation="nearest",
            origin="upper",
            aspect="equal",
        )

        step = 1 if max(rows, cols) <= 12 else 5

        ax.set_xticks(np.arange(0, cols, step))
        ax.set_yticks(np.arange(0, rows, step))

        ax.set_xticks(
            np.arange(-0.5, cols, 1),
            minor=True,
        )

        ax.set_yticks(
            np.arange(-0.5, rows, 1),
            minor=True,
        )

        ax.grid(
            which="minor",
            color=GRID,
            linewidth=0.9,
        )

        ax.tick_params(
            which="minor",
            length=0,
        )

    # ---------------------------------------
    # 4. Matrizes grandes
    # ---------------------------------------
    else:
        rr, cc = A.nonzero()

        dimension = max(rows, cols)

        if dimension <= 1000:
            marker_size = 1.4
        elif dimension <= 3000:
            marker_size = 0.7
        else:
            marker_size = 0.28

        ax.scatter(
            cc,
            rr,
            s=marker_size,
            marker="s",
            color=BLUE,
            linewidths=0,
            rasterized=True,
        )

        ax.set_xticks(
            np.linspace(
                0,
                cols - 1,
                num=6,
                dtype=int,
            )
        )

        ax.set_yticks(
            np.linspace(
                0,
                rows - 1,
                num=6,
                dtype=int,
            )
        )

    # ---------------------------------------
    # 5. Configuração dos eixos
    # ---------------------------------------
    ax.set_xlim(-0.5, cols - 0.5)
    ax.set_ylim(rows - 0.5, -0.5)

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    ax.tick_params(
        axis="both",
        which="major",
        top=False,
        labeltop=False,
        bottom=True,
        labelbottom=True,
        length=0,
        pad=7,
        labelsize=10,
        labelcolor=MUTED,
    )

    ax.set_xlabel(
        "Coluna (j)",
        fontsize=12,
        color=TEXT,
        labelpad=12,
    )

    ax.set_ylabel(
        "Linha (i)",
        fontsize=12,
        color=TEXT,
        labelpad=12,
    )

    for spine in ax.spines.values():
        spine.set_color(GRID)
        spine.set_linewidth(0.8)

    # ---------------------------------------
    # 6. Salvamento em PNG e PDF
    # ---------------------------------------
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    base = output_dir / f"{mtx_path.stem}_spy_bonito"

    for extension in ("png", "pdf"):
        output = base.with_suffix(f".{extension}")

        fig.savefig(
            output,
            dpi=dpi,
            bbox_inches="tight",
            pad_inches=0.20,
            facecolor="white",
        )

        print(f"Figura salva em: {output}")

    plt.close(fig)


def main():

    parser = argparse.ArgumentParser(description="Gera spy plots em PNG e PDF.")

    parser.add_argument(
        "matrix",
        type=Path,
        help="Arquivo Matrix Market (.mtx)",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Diretório de saída",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolução da imagem PNG",
    )

    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir is not None else args.matrix.parent

    generate_spy_plot(
        args.matrix,
        output_dir,
        args.dpi,
    )


if __name__ == "__main__":
    main()
