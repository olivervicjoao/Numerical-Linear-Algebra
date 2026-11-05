#!/usr/bin/env bash

# Final experiments - Computational Linear Algebra
#
# Default execution:
#   - Thomas vs LAPACKE: 20 repetitions
#   - Dense LU and Cholesky vs LAPACKE: 5 repetitions
#   - Four selected real matrices: 20 repetitions
#
# Usage:
#   bash scripts/run_experimentos_finais.bash
#
# Individual experiments:
#   bash scripts/run_experimentos_finais.bash --only-thomas
#   bash scripts/run_experimentos_finais.bash --only-dense
#   bash scripts/run_experimentos_finais.bash --only-selected
#   bash scripts/run_experimentos_finais.bash --only-nist

set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# ============================================================
# THREAD CONFIGURATION
# ============================================================

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ============================================================
# BENCHMARK CONFIGURATION
# ============================================================

THOMAS_REPEATS=20
DENSE_REPEATS=5
SELECTED_REPEATS=20
NIST_REPEATS=5

RUN_THOMAS=1
RUN_DENSE=1
RUN_SELECTED=1
RUN_NIST=0
VALIDATE=1

usage() {
    cat <<'USAGE'
Usage:
  bash scripts/run_experimentos_finais.bash [options]

Options:
  --only-thomas      Run only Thomas scaling benchmarks
  --only-dense       Run only dense LU and Cholesky scaling
  --only-selected    Run only four selected real matrices
  --only-nist        Run all 19 NIST matrices
  --repeat-nist      Include all 19 NIST matrices
  --skip-validation  Skip CSV validation
  --help             Show this message

Default execution does not repeat the 19 NIST matrices.

Output directories:
  tests/thomas_scaling/
  tests/dense_scaling/
  tests/nist_final_20/
  tests/nist_repeticao_5/ (optional)
USAGE
}

# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

ONLY_USED=0

while (($#)); do
    case "$1" in
        --repeat-nist)
            RUN_NIST=1
            ;;

        --only-nist|--only-thomas|--only-dense|--only-selected)
            if (( ONLY_USED )); then
                echo "Error: use only one --only-* option." >&2
                exit 2
            fi

            ONLY_USED=1

            RUN_THOMAS=0
            RUN_DENSE=0
            RUN_SELECTED=0
            RUN_NIST=0

            case "$1" in
                --only-nist)
                    RUN_NIST=1
                    ;;
                --only-thomas)
                    RUN_THOMAS=1
                    ;;
                --only-dense)
                    RUN_DENSE=1
                    ;;
                --only-selected)
                    RUN_SELECTED=1
                    ;;
            esac
            ;;

        --skip-validation)
            VALIDATE=0
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac

    shift
done

# ============================================================
# CHECK REQUIRED FILES
# ============================================================

if [[ ! -x ./nla ]]; then
    echo "Error: executable ./nla not found." >&2
    echo "Compile the project before running benchmarks." >&2
    exit 1
fi

if [[ ! -f scripts/run_all_tests.sh ]]; then
    echo "Error: scripts/run_all_tests.sh not found." >&2
    exit 1
fi

# ============================================================
# SELECTED MATRICES
# ============================================================

SELECTED_MATRICES=(
    matrix/SPD/bcsstk11.mtx
    matrix/SPD/bcsstk14.mtx
    matrix/NS/sherman4.mtx
    matrix/NS/orsreg_1.mtx
)

# ============================================================
# COMPLETE NIST DATASET
# ============================================================

NIST_MATRICES=(
    matrix/SPD/bcsstk09.mtx
    matrix/SPD/bcsstk11.mtx
    matrix/SPD/bcsstk12.mtx
    matrix/SPD/bcsstk14.mtx
    matrix/SPD/bcsstk15.mtx
    matrix/SPD/bcsstk17.mtx
    matrix/SPD/bcsstk18.mtx
    matrix/NS/olm500.mtx
    matrix/NS/fidap022.mtx
    matrix/NS/olm1000.mtx
    matrix/NS/sherman4.mtx
    matrix/NS/e30r0100.mtx
    matrix/NS/utm5940.mtx
    matrix/NS/lns_3937.mtx
    matrix/NS/fidap007.mtx
    matrix/NS/orsreg_1.mtx
    matrix/NS/fidap015.mtx
    matrix/NS/e20r0100.mtx
    matrix/NS/fidap032.mtx
)

# ============================================================
# MATRIX VALIDATION
# ============================================================

check_matrices() {
    local file
    local missing=0

    for file in "$@"; do
        if [[ ! -f "$file" ]]; then
            echo "Missing matrix: $file" >&2
            missing=1
        fi
    done

    if (( missing )); then
        echo "Fix missing matrix paths before running." >&2
        exit 1
    fi
}

if (( RUN_SELECTED )); then
    check_matrices "${SELECTED_MATRICES[@]}"
fi

if (( RUN_NIST )); then
    check_matrices "${NIST_MATRICES[@]}"
fi

# ============================================================
# INITIALIZATION
# ============================================================

mkdir -p tests

echo "============================================================"
echo " COMPUTATIONAL LINEAR ALGEBRA - FINAL EXPERIMENTS"
echo "============================================================"
echo "Start: $(date -Is)"
echo "Project: $ROOT"
echo "OpenBLAS threads: $OPENBLAS_NUM_THREADS"
echo "OpenMP threads: $OMP_NUM_THREADS"
echo "Thomas repetitions: $THOMAS_REPEATS"
echo "Dense repetitions: $DENSE_REPEATS"
echo "Selected matrix repetitions: $SELECTED_REPEATS"
echo "============================================================"

# ============================================================
# SAVE COMPUTATIONAL ENVIRONMENT
# ============================================================

{
    echo "=== DATE ==="
    date -Is

    echo "=== SYSTEM ==="
    uname -a

    echo "=== CPU ==="
    lscpu 2>/dev/null || true

    echo "=== MEMORY ==="
    free -h 2>/dev/null || true

    echo "=== COMPILER ==="
    g++ --version 2>/dev/null || true

    echo "=== LINKED LIBRARIES ==="
    ldd ./nla 2>/dev/null || true

    echo "=== MAKEFILE ==="
    if [[ -f Makefile ]]; then
        cat Makefile
    fi

    echo "=== THREAD VARIABLES ==="
    env | grep -E \
        '^(OPENBLAS|OMP|MKL|BLIS|VECLIB|NUMEXPR)_.*THREADS=' \
        || true

} > tests/ambiente_computacional.txt

# ============================================================
# BENCHMARK RUNNER
# ============================================================

run_section() {
    local label="$1"
    shift

    echo
    echo "==================== $label ===================="

    "$@"

    echo "==================== [OK] $label ===================="
}

# ============================================================
# EXPERIMENT 1 - THOMAS SCALING
# ============================================================

if (( RUN_THOMAS )); then
    run_section "THOMAS VS LAPACKE" \
        bash scripts/run_all_tests.sh \
        --output tests/thomas_scaling \
        --repeats "$THOMAS_REPEATS" \
        --thomas-scale
fi

# ============================================================
# EXPERIMENT 2 - DENSE LU / CHOLESKY SCALING
# ============================================================

if (( RUN_DENSE )); then
    run_section "DENSE LU AND CHOLESKY VS LAPACKE" \
        bash scripts/run_all_tests.sh \
        --output tests/dense_scaling \
        --repeats "$DENSE_REPEATS" \
        --dense-scale
fi

# ============================================================
# EXPERIMENT 3 - SELECTED REAL MATRICES
# ============================================================

if (( RUN_SELECTED )); then
    run_section "SELECTED MATRICES - 20 REPETITIONS" \
        bash scripts/run_all_tests.sh \
        --output tests/nist_final_20 \
        --repeats "$SELECTED_REPEATS" \
        --spy \
        "${SELECTED_MATRICES[@]}"
fi

# ============================================================
# OPTIONAL - COMPLETE NIST DATASET
# ============================================================

if (( RUN_NIST )); then
    run_section "COMPLETE NIST DATASET" \
        bash scripts/run_all_tests.sh \
        --output tests/nist_repeticao_5 \
        --repeats "$NIST_REPEATS" \
        --spy \
        "${NIST_MATRICES[@]}"
fi

# ============================================================
# NUMERICAL AND STATISTICAL VALIDATION
# ============================================================

if (( VALIDATE )); then
    if [[ -f scripts/validate_benchmarks.py ]]; then

        CHECK_DIRS=()

        for dir in \
            tests/nist_all \
            tests/thomas_scaling \
            tests/dense_scaling \
            tests/nist_final_20 \
            tests/nist_repeticao_5
        do
            if [[ -f "$dir/all_summary.csv" ]]; then
                CHECK_DIRS+=("$dir")
            fi
        done

        if ((${#CHECK_DIRS[@]})); then
            echo
            echo "==================== CSV VALIDATION ===================="

            python3 scripts/validate_benchmarks.py \
                "${CHECK_DIRS[@]}" || {
                echo "Warning: CSV validation failed." >&2
                echo "Benchmark results were preserved." >&2
            }
        fi

    else
        echo "Warning: validate_benchmarks.py not found." >&2
        echo "Skipping additional validation." >&2
    fi
fi

# ============================================================
# FINAL REPORT
# ============================================================

echo
echo "============================================================"
echo "FINAL EXPERIMENTS COMPLETED"
echo "End: $(date -Is)"
echo "============================================================"

for dir in \
    tests/thomas_scaling \
    tests/dense_scaling \
    tests/nist_final_20 \
    tests/nist_repeticao_5
do
    if [[ -f "$dir/all_summary.csv" ]]; then
        echo "Results: $dir/all_summary.csv"
    fi
done

echo "Environment: tests/ambiente_computacional.txt"
echo "============================================================"
