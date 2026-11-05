#!/usr/bin/env bash
set -Eeuo pipefail

# Examples:
# ./scripts/run_all_tests.sh --output tests/full --repeats 5 --spy matrix/SPD/A.mtx matrix/NS/B.mtx
# ./scripts/run_all_tests.sh --output tests/scaling --thomas-scale --dense-scale

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUTPUT="tests/benchmarks"
REPEATS=5
SPY=0
THOMAS_SCALE=0
DENSE_SCALE=0
MATRICES=()
THOMAS_SIZES=(1000 10000 100000 1000000)
DENSE_SIZES=(250 500 1000 2000 4000)

while (($#)); do
    case "$1" in
        --output) OUTPUT="$2"; shift 2;;
        --repeats) REPEATS="$2"; shift 2;;
        --spy) SPY=1; shift;;
        --thomas-scale) THOMAS_SCALE=1; shift;;
        --dense-scale) DENSE_SCALE=1; shift;;
        --help|-h)
            echo "Usage: $0 [--output DIR] [--repeats N>=5] [--spy] [--thomas-scale] [--dense-scale] [matrix1.mtx ...]"
            exit 0;;
        --*) echo "Unknown option: $1" >&2; exit 1;;
        *) MATRICES+=("$1"); shift;;
    esac
done

if [[ ! -x ./nla ]]; then
    echo "Executable ./nla not found; compile first with make." >&2
    exit 1
fi

mkdir -p "$OUTPUT"
OUTPUT="$(realpath -m "$OUTPUT")"

if (( REPEATS < 5 )); then
    echo "At least 5 repeats required." >&2
    exit 1
fi

# Matrix names may repeat; numeric prefix keeps the output directories unique.
index=0
for matrix_file in "${MATRICES[@]}"; do
    ((index+=1))
    if [[ ! -f "$matrix_file" ]]; then
        echo "Matrix not found: $matrix_file" >&2
        exit 1
    fi
    name="$(basename "$matrix_file" .mtx)"
    out="$OUTPUT/matrix_${index}_${name}"
    args=(--matrix "$matrix_file" --output "$out" --repeats "$REPEATS")
    (( SPY )) && args+=(--spy)
    echo "=== Matrix: $matrix_file ==="
    ./nla "${args[@]}" | tee "$out.log"
done

if (( THOMAS_SCALE )); then
    for n in "${THOMAS_SIZES[@]}"; do
        out="$OUTPUT/thomas_scale_${n}"
        echo "=== Thomas scaling: n=$n ==="
        ./nla --thomas-scale "$n" --output "$out" --repeats "$REPEATS" | tee "$out.log"
    done
fi

if (( DENSE_SCALE )); then
    for n in "${DENSE_SIZES[@]}"; do
        out="$OUTPUT/dense_scale_${n}"
        echo "=== Dense scaling: n=$n ==="
        ./nla --dense-scale "$n" --output "$out" --repeats "$REPEATS" | tee "$out.log"
    done
fi

# Aggregate all CSVs, inserting the source test directory for traceability.
python3 - "$OUTPUT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
for name in ("summary.csv", "raw.csv", "matrix_metadata.csv", "part3_3.csv", "part3_comparison.csv", "preprocessing.csv", "preprocessing_raw.csv"):
    files = sorted(p for p in root.glob("*/" + name))
    if not files:
        continue
    with (root / ("all_" + name)).open("w", newline="", encoding="utf-8") as dst:
        writer = None
        for p in files:
            with p.open(newline="", encoding="utf-8") as src:
                reader = csv.DictReader(src)
                if writer is None:
                    writer = csv.DictWriter(dst, fieldnames=["test"] + (reader.fieldnames or []))
                    writer.writeheader()
                for row in reader:
                    writer.writerow({"test": p.parent.name, **row})
print("Aggregated CSV files in", root)
PY

if [[ -f scripts/plot_results.py ]]; then
    python3 scripts/plot_results.py "$OUTPUT" || echo "Warning: could not generate report plots" >&2
fi

echo "Completed. Results: $OUTPUT"
