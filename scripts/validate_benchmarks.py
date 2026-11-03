#!/usr/bin/env python3

"""
Validate numerical and statistical benchmark results.

Usage:
    python3 scripts/validate_benchmarks.py tests/nist_all

    python3 scripts/validate_benchmarks.py \
        tests/nist_all \
        tests/thomas_scaling \
        tests/dense_scaling \
        tests/nist_final_20

Generated files:
    all_summary_checked.csv
    timing_variability.csv
    scaling_slopes.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics

from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# CSV UTILITIES
# ============================================================

def as_number(value):
    try:
        number = float(value)

        if math.isfinite(number):
            return number

        return None

    except (ValueError, TypeError):
        return None


def read_csv(path):
    with path.open(
        newline="",
        encoding="utf-8"
    ) as file:
        return list(csv.DictReader(file))


def write_csv(path, rows, fieldnames):
    with path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# STATISTICAL UTILITIES
# ============================================================

def quartile(sorted_values, proportion):
    position = (len(sorted_values) - 1) * proportion

    lower = math.floor(position)
    upper = math.ceil(position)

    return (
        sorted_values[lower]
        + (position - lower)
        * (sorted_values[upper] - sorted_values[lower])
    )


def least_squares_slope(points):
    """
    Estimate alpha in:

        T(n) = C * n^alpha

    By fitting:

        log(T) = alpha * log(n) + log(C)
    """

    xs = [
        math.log(n)
        for n, time in points
    ]

    ys = [
        math.log(time)
        for n, time in points
    ]

    mean_x = statistics.mean(xs)
    mean_y = statistics.mean(ys)

    denominator = sum(
        (x - mean_x) ** 2
        for x in xs
    )

    if denominator == 0:
        return math.nan

    numerator = sum(
        (x - mean_x) * (y - mean_y)
        for x, y in zip(xs, ys)
    )

    return numerator / denominator


# ============================================================
# NUMERICAL VALIDATION
# ============================================================

def analyze(root, backward_tol, forward_tol):

    summary_path = root / "all_summary.csv"

    if not summary_path.exists():
        summary_path = root / "summary.csv"

    if not summary_path.exists():
        print(
            f"[SKIP] {root}: summary.csv not found"
        )
        return

    rows = read_csv(summary_path)

    if not rows:
        print(
            f"[SKIP] {root}: empty summary file"
        )
        return

    counts = Counter()

    for row in rows:

        status = row.get("status", "")

        backward = as_number(
            row.get("backward_error_inf")
        )

        forward = as_number(
            row.get("forward_error_inf")
        )

        kappa_u = as_number(
            row.get("kappa_u")
        )

        # --------------------------------------------
        # Classify numerical quality
        # --------------------------------------------

        if status != "ok":

            quality = "factorization_failed"

        elif backward is None or forward is None:

            quality = "missing_error_metric"

        elif backward > backward_tol:

            quality = "large_backward_error"

        elif forward > forward_tol:

            if (
                kappa_u is not None
                and kappa_u >= forward_tol
            ):

                quality = (
                    "ill_conditioned_large_forward_error"
                )

            else:

                quality = "large_forward_error"

        else:

            quality = "numerically_acceptable"

        row["numerical_quality"] = quality

        # A large forward error may be explained
        # by ill-conditioning. These results should
        # be reviewed separately.
        valid_categories = (
            "numerically_acceptable",
            "ill_conditioned_large_forward_error",
        )

        row["include_in_performance_comparison"] = (
            "1" if quality in valid_categories else "0"
        )

        counts[quality] += 1

    # --------------------------------------------
    # Save numerical validation
    # --------------------------------------------

    output_path = root / "all_summary_checked.csv"

    write_csv(
        output_path,
        rows,
        list(rows[0])
    )

    print(f"\nDirectory: {root}")
    print(f"Total method executions: {len(rows)}")

    for quality, count in counts.most_common():
        print(f"  {quality}: {count}")

    print(f"Numerical validation: {output_path}")

    # ========================================================
    # TIMING VARIABILITY
    # ========================================================

    raw_path = root / "all_raw.csv"

    if not raw_path.exists():
        raw_path = root / "raw.csv"

    if raw_path.exists():

        raw_rows = read_csv(raw_path)

        groups = defaultdict(list)

        for row in raw_rows:

            seconds = as_number(
                row.get("factor_s")
            )

            measurement = "factor_s"

            # Combined factorization + solution
            # for routines such as LAPACKE dgtsv.
            if seconds is None:

                seconds = as_number(
                    row.get("total_s")
                )

                measurement = "total_s"

            if seconds is not None and seconds > 0:

                key = (
                    row.get("test", ""),
                    row.get("implementation", ""),
                    row.get("method", ""),
                    row.get("variant", ""),
                    measurement,
                )

                groups[key].append(seconds)

        statistics_rows = []

        for key, values in sorted(groups.items()):

            ordered = sorted(values)

            median = statistics.median(ordered)

            q1 = quartile(ordered, 0.25)
            q3 = quartile(ordered, 0.75)

            iqr = q3 - q1

            relative_iqr = (
                100.0 * iqr / median
            )

            result = dict(
                zip(
                    [
                        "test",
                        "implementation",
                        "method",
                        "variant",
                        "measurement",
                    ],
                    key,
                )
            )

            result.update(
                repeats=len(values),
                median_s=f"{median:.17g}",
                iqr_s=f"{iqr:.17g}",
                relative_iqr_percent=(
                    f"{relative_iqr:.9g}"
                ),
            )

            statistics_rows.append(result)

        if statistics_rows:

            output_path = (
                root / "timing_variability.csv"
            )

            write_csv(
                output_path,
                statistics_rows,
                list(statistics_rows[0])
            )

            print(
                f"Timing statistics: {output_path}"
            )

    # ========================================================
    # SCALING ANALYSIS
    # ========================================================

    scaling_groups = defaultdict(dict)

    for row in rows:

        invalid_categories = (
            "factorization_failed",
            "large_backward_error",
            "missing_error_metric",
        )

        if row["numerical_quality"] in invalid_categories:
            continue

        variant = row.get("variant", "")

        if variant not in (
            "dense_scale",
            "tridiag_scale",
        ):
            continue

        n = as_number(
            row.get("n")
        )

        if n is None or n <= 0:
            continue

        # Thomas: compare complete solution time.
        # Dense LU/Cholesky: compare factorization time.
        if variant == "tridiag_scale":

            measurement = "total_s"

        else:

            measurement = "factor_s"

        time = as_number(
            row.get(measurement)
        )

        if time is None or time <= 0:
            continue

        key = (
            variant,
            row.get("implementation", ""),
            row.get("method", ""),
            measurement,
        )

        scaling_groups[key][n] = time

    scaling_results = []

    for key, points in sorted(scaling_groups.items()):

        variant, implementation, method, measurement = key

        if len(points) < 3:
            continue

        ordered_points = sorted(
            points.items()
        )

        exponent = least_squares_slope(
            ordered_points
        )

        if variant == "tridiag_scale":
            expected_exponent = 1
        else:
            expected_exponent = 3

        scaling_results.append(
            {
                "variant": variant,
                "implementation": implementation,
                "method": method,
                "measurement": measurement,
                "n_points": len(points),
                "n_min": int(min(points)),
                "n_max": int(max(points)),
                "empirical_exponent": (
                    f"{exponent:.9g}"
                ),
                "expected_exponent": expected_exponent,
            }
        )

    if scaling_results:

        output_path = (
            root / "scaling_slopes.csv"
        )

        write_csv(
            output_path,
            scaling_results,
            list(scaling_results[0])
        )

        print(
            f"Scaling slopes: {output_path}"
        )

        for result in scaling_results:

            print(
                f"  {result['method']}: "
                f"alpha={result['empirical_exponent']} "
                f"(expected={result['expected_exponent']})"
            )


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Validate numerical and statistical "
            "benchmark results."
        )
    )

    parser.add_argument(
        "directories",
        nargs="+",
        type=Path,
        help="Directories containing benchmark CSV files",
    )

    parser.add_argument(
        "--backward-tol",
        type=float,
        default=1e-10,
        help="Backward error screening threshold",
    )

    parser.add_argument(
        "--forward-tol",
        type=float,
        default=1e-6,
        help="Forward error screening threshold",
    )

    args = parser.parse_args()

    for folder in args.directories:

        analyze(
            folder,
            args.backward_tol,
            args.forward_tol,
        )


if __name__ == "__main__":
    main()
