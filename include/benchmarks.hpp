#ifndef BENCHMARKS_HPP
#define BENCHMARKS_HPP

#include "defines.hpp"
#include "graph.hpp"
#include "lapack_tests.hpp"
#include "lss.hpp"
#include "matrix.hpp"
#include "rcm.hpp"

namespace benchmarks
{
struct Options
{
    std::size_t repeats   = 5;
    bool        spy_plots = false;
};

// Parts 1-4: classification, RCM, five custom solvers, LAPACK and CSV.
void run_matrix(const std::filesystem::path& matrix_file, const std::filesystem::path& output_dir, const Options& options);

// Required scale studies: Thomas vs LAPACK dgtsv, and dense LU/Cholesky.
void run_thomas_scale(std::size_t n, const std::filesystem::path& output_dir, const Options& options);

void run_dense_scale(std::size_t n, const std::filesystem::path& output_dir, const Options& options);
}  // namespace benchmarks

#endif
