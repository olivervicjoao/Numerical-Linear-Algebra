#ifndef LAPACK_TESTS_HPP
#define LAPACK_TESTS_HPP

#include "defines.hpp"
#include "matrix.hpp"

class matrix;

// All LAPACKE calls and LAPACK-specific layouts live in lapack_tests.cc.
namespace lapack_tests
{
struct Result
{
    // For dgtsv, factorization/solve cannot be separated: factor_s and
    // solve_s are NaN and total_s measures the combined LAPACK call.
    double              factor_s;
    double              solve_s;
    double              total_s;
    double              growth_factor;
    bool                combined;
    std::vector<double> x;
};

// Dense input is column-major, of size n*n.
Result dense_lu(const std::vector<double>& a, const std::vector<double>& b, std::size_t n, double max_original);

Result dense_cholesky(const std::vector<double>& a, const std::vector<double>& b, std::size_t n);

// Input matrices use the project's matrix interface; these functions
// pack them into LAPACK's own band formats outside the measured times.
Result band_lu(const matrix& a, const std::vector<double>& b, double max_original);

Result band_cholesky(const matrix& a, const std::vector<double>& b);

Result tridiagonal(const std::vector<double>& lower, const std::vector<double>& diagonal, const std::vector<double>& upper, const std::vector<double>& b);

// ||A||_1 estimated with dlange; kappa1 via dgetrf + dgecon.
double estimate_condition_number_1(const std::vector<double>& a, std::size_t n);
}  // namespace lapack_tests

#endif
