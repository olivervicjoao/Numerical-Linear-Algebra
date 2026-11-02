#include "../include/lapack_tests.hpp"

namespace
{
using Clock = std::chrono::steady_clock;

double elapsed(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration<double>(b - a).count();
}

lapack_int checked_int(std::size_t n)
{
    if (n > static_cast<std::size_t>(std::numeric_limits<lapack_int>::max()))
        throw std::invalid_argument("Dimension exceeds LAPACK integer range!");
    return static_cast<lapack_int>(n);
}

void check_info(const char *method, lapack_int info)
{
    if (info != 0)
        throw std::runtime_error(std::string(method) + " failed, info=" + std::to_string(info));
}

double max_upper_dense(const std::vector<double>& lu, std::size_t n)
{
    double value = 0.0;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i <= j; ++i)
            value = std::max(value, std::abs(lu[j * n + i]));
    return value;
}

double max_upper_band(const std::vector<double>& ab, std::size_t n, std::size_t kl, std::size_t ku)
{
    const std::size_t ldab  = 2 * kl + ku + 1;
    double            value = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t i_min = (j > kl + ku) ? j - kl - ku : 0;
        for (std::size_t i = i_min; i <= j; ++i)
            value = std::max(value, std::abs(ab[j * ldab + kl + ku + i - j]));
    }
    return value;
}

std::vector<double> pack_general_band(const matrix& a)
{
    const std::size_t   n = a.rows(), kl = a.bi(), ku = a.bs();
    const std::size_t   ldab = 2 * kl + ku + 1;
    std::vector<double> ab(ldab * n, 0.0);
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t i_min = j > ku ? j - ku : 0;
        const std::size_t i_max = j + std::min(kl, n - 1 - j);
        for (std::size_t i = i_min; i <= i_max; ++i)
            ab[j * ldab + kl + ku + i - j] = a.get(i, j);
    }
    return ab;
}

std::vector<double> pack_cholesky_band(const matrix& a)
{
    const std::size_t   n = a.rows(), kd = a.bi(), ldab = kd + 1;
    std::vector<double> ab(ldab * n, 0.0);
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t i_max = j + std::min(kd, n - 1 - j);
        for (std::size_t i = j; i <= i_max; ++i)
            ab[j * ldab + i - j] = a.get(i, j);
    }
    return ab;
}

lapack_tests::Result separate_result(double factor_s, double solve_s, double growth, std::vector<double> x)
{
    return { factor_s, solve_s, factor_s + solve_s, growth, false, std::move(x) };
}
}  // namespace

namespace lapack_tests
{
Result dense_lu(const std::vector<double>& a, const std::vector<double>& b, std::size_t n, double max_original)
{
    if (n == 0 || a.size() != n * n || b.size() != n)
        throw std::invalid_argument("Invalid dimensions for LAPACK dense LU!");
    const lapack_int        N    = checked_int(n);
    auto                    work = a, x = b;
    std::vector<lapack_int> ipiv(n);

    const auto t0 = Clock::now();
    check_info("dgetrf", LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, work.data(), N, ipiv.data()));
    const auto t1 = Clock::now();
    check_info("dgetrs", LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', N, 1, work.data(), N, ipiv.data(), x.data(), N));
    const auto   t2     = Clock::now();
    const double growth = max_original > 0.0 ? max_upper_dense(work, n) / max_original : std::numeric_limits<double>::quiet_NaN();
    return separate_result(elapsed(t0, t1), elapsed(t1, t2), growth, std::move(x));
}

Result dense_cholesky(const std::vector<double>& a, const std::vector<double>& b, std::size_t n)
{
    if (n == 0 || a.size() != n * n || b.size() != n)
        throw std::invalid_argument("Invalid dimensions for LAPACK dense Cholesky!");
    const lapack_int N    = checked_int(n);
    auto             work = a, x = b;

    const auto t0 = Clock::now();
    check_info("dpotrf", LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', N, work.data(), N));
    const auto t1 = Clock::now();
    check_info("dpotrs", LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', N, 1, work.data(), N, x.data(), N));
    const auto t2 = Clock::now();
    return separate_result(elapsed(t0, t1), elapsed(t1, t2), std::numeric_limits<double>::quiet_NaN(), std::move(x));
}

Result band_lu(const matrix& a, const std::vector<double>& b, double max_original)
{
    const std::size_t n = a.rows(), kl = a.bi(), ku = a.bs();
    if (n == 0 || a.cols() != n || b.size() != n)
        throw std::invalid_argument("Invalid dimensions for LAPACK band LU!");
    const lapack_int        N  = checked_int(n);
    const lapack_int        KL = checked_int(kl), KU = checked_int(ku);
    const lapack_int        LDAB = checked_int(2 * kl + ku + 1);
    auto                    work = pack_general_band(a), x = b;
    std::vector<lapack_int> ipiv(n);

    const auto t0 = Clock::now();
    check_info("dgbtrf", LAPACKE_dgbtrf(LAPACK_COL_MAJOR, N, N, KL, KU, work.data(), LDAB, ipiv.data()));
    const auto t1 = Clock::now();
    check_info("dgbtrs", LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', N, KL, KU, 1, work.data(), LDAB, ipiv.data(), x.data(), N));
    const auto   t2     = Clock::now();
    const double growth = max_original > 0.0 ? max_upper_band(work, n, kl, ku) / max_original : std::numeric_limits<double>::quiet_NaN();
    return separate_result(elapsed(t0, t1), elapsed(t1, t2), growth, std::move(x));
}

Result band_cholesky(const matrix& a, const std::vector<double>& b)
{
    const std::size_t n = a.rows(), kd = a.bi();
    if (n == 0 || a.cols() != n || b.size() != n || a.bs() != kd)
        throw std::invalid_argument("Invalid dimensions/bands for LAPACK band Cholesky!");
    const lapack_int N = checked_int(n), KD = checked_int(kd);
    const lapack_int LDAB = checked_int(kd + 1);
    auto             work = pack_cholesky_band(a), x = b;

    const auto t0 = Clock::now();
    check_info("dpbtrf", LAPACKE_dpbtrf(LAPACK_COL_MAJOR, 'L', N, KD, work.data(), LDAB));
    const auto t1 = Clock::now();
    check_info("dpbtrs", LAPACKE_dpbtrs(LAPACK_COL_MAJOR, 'L', N, KD, 1, work.data(), LDAB, x.data(), N));
    const auto t2 = Clock::now();
    return separate_result(elapsed(t0, t1), elapsed(t1, t2), std::numeric_limits<double>::quiet_NaN(), std::move(x));
}

Result tridiagonal(const std::vector<double>& lower, const std::vector<double>& diagonal, const std::vector<double>& upper, const std::vector<double>& b)
{
    const std::size_t n = diagonal.size();
    if (n == 0 || lower.size() != n - 1 || upper.size() != n - 1 || b.size() != n)
        throw std::invalid_argument("Invalid dimensions for LAPACK dgtsv!");
    const lapack_int N  = checked_int(n);
    auto             dl = lower, d = diagonal, du = upper, x = b;
    const auto       t0 = Clock::now();
    check_info("dgtsv", LAPACKE_dgtsv(LAPACK_COL_MAJOR, N, 1, dl.data(), d.data(), du.data(), x.data(), N));
    const auto   t1 = Clock::now();
    const double dt = elapsed(t0, t1);
    const double na = std::numeric_limits<double>::quiet_NaN();
    return { na, na, dt, na, true, std::move(x) };
}

double estimate_condition_number_1(const std::vector<double>& a, std::size_t n)
{
    if (n == 0 || a.size() != n * n)
        throw std::invalid_argument("Invalid dimensions for condition estimate!");
    const lapack_int        N     = checked_int(n);
    const double            anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', N, N, a.data(), N);
    auto                    work  = a;
    std::vector<lapack_int> ipiv(n);
    const lapack_int        info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, work.data(), N, ipiv.data());
    if (info < 0)
        check_info("dgetrf (condition)", info);
    if (info > 0)
        return std::numeric_limits<double>::infinity();
    double rcond = 0.0;
    check_info("dgecon", LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', N, work.data(), N, anorm, &rcond));
    return rcond > 0.0 ? 1.0 / rcond : std::numeric_limits<double>::infinity();
}
}  // namespace lapack_tests
