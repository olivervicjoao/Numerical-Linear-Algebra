#include "../include/benchmarks.hpp"

namespace
{
using Clock     = std::chrono::steady_clock;
using Scalar    = real;
const double NA = std::numeric_limits<double>::quiet_NaN();

struct Trial
{
    double              factor_s = 0.0;
    double              solve_s  = 0.0;
    double              growth   = NA;
    bool                combined = false;
    std::vector<Scalar> x;
};

struct Row
{
    std::string         method, variant, status = "ok";
    bool                combined = false;
    std::size_t         n = 0, p = 0, q = 0, repeats = 0;
    double              factor_s = NA, solve_s = NA, total_s = NA;
    double              factor_flops = NA, total_flops = NA, gflops = NA;
    double              backward_error = NA, forward_error = NA, lapack_diff = NA;
    double              growth = NA, kappa1 = NA, kappa_u = NA;
    std::vector<Scalar> x_original;
};

struct Raw
{
    std::string method, variant;
    std::size_t repeat;
    double      factor_s, solve_s, total_s;
    bool        combined;
};

Trial as_trial(lapack_tests::Result r)
{
    // For combined dgtsv store its sole elapsed time internally in solve_s;
    // add_benchmark uses combined=true to expose it only in total_s.
    return { r.combined ? 0.0 : r.factor_s, r.combined ? r.total_s : r.solve_s, r.growth_factor, r.combined, std::move(r.x) };
}

std::string implementation(const std::string& method)
{
    return method.compare(0, 7, "lapack_") == 0 ? "lapacke" : "propria";
}

double seconds(Clock::time_point start, Clock::time_point stop)
{
    return std::chrono::duration<double>(stop - start).count();
}

double median(std::vector<double> a)
{
    if (a.empty())
        return NA;
    std::sort(a.begin(), a.end());
    const std::size_t m = a.size() / 2;
    return (a.size() % 2) ? a[m] : (a[m - 1] + a[m]) / 2.0;
}

void write_header(std::ofstream& f)
{
    f << "implementation,method,variant,status,n,p,q,repeats,combined,factor_s,solve_s,total_s,"
         "factor_flops_est,total_flops_est,factor_gflops,backward_error_inf,"
         "forward_error_inf,diff_from_lapack_inf,lu_growth,kappa1_est,kappa_u\n";
    f << std::setprecision(17);
}

void write_row(std::ofstream& f, const Row& r)
{
    f << implementation(r.method) << ',' << r.method << ',' << r.variant << ',' << r.status << ',' << r.n << ',' << r.p << ',' << r.q << ',' << r.repeats << ','
      << (r.combined ? 1 : 0) << ',' << r.factor_s << ',' << r.solve_s << ',' << r.total_s << ',' << r.factor_flops << ',' << r.total_flops << ',' << r.gflops
      << ',' << r.backward_error << ',' << r.forward_error << ',' << r.lapack_diff << ',' << r.growth << ',' << r.kappa1 << ',' << r.kappa_u << '\n';
}

std::ofstream new_csv(const std::filesystem::path& path)
{
    std::ofstream file(path);
    if (!file)
        throw std::runtime_error("Cannot create file: " + path.string());
    return file;
}

std::vector<Scalar> matrix_times_ones(const matrix& A)
{
    const std::size_t   n = A.rows();
    std::vector<Scalar> b(n, 0.0);
    const std::size_t   p = A.bi(), q = A.bs();
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imin = j > q ? j - q : 0;
        const std::size_t imax = j + std::min(p, n - 1 - j);
        for (std::size_t i = imin; i <= imax; ++i)
            b[i] += A.get(static_cast<unsigned>(i), static_cast<unsigned>(j));
    }
    return b;
}

// Backward error per assignment: ||b-A*x||inf / (||A||inf * ||x||inf).
double backward_error(const matrix& A, const std::vector<Scalar>& x, const std::vector<Scalar>& b)
{
    const std::size_t   n = A.rows();
    std::vector<double> r(n, 0.0), rowsum(n, 0.0);
    const std::size_t   p = A.bi(), q = A.bs();
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imin = j > q ? j - q : 0;
        const std::size_t imax = j + std::min(p, n - 1 - j);
        for (std::size_t i = imin; i <= imax; ++i)
        {
            const double a = A.get(static_cast<unsigned>(i), static_cast<unsigned>(j));
            r[i] += a * x[j];
            rowsum[i] += std::abs(a);
        }
    }
    double nr = 0.0, na = 0.0, nx = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        nr = std::max(nr, std::abs(b[i] - r[i]));
        na = std::max(na, rowsum[i]);
        nx = std::max(nx, std::abs(x[i]));
    }
    return (na * nx > 0.0) ? nr / (na * nx) : NA;
}

double forward_error(const std::vector<Scalar>& x)
{
    double err = 0.0;
    for (double value : x)
        err = std::max(err, std::abs(value - 1.0));
    return err;
}

double difference_inf(const std::vector<Scalar>& a, const std::vector<Scalar>& b)
{
    if (a.size() != b.size())
        return NA;
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        diff = std::max(diff, std::abs(a[i] - b[i]));
    return diff;
}

double max_abs_matrix(const matrix& A)
{
    const std::size_t n = A.rows(), p = A.bi(), q = A.bs();
    double            result = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imin = j > q ? j - q : 0;
        const std::size_t imax = j + std::min(p, n - 1 - j);
        for (std::size_t i = imin; i <= imax; ++i)
            result = std::max(result, std::abs(A.get(static_cast<unsigned>(i), static_cast<unsigned>(j))));
    }
    return result;
}

double growth_matrix(const matrix& LU, double max_original, std::size_t ku)
{
    if (max_original == 0.0)
        return NA;
    double            max_upper = 0.0;
    const std::size_t n         = LU.rows();
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imin = j > ku ? j - ku : 0;
        for (std::size_t i = imin; i <= j; ++i)
            max_upper = std::max(max_upper, std::abs(LU.get(static_cast<unsigned>(i), static_cast<unsigned>(j))));
    }
    return max_upper / max_original;
}

double lu_flops(std::size_t n, std::size_t p, std::size_t q)
{
    long double count = 0.0;
    for (std::size_t k = 0; k < n; ++k)
    {
        const long double l = std::min(p, n - 1 - k);
        const long double u = std::min(q, n - 1 - k);
        count += l + 2 * l * u;
    }
    return static_cast<double>(count);
}

double chol_flops(std::size_t n, std::size_t p)
{
    long double count = 0.0;
    for (std::size_t k = 0; k < n; ++k)
    {
        const long double m = std::min(p, n - 1 - k);
        count += m * m + 2 * m;
    }
    return static_cast<double>(count);
}

std::vector<Scalar> dense_data(const matrix& A)
{
    const std::size_t   n = A.rows();
    std::vector<Scalar> data(n * n);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
            data[j * n + i] = A.get(static_cast<unsigned>(i), static_cast<unsigned>(j));
    return data;
}

bool symmetric(const matrix& A)
{
    if (A.rows() != A.cols() || A.bi() != A.bs())
        return false;
    const std::size_t n = A.rows(), p = A.bi();
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imax = j + std::min(p, n - 1 - j);
        for (std::size_t i = j + 1; i <= imax; ++i)
        {
            const double a   = A.get(static_cast<unsigned>(i), static_cast<unsigned>(j));
            const double b   = A.get(static_cast<unsigned>(j), static_cast<unsigned>(i));
            const double tol = 1e-13 * std::max({ 1.0, std::abs(a), std::abs(b) });
            if (std::abs(a - b) > tol)
                return false;
        }
    }
    return true;
}

std::size_t count_nonzero(const matrix& A)
{
    const std::size_t n = A.rows(), p = A.bi(), q = A.bs();
    std::size_t       count = 0;
    for (std::size_t j = 0; j < n; ++j)
    {
        const std::size_t imin = j > q ? j - q : 0;
        const std::size_t imax = j + std::min(p, n - 1 - j);
        for (std::size_t i = imin; i <= imax; ++i)
            if (A.get(static_cast<unsigned>(i), static_cast<unsigned>(j)) != 0.0)
                ++count;
    }
    return count;
}

void save_rows(const std::filesystem::path& dir, const std::vector<Row>& rows, const std::vector<Raw>& raw)
{
    auto file = new_csv(dir / "summary.csv");
    write_header(file);
    for (const auto& row : rows)
        write_row(file, row);

    auto detail = new_csv(dir / "raw.csv");
    detail << "implementation,method,variant,repeat,combined,factor_s,solve_s,total_s\n" << std::setprecision(17);
    for (const auto& r : raw)
        detail << implementation(r.method) << ',' << r.method << ',' << r.variant << ',' << r.repeat << ',' << (r.combined ? 1 : 0) << ',' << r.factor_s << ','
               << r.solve_s << ',' << r.total_s << '\n';
}

std::string csv_status(std::string text)
{
    for (char& c : text)
        if (c == ',' || c == '\n' || c == '\r')
            c = ';';
    return text;
}

// Each measured run owns fresh inputs; setup / copying isn't timed.
template <typename MakeTrial>
void add_benchmark(std::vector<Row>& rows, std::vector<Raw>& raw, const std::string& name, const std::string& variant, std::size_t n, std::size_t p,
                   std::size_t q, double factor_flops, double total_flops, const matrix *A_for_error, const std::vector<Scalar> *b_for_error,
                   const std::vector<std::size_t> *permutation, std::size_t repeats, MakeTrial make_trial)
{
    Row result;
    result.method       = name;
    result.variant      = variant;
    result.n            = n;
    result.p            = p;
    result.q            = q;
    result.repeats      = repeats;
    result.factor_flops = factor_flops;
    result.total_flops  = total_flops;

    std::vector<double> factor_times, solve_times, total_times;
    try
    {
        for (std::size_t rep = 0; rep < repeats; ++rep)
        {
            Trial t = make_trial();
            if (!t.combined)
            {
                factor_times.push_back(t.factor_s);
                solve_times.push_back(t.solve_s);
            }
            else
                result.combined = true;
            total_times.push_back(t.factor_s + t.solve_s);
            raw.push_back({ name, variant, rep, t.combined ? NA : t.factor_s, t.combined ? NA : t.solve_s, t.factor_s + t.solve_s, t.combined });
            if (rep == 0)
            {
                result.growth = t.growth;
                if (permutation)
                {
                    if (permutation->size() != n || t.x.size() != n)
                        throw std::runtime_error("Invalid reordered solution dimensions");
                    result.x_original.assign(n, 0.0);
                    for (std::size_t i = 0; i < n; ++i)
                        result.x_original[permutation->at(i)] = t.x[i];
                }
                else
                    result.x_original = std::move(t.x);
            }
        }
        result.factor_s         = median(factor_times);
        result.solve_s          = median(solve_times);
        result.total_s          = median(total_times);
        const double timed_s    = factor_flops > 0.0 ? result.factor_s : result.total_s;
        const double flop_count = factor_flops > 0.0 ? factor_flops : total_flops;
        result.gflops           = timed_s > 0.0 ? flop_count / (timed_s * 1.0e9) : NA;
        if (A_for_error && b_for_error)
        {
            result.backward_error = backward_error(*A_for_error, result.x_original, *b_for_error);
            result.forward_error  = forward_error(result.x_original);
        }
    }
    catch (const std::exception& e)
    {
        result.status = "failed: " + csv_status(e.what());
        result.x_original.clear();
        std::cerr << "[skip] " << variant << '/' << name << ": " << e.what() << '\n';
    }
    rows.push_back(std::move(result));
}

void compare_to_reference(std::vector<Row>& rows)
{
    const Row *lapack = nullptr;
    for (const auto& row : rows)
        if (row.method == "lapack_dgetrf" && row.variant == "original" && row.status == "ok")
        {
            lapack = &row;
            break;
        }
    if (!lapack)
        return;
    for (auto& row : rows)
        if (row.status == "ok")
            row.lapack_diff = difference_inf(row.x_original, lapack->x_original);
}

std::string quote_shell(const std::string& s)
{
    std::string r = "'";
    for (char c : s)
        r += c == '\'' ? "'\\''" : std::string(1, c);
    return r + "'";
}

void plot(const std::filesystem::path& matrix_file, const std::filesystem::path& output_dir)
{
    const std::string cmd = "python3 scripts/spy_plot.py " + quote_shell(matrix_file.string()) + " -o " + quote_shell(output_dir.string());
    if (std::system(cmd.c_str()) != 0)
        std::cerr << "Warning: failed generating spy plot for " << matrix_file << '\n';
}

// Dense Cholesky and LU benchmarks on an in-memory column-major matrix.
// Benchmarked both for SuiteSparse matrices and dense scaling series.
void add_dense_methods(std::vector<Row>& rows, std::vector<Raw>& raw, const matrix& A_dense, const matrix& A_original_band, const std::vector<Scalar>& b,
                       std::size_t repeats, bool try_cholesky)
{
    const std::size_t         n  = A_dense.rows();
    const std::size_t         p0 = A_original_band.bi(), q0 = A_original_band.bs();
    const double              amax = max_abs_matrix(A_original_band);
    const std::vector<Scalar> base = dense_data(A_dense);

    add_benchmark(rows,
                  raw,
                  "custom_lu_dense",
                  "original",
                  n,
                  p0,
                  q0,
                  lu_flops(n, n - 1, n - 1),
                  NA,
                  &A_original_band,
                  &b,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      matrix                work = A_dense;
                      std::vector<unsigned> pivots;
                      auto                  t0 = Clock::now();
                      lss::direct::lu::factorize(work, pivots);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::lu::solve(work, pivots, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), growth_matrix(work, amax, n - 1), false, std::move(x) };
                  });

    add_benchmark(rows,
                  raw,
                  "lapack_dgetrf",
                  "original",
                  n,
                  p0,
                  q0,
                  lu_flops(n, n - 1, n - 1),
                  NA,
                  &A_original_band,
                  &b,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::dense_lu(base, b, n, amax)); });

    if (!try_cholesky)
        return;
    add_benchmark(rows,
                  raw,
                  "custom_cholesky_dense",
                  "original",
                  n,
                  p0,
                  q0,
                  chol_flops(n, n - 1),
                  NA,
                  &A_original_band,
                  &b,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      matrix work = A_dense;
                      auto   t0   = Clock::now();
                      lss::direct::cholesky::factorize(work);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::cholesky::solve(work, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), NA, false, std::move(x) };
                  });
    add_benchmark(rows,
                  raw,
                  "lapack_dpotrf",
                  "original",
                  n,
                  p0,
                  q0,
                  chol_flops(n, n - 1),
                  NA,
                  &A_original_band,
                  &b,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::dense_cholesky(base, b, n)); });
}

// Original or RCM band methods; each method uses its matching b and matrix.
void add_band_methods(std::vector<Row>& rows, std::vector<Raw>& raw, const matrix& band, const matrix& original_for_error, const std::vector<Scalar>& rhs,
                      const std::vector<Scalar>& b_original, const std::vector<std::size_t> *permutation, const std::string& variant, std::size_t repeats,
                      bool try_cholesky)
{
    const std::size_t n = band.rows(), kl = band.bi(), ku = band.bs();
    const double      amax = max_abs_matrix(original_for_error);

    add_benchmark(rows,
                  raw,
                  "custom_lu_band",
                  variant,
                  n,
                  kl,
                  ku,
                  lu_flops(n, kl, ku),
                  NA,
                  &original_for_error,
                  &b_original,
                  permutation,
                  repeats,
                  [&]()
                  {
                      matrix work = band;
                      auto   t0   = Clock::now();
                      lss::direct::lu::band::factorize(work);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::lu::band::solve(work, rhs);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), growth_matrix(work, amax, ku), false, std::move(x) };
                  });

    add_benchmark(rows,
                  raw,
                  "lapack_dgbtrf",
                  variant,
                  n,
                  kl,
                  ku,
                  lu_flops(n, kl, ku),
                  NA,
                  &original_for_error,
                  &b_original,
                  permutation,
                  repeats,
                  [&]() { return as_trial(lapack_tests::band_lu(band, rhs, amax)); });

    if (!try_cholesky)
        return;

    add_benchmark(rows,
                  raw,
                  "custom_cholesky_band",
                  variant,
                  n,
                  kl,
                  ku,
                  chol_flops(n, kl),
                  NA,
                  &original_for_error,
                  &b_original,
                  permutation,
                  repeats,
                  [&]()
                  {
                      matrix work = band;
                      auto   t0   = Clock::now();
                      lss::direct::cholesky::band::factorize(work);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::cholesky::band::solve(work, rhs);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), NA, false, std::move(x) };
                  });

    add_benchmark(rows,
                  raw,
                  "lapack_dpbtrf",
                  variant,
                  n,
                  kl,
                  ku,
                  chol_flops(n, kl),
                  NA,
                  &original_for_error,
                  &b_original,
                  permutation,
                  repeats,
                  [&]() { return as_trial(lapack_tests::band_cholesky(band, rhs)); });
}

void add_thomas(std::vector<Row>& rows, std::vector<Raw>& raw, const matrix& A, const std::vector<Scalar>& b, std::size_t repeats)
{
    const std::size_t n = A.rows();
    if (A.bi() > 1 || A.bs() > 1)
        return;
    std::vector<Scalar> l(n - 1), d(n), u(n - 1);
    for (std::size_t i = 0; i < n; ++i)
        d[i] = A.get(static_cast<unsigned>(i), static_cast<unsigned>(i));
    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        l[i] = A.get(static_cast<unsigned>(i + 1), static_cast<unsigned>(i));
        u[i] = A.get(static_cast<unsigned>(i), static_cast<unsigned>(i + 1));
    }
    add_benchmark(rows,
                  raw,
                  "custom_thomas",
                  "original",
                  n,
                  A.bi(),
                  A.bs(),
                  0.0,
                  8.0 * n,
                  &A,
                  &b,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      auto lf = l, df = d;
                      auto t0 = Clock::now();
                      lss::direct::thomas::factorize(lf, df, u);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::thomas::solve(lf, df, u, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), NA, false, std::move(x) };
                  });
    add_benchmark(rows,
                  raw,
                  "lapack_dgtsv",
                  "original",
                  n,
                  A.bi(),
                  A.bs(),
                  0.0,
                  8.0 * n,
                  &A,
                  &b,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::tridiagonal(l, d, u, b)); });
}

void print_results(const std::vector<Row>& rows)
{
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\nMethod                       Variant      Fact(s)      Solve(s)     "
                 "Backward        Forward          Status\n";
    std::cout << "(dgtsv: factor/solve separately unavailable; see total_s in CSV)\n";
    for (const auto& r : rows)
        std::cout << std::left << std::setw(28) << r.method << std::setw(13) << r.variant << std::right << std::setw(12) << r.factor_s << std::setw(13)
                  << r.solve_s << std::setw(16) << r.backward_error << std::setw(16) << r.forward_error << "  " << r.status << '\n';
}

// Scaled Thomas test: zero dense matrices are allocated.
void thomas_scale(std::size_t n, const std::filesystem::path& dir, std::size_t repeats)
{
    if (n == 0)
        throw std::invalid_argument("n must be > 0");
    std::vector<Scalar> lower(n - 1, -1.0), diagonal(n, 2.0), upper(n - 1, -1.0);
    std::vector<Scalar> b(n, 0.0);
    b[0] += 1.0;
    b[n - 1] += 1.0;
    std::vector<Row> rows;
    std::vector<Raw> raw;

    add_benchmark(rows,
                  raw,
                  "custom_thomas",
                  "tridiag_scale",
                  n,
                  1,
                  1,
                  0.0,
                  8.0 * n,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      auto lf = lower, df = diagonal;
                      auto t0 = Clock::now();
                      lss::direct::thomas::factorize(lf, df, upper);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::thomas::solve(lf, df, upper, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), NA, false, std::move(x) };
                  });
    add_benchmark(rows,
                  raw,
                  "lapack_dgtsv",
                  "tridiag_scale",
                  n,
                  1,
                  1,
                  0.0,
                  8.0 * n,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::tridiagonal(lower, diagonal, upper, b)); });
    for (auto& row : rows)
        if (row.status == "ok")
        {
            row.forward_error = forward_error(row.x_original);
            // residual tridiag(-1,2,-1), independently of matrix class.
            double      normr = 0.0, normx = 0.0;
            const auto& x = row.x_original;
            for (std::size_t i = 0; i < n; ++i)
            {
                double ax = 2.0 * x[i];
                if (i > 0)
                    ax -= x[i - 1];
                if (i + 1 < n)
                    ax -= x[i + 1];
                normr = std::max(normr, std::abs(b[i] - ax));
                normx = std::max(normx, std::abs(x[i]));
            }
            row.backward_error = normx > 0 ? normr / (4.0 * normx) : NA;
        }
    compare_to_reference(rows);  // No dense reference in this scale test.
    for (auto& row : rows)
        if (row.method == "custom_thomas" && row.status == "ok")
            for (const auto& lap : rows)
                if (lap.method == "lapack_dgtsv" && lap.status == "ok")
                    row.lapack_diff = difference_inf(row.x_original, lap.x_original);
    save_rows(dir, rows, raw);
    print_results(rows);
}

// Generate symmetric strictly diagonally dominant dense matrix: SPD.
// O(n^2) setup is outside the timed region.
void dense_scale(std::size_t n, const std::filesystem::path& dir, std::size_t repeats)
{
    if (n == 0)
        throw std::invalid_argument("n must be > 0");
    std::vector<Scalar>                    base(n * n, 0.0);
    std::vector<double>                    sum(n, 0.0);
    std::mt19937_64                        rng(20261101);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = j + 1; i < n; ++i)
        {
            const double v  = distribution(rng) / static_cast<double>(n);
            base[j * n + i] = v;
            base[i * n + j] = v;
            sum[i] += std::abs(v);
            sum[j] += std::abs(v);
        }
    for (std::size_t i = 0; i < n; ++i)
        base[i * n + i] = sum[i] + 1.0;

    std::vector<Scalar> b(n, 0.0);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
            b[i] += base[j * n + i];

    const double     flops_lu = lu_flops(n, n - 1, n - 1), flops_chol = chol_flops(n, n - 1);
    std::vector<Row> rows;
    std::vector<Raw> raw;

    // Initialize a dense matrix without parsing n*n Matrix Market entries:
    // read an n-entry diagonal skeleton, then assign every dense coefficient
    // through the project's operator()(i,j). Generation is outside timing.
    // The reader's printed metadata refer to the skeleton, not the final
    // in-memory random SPD matrix.
    const auto path = dir / "generated_dense_skeleton.mtx";
    {
        auto f = new_csv(path);
        f << "%%MatrixMarket matrix coordinate real general\n";
        f << n << " " << n << " " << n << "\n";
        for (std::size_t i = 0; i < n; ++i)
            f << i + 1 << " " << i + 1 << " 1.0\n";
    }
    matrix A(true);
    A.read_from_file(path.string());
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
            A(i, j) = base[j * n + i];
    // Keep generated matrix available for reproducibility.
    const double amax = *std::max_element(sum.begin(), sum.end()) + 1.0;

    add_benchmark(rows,
                  raw,
                  "custom_lu_dense",
                  "dense_scale",
                  n,
                  n - 1,
                  n - 1,
                  flops_lu,
                  NA,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      matrix                work = A;
                      std::vector<unsigned> pivots;
                      auto                  t0 = Clock::now();
                      lss::direct::lu::factorize(work, pivots);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::lu::solve(work, pivots, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), growth_matrix(work, amax, n - 1), false, std::move(x) };
                  });
    add_benchmark(rows,
                  raw,
                  "lapack_dgetrf",
                  "dense_scale",
                  n,
                  n - 1,
                  n - 1,
                  flops_lu,
                  NA,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::dense_lu(base, b, n, amax)); });
    add_benchmark(rows,
                  raw,
                  "custom_cholesky_dense",
                  "dense_scale",
                  n,
                  n - 1,
                  n - 1,
                  flops_chol,
                  NA,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]()
                  {
                      matrix work = A;
                      auto   t0   = Clock::now();
                      lss::direct::cholesky::factorize(work);
                      auto t1 = Clock::now();
                      auto x  = lss::direct::cholesky::solve(work, b);
                      auto t2 = Clock::now();
                      return Trial{ seconds(t0, t1), seconds(t1, t2), NA, false, std::move(x) };
                  });
    add_benchmark(rows,
                  raw,
                  "lapack_dpotrf",
                  "dense_scale",
                  n,
                  n - 1,
                  n - 1,
                  flops_chol,
                  NA,
                  nullptr,
                  nullptr,
                  nullptr,
                  repeats,
                  [&]() { return as_trial(lapack_tests::dense_cholesky(base, b, n)); });

    const double kappa1 = lapack_tests::estimate_condition_number_1(base, n);
    for (auto& row : rows)
        if (row.status == "ok")
        {
            row.forward_error = forward_error(row.x_original);
            row.kappa1        = kappa1;
            row.kappa_u       = kappa1 * 1.1e-16;
            // residual computed against the generated dense array.
            std::vector<double> r(n, 0.0), row_sum(n, 0.0);
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = 0; i < n; ++i)
                {
                    const double a = base[j * n + i];
                    r[i] += a * row.x_original[j];
                    row_sum[i] += std::abs(a);
                }
            double err = 0.0, na = 0.0, nx = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                err = std::max(err, std::abs(r[i] - b[i]));
                na  = std::max(na, row_sum[i]);
                nx  = std::max(nx, std::abs(row.x_original[i]));
            }
            row.backward_error = (na * nx > 0) ? err / (na * nx) : NA;
        }
    for (auto& row : rows)
        if (row.method != "lapack_dgetrf" && row.status == "ok")
            for (const auto& ref : rows)
                if (ref.method == "lapack_dgetrf" && ref.status == "ok")
                    row.lapack_diff = difference_inf(row.x_original, ref.x_original);
    save_rows(dir, rows, raw);
    print_results(rows);
}
}  // namespace

namespace benchmarks
{
void run_matrix(const std::filesystem::path& matrix_file, const std::filesystem::path& output_dir, const Options& options)
{
    if (options.repeats < 5)
        throw std::invalid_argument("At least five repeats are required by the assignment");
    std::filesystem::create_directories(output_dir);

    matrix A_original_band(false);
    A_original_band.read_from_file(matrix_file.string());
    const std::size_t n = A_original_band.rows();
    if (n == 0 || A_original_band.cols() != n)
        throw std::invalid_argument("Matrix must be nonempty and square");

    const bool        is_symmetric = symmetric(A_original_band);
    const std::size_t nnz          = count_nonzero(A_original_band);
    const auto        b            = matrix_times_ones(A_original_band);

    // Report median graph and RCM times, using the same number of
    // repetitions required for the solver benchmarks.
    std::vector<double> graph_samples, rcm_samples;
    graph_samples.reserve(options.repeats);
    rcm_samples.reserve(options.repeats);
    reordering::Graph        graph;
    std::vector<std::size_t> permutation;
    for (std::size_t rep = 0; rep < options.repeats; ++rep)
    {
        const auto t0 = Clock::now();
        graph         = reordering::build_graph(A_original_band);
        const auto t1 = Clock::now();
        permutation   = reordering::reverse_cuthill_mckee(graph);
        const auto t2 = Clock::now();
        graph_samples.push_back(seconds(t0, t1));
        rcm_samples.push_back(seconds(t1, t2));
    }
    if (permutation.size() != n)
        throw std::runtime_error("Invalid RCM size");
    const double graph_s = median(graph_samples);
    const double rcm_s   = median(rcm_samples);

    const auto rcm_path       = output_dir / "matrix_rcm.mtx";
    const auto t_export_start = Clock::now();
    reordering::save_permuted_matrix_mtx(A_original_band, permutation, rcm_path);
    const auto t_export_end = Clock::now();

    matrix     A_rcm(false);
    const auto t_read_start = Clock::now();
    A_rcm.read_from_file(rcm_path.string());
    const auto   t_read_end           = Clock::now();
    const double permutation_export_s = seconds(t_export_start, t_export_end);
    const double reload_s             = seconds(t_read_start, t_read_end);
    const double reorder_s            = graph_s + rcm_s + permutation_export_s + reload_s;
    auto         prep                 = new_csv(output_dir / "preprocessing.csv");
    prep << "stage,repeats,time_s,aggregation\n"
         << std::setprecision(17) << "graph," << options.repeats << "," << graph_s << ",median\n"
         << "rcm_ordering," << options.repeats << "," << rcm_s << ",median\n"
         << "permutation_export,1," << permutation_export_s << ",single\n"
         << "reload_rcm_matrix,1," << reload_s << ",single\n"
         << "preprocessing_total_estimate,0," << reorder_s << ",mixed_estimate\n";
    prep.close();
    auto prep_raw = new_csv(output_dir / "preprocessing_raw.csv");
    prep_raw << "stage,repeat,time_s\n" << std::setprecision(17);
    for (std::size_t rep = 0; rep < options.repeats; ++rep)
    {
        prep_raw << "graph," << rep << ',' << graph_samples[rep] << '\n';
        prep_raw << "rcm_ordering," << rep << ',' << rcm_samples[rep] << '\n';
    }
    prep_raw.close();
    matrix A_dense(true);
    A_dense.read_from_file(matrix_file.string());
    const double        kappa1 = lapack_tests::estimate_condition_number_1(dense_data(A_dense), n);
    std::vector<Scalar> b_rcm(n);
    for (std::size_t i = 0; i < n; ++i)
        b_rcm[i] = b[permutation[i]];

    const std::size_t p0 = A_original_band.bi(), q0 = A_original_band.bs();
    const std::size_t p1 = A_rcm.bi(), q1 = A_rcm.bs();
    const std::size_t mem0 = n * (p0 + q0 + 1) * sizeof(Scalar);
    const std::size_t mem1 = n * (p1 + q1 + 1) * sizeof(Scalar);

    auto metadata = new_csv(output_dir / "matrix_metadata.csv");
    metadata << "n,nnz,density,symmetric,lower_before,upper_before,"
                "lower_after,upper_after,band_bytes_before,band_bytes_after,dense_bytes,kappa1_est,kappa_u\n";
    metadata << n << ',' << nnz << ',' << std::setprecision(17) << double(nnz) / double(n * n) << ',' << (is_symmetric ? 1 : 0) << ',' << p0 << ',' << q0 << ','
             << p1 << ',' << q1 << ',' << mem0 << ',' << mem1 << ',' << n * n * sizeof(Scalar) << ',' << kappa1 << ',' << kappa1 * 1.1e-16 << '\n';
    metadata.close();

    if (options.spy_plots)
    {
        plot(matrix_file, output_dir);
        plot(rcm_path, output_dir);
    }

    std::cout << "\nMATRIX: " << matrix_file << "\n"
              << "n=" << n << " nnz=" << nnz << " symmetric=" << is_symmetric << '\n'
              << "band before (" << p0 << "," << q0 << ") -> after (" << p1 << "," << q1 << ")\n"
              << "band storage bytes: " << mem0 << " -> " << mem1 << '\n';

    std::vector<Row> rows;
    std::vector<Raw> raw;
    add_dense_methods(rows, raw, A_dense, A_original_band, b, options.repeats, is_symmetric);
    add_band_methods(rows, raw, A_original_band, A_original_band, b, b, nullptr, "original", options.repeats, is_symmetric);
    add_band_methods(rows, raw, A_rcm, A_original_band, b_rcm, b, &permutation, "rcm", options.repeats, is_symmetric);
    add_thomas(rows, raw, A_original_band, b, options.repeats);
    compare_to_reference(rows);
    for (auto& row : rows)
    {
        row.kappa1  = kappa1;
        row.kappa_u = kappa1 * 1.1e-16;
    }
    save_rows(output_dir, rows, raw);

    auto p3 = new_csv(output_dir / "part3_3.csv");
    write_header(p3);
    // Cholesky only when BOTH dense and RCM band factorizations succeeded.
    auto success = [&](const std::string& method, const std::string& variant)
    {
        for (const auto& row : rows)
            if (row.method == method && row.variant == variant && row.status == "ok")
                return true;
        return false;
    };
    const bool        use_cholesky = is_symmetric && success("custom_cholesky_dense", "original") && success("custom_cholesky_band", "rcm");
    const std::string dense_name   = use_cholesky ? "custom_cholesky_dense" : "custom_lu_dense";
    const std::string band_name    = use_cholesky ? "custom_cholesky_band" : "custom_lu_band";
    const Row        *dense3       = nullptr;
    const Row        *band3        = nullptr;
    for (const auto& row : rows)
    {
        if (row.method == dense_name && row.variant == "original")
        {
            write_row(p3, row);
            dense3 = &row;
        }
        if (row.method == band_name && row.variant == "rcm")
        {
            write_row(p3, row);
            band3 = &row;
        }
    }
    auto comparison = new_csv(output_dir / "part3_comparison.csv");
    comparison << "n,method,reorder_s,graph_s,rcm_ordering_s,export_s,reload_s,dense_factor_s,rcm_band_factor_s,"
                  "measured_factor_speedup,predicted_flop_ratio,dense_flops,band_flops,"
                  "band_bytes_before,band_bytes_after\n"
               << std::setprecision(17);
    if (dense3 && band3 && dense3->status == "ok" && band3->status == "ok")
    {
        comparison << n << ',' << (use_cholesky ? "cholesky" : "lu") << ',' << reorder_s << ',' << graph_s << ',' << rcm_s << ',' << permutation_export_s << ','
                   << reload_s << ',' << dense3->factor_s << ',' << band3->factor_s << ',' << (band3->factor_s > 0 ? dense3->factor_s / band3->factor_s : NA)
                   << ',' << (band3->factor_flops > 0 ? dense3->factor_flops / band3->factor_flops : NA) << ',' << dense3->factor_flops << ','
                   << band3->factor_flops << ',' << mem0 << ',' << mem1 << '\n';
        std::cout << "RCM preprocessing time (graph + BFS + export + reload): " << reorder_s << " s\n"
                  << "Measured factor speedup (dense/band): " << (band3->factor_s > 0 ? dense3->factor_s / band3->factor_s : NA) << "\n"
                  << "Predicted factor work ratio (dense/band): " << (band3->factor_flops > 0 ? dense3->factor_flops / band3->factor_flops : NA) << "\n";
    }
    else
        std::cout << "Part 3.3 comparison unavailable: a selected factorization failed.\n";

    print_results(rows);
    std::cout << "CSV files: " << output_dir / "summary.csv" << ", " << output_dir / "raw.csv" << ", " << output_dir / "part3_3.csv" << '\n';
}

void run_thomas_scale(std::size_t n, const std::filesystem::path& output_dir, const Options& options)
{
    if (options.repeats < 5)
        throw std::invalid_argument("repeats must be at least 5");
    std::filesystem::create_directories(output_dir);
    thomas_scale(n, output_dir, options.repeats);
}

void run_dense_scale(std::size_t n, const std::filesystem::path& output_dir, const Options& options)
{
    if (options.repeats < 5)
        throw std::invalid_argument("repeats must be at least 5");
    std::filesystem::create_directories(output_dir);
    dense_scale(n, output_dir, options.repeats);
}
}  // namespace benchmarks
