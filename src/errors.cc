#include "../include/errors.hpp"

real lss::validation::relative_residual_inf(const matrix& A, const std::vector<real>& x, const std::vector<real>& b)
{
    const std::size_t n = A.rows();

    if (n == 0 || A.cols() != n || x.size() != n || b.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    std::vector<real> residual(n, 0.0);
    std::vector<real> row_sums(n, 0.0);

    // Calcular Ax e as somas absolutas das linhas de A
    // Percurso por colunas (column-major)
    for (std::size_t j = 0; j < n; j++)
    {
        for (std::size_t i = 0; i < n; i++)
        {
            const real aij = A(i, j);

            residual[i] += aij * x[j];
            row_sums[i] += std::abs(aij);
        }
    }

    real norm_r = 0.0;
    real norm_A = 0.0;
    real norm_x = 0.0;
    real norm_b = 0.0;

    for (std::size_t i = 0; i < n; i++)
    {
        // r = Ax - b
        residual[i] -= b[i];

        norm_r = std::max(norm_r, std::abs(residual[i]));
        norm_A = std::max(norm_A, row_sums[i]);
        norm_x = std::max(norm_x, std::abs(x[i]));
        norm_b = std::max(norm_b, std::abs(b[i]));
    }

    const real denominator = norm_A * norm_x + norm_b;

    if (denominator == 0.0)
    {
        if (norm_r == 0.0)
            return 0.0;

        throw std::runtime_error("Undefined relative residual!");
    }

    return norm_r / denominator;
}
