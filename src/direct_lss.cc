#include "../include/lss.hpp"
#include <cmath>
#include <cstddef>
#include <vector>
namespace
{
std::size_t find_pivot(const matrix& A, std::size_t k)
{
    if (k >= A.rows())
        throw std::out_of_range("Invalid pivot column");

    std::size_t index     = k;
    real        max_value = std::abs(A(k, k));

    for (std::size_t i = k + 1; i < A.rows(); i++)
    {
        real value = std::abs(A(i, k));

        if (value > max_value)
        {
            max_value = value;
            index     = i;
        }
    }

    return index;
}
}  // namespace

std::vector<real> lss::direct::lu::solve_LU(matrix& A, const std::vector<real>& b)
{
    std::vector<unsigned> pivots;
    lss::direct::lu::factorize(A, pivots);
    std::vector<real> x = lss::direct::lu::solve(A, pivots, b);

    return x;
}

void lss::direct::lu::factorize(matrix& A, std::vector<unsigned>& pivots)
{
    if (A.rows() == 0 || A.rows() != A.cols())
    {
        throw std::invalid_argument("Matrix must be nonempty and square");
    }
    pivots.resize(A.rows());
    std::iota(pivots.begin(), pivots.end(), 0u);

    for (std::size_t k = 0; k < A.rows() - 1; k++)
    {
        // Passo 1: achar o maior pivô
        std::size_t pivot_row = find_pivot(A, k);

        // Checando singularidade
        if (A(pivot_row, k) == 0)
        {
            throw std::runtime_error("Singular matrix!");
        }

        // Passo 2: trocar as linhas
        if (pivot_row != k)
        {
            for (std::size_t j = 0; j < A.cols(); j++)
            {
                std::swap(A(k, j), A(pivot_row, j));
            }
            std::swap(pivots[k], pivots[pivot_row]);
        }

        // Passo 3: calcular os multiplicadores
        const real pivot = A(k, k);

        for (std::size_t i = k + 1; i < A.rows(); i++)
        {
            A(i, k) /= pivot;
        }

        // Passo 4: atualizar a submatriz
        for (std::size_t j = k + 1; j < A.cols(); j++)
        {
            const real ukj = A(k, j);

            for (std::size_t i = k + 1; i < A.rows(); i++)
            {
                A(i, j) -= A(i, k) * ukj;
            }
        }
    }

    // Verificar também o último pivô
    if (A(A.rows() - 1, A.cols() - 1) == 0)
    {
        throw std::runtime_error("Singular matrix!");
    }
}

std::vector<real> lss::direct::lu::solve(const matrix& LU, const std::vector<unsigned>& pivots, const std::vector<real>& b)
{
    const std::size_t n = LU.rows();

    if (n == 0 || LU.cols() != n || b.size() != n || pivots.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    // Passo 1: permutação do vetor b (b' = Pb)
    std::vector<real> x(n);

    for (std::size_t i = 0; i < n; i++)
    {
        if (pivots[i] >= n)
        {
            throw std::out_of_range("Invalid pivot index!");
        }

        x[i] = b[pivots[i]];
    }

    // Passo 2: substituição progressiva (Ly = b')
    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = 0; j < i; j++)
        {
            x[i] -= LU(i, j) * x[j];
        }
    }

    // Passo 3: substituição regressiva (Ux = y)
    for (std::size_t ii = n; ii > 0; --ii)
    {
        const std::size_t i = ii - 1;

        for (std::size_t j = i + 1; j < n; j++)
        {
            x[i] -= LU(i, j) * x[j];
        }

        if (LU(i, i) == 0)
        {
            throw std::runtime_error("Zero diagonal in U!");
        }

        x[i] /= LU(i, i);
    }

    return x;
}

void lss::direct::lu::band::factorize(matrix& A)
{
    std::size_t n = A.rows();
}

std::vector<real> lss::direct::lu::band::solve(const matrix& LU, const std::vector<real>& b)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void lss::direct::cholesky::factorize(matrix& A)
{
    if (A.rows() == 0 || A.rows() != A.cols())
    {
        throw std::invalid_argument("Cholesky requires a non-empty square matrix!");
    }

    for (std::size_t k = 0; k < A.cols(); k++)
    {
        // Detecção de pivô negativo
        if (!(A(k, k) > 0.0) || !std::isfinite(A(k, k)))
        {
            throw std::runtime_error("Cholesky factorization failed: invalid pivot!");
        }

        // Elemento da diagonal principal
        A(k, k) = std::sqrt(A(k, k));

        // Calcular coluna k abaixo da diagonal
        for (std::size_t i = k + 1; i < A.cols(); i++)
        {
            A(i, k) /= A(k, k);
        }

        // Atualizar a submatriz restante
        for (std::size_t j = k + 1; j < A.cols(); j++)
        {
            for (std::size_t i = j; i < A.cols(); i++)
            {
                A(i, j) -= A(i, k) * A(j, k);
            }
        }
    }
}

std::vector<real> lss::direct::cholesky::solve(const matrix& L, const std::vector<real>& b)
{
    const std::size_t n = L.rows();

    if (n == 0 || L.cols() != n || b.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    std::vector<real> x = b;

    // Substituição progressiva: Ly = b
    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = 0; j < i; j++)
        {
            x[i] -= L(i, j) * x[j];
        }

        x[i] /= L(i, i);
    }

    // Substituição regressiva: Lᵀx = y
    for (std::size_t ii = n; ii > 0; --ii)
    {
        const std::size_t i = ii - 1;

        for (std::size_t j = i + 1; j < n; j++)
        {
            x[i] -= L(j, i) * x[j];
        }

        x[i] /= L(i, i);
    }

    return x;
}

std::vector<real> lss::direct::cholesky::solve_cholesky(matrix& A, const std::vector<real>& b)
{
    lss::direct::cholesky::factorize(A);

    return lss::direct::cholesky::solve(A, b);
}

void lss::direct::cholesky::band::factorize(matrix& A)
{
}

std::vector<real> lss::direct::cholesky::band::solve(const matrix& L, const std::vector<real>& b)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<real> lss::direct::thomas::solve_thomas(matrix& A, const std::vector<real>& b)
{
    const std::size_t n = A.rows();

    if (n == 0 || A.cols() != n || b.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    // Extrair as três diagonais
    std::vector<real> lower(n - 1);
    std::vector<real> diagonal(n);
    std::vector<real> upper(n - 1);

    // Diagonal principal
    for (std::size_t i = 0; i < n; ++i)
    {
        diagonal[i] = A(i, i);
    }

    // Diagonais inferior e superior
    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        lower[i] = A(i + 1, i);
        upper[i] = A(i, i + 1);
    }

    // Fatoração de Thomas
    lss::direct::thomas::factorize(lower, diagonal, upper);

    // Resolver o sistema
    return lss::direct::thomas::solve(lower, diagonal, upper, b);
}

void lss::direct::thomas::factorize(std::vector<real>& lower, std::vector<real>& diagonal, const std::vector<real>& upper)
{
    const std::size_t n = diagonal.size();

    if (n == 0 || lower.size() != n - 1 || upper.size() != n - 1)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    for (std::size_t i = 0; i < n - 1; i++)
    {
        if (diagonal[i] == 0.0)
        {
            throw std::runtime_error("Zero pivot in Thomas factorization!");
        }

        // Multiplicador de L
        lower[i] /= diagonal[i];

        // Atualização da diagonal de U
        diagonal[i + 1] -= lower[i] * upper[i];
    }

    if (diagonal[n - 1] == 0.0)
    {
        throw std::runtime_error("Zero last pivot in Thomas factorization!");
    }
}

std::vector<real>
lss::direct::thomas::solve(const std::vector<real>& lower, const std::vector<real>& diagonal, const std::vector<real>& upper, const std::vector<real>& b)
{
    const std::size_t n = diagonal.size();

    if (n == 0 || lower.size() != n - 1 || upper.size() != n - 1 || b.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    std::vector<real> x = b;

    // Substituição progressiva: Ly = b
    for (std::size_t i = 1; i < n; ++i)
    {
        x[i] -= lower[i - 1] * x[i - 1];
    }

    // Substituição regressiva: Ux = y
    x[n - 1] /= diagonal[n - 1];

    for (std::size_t ii = n - 1; ii > 0; --ii)
    {
        const std::size_t i = ii - 1;

        x[i] = (x[i] - upper[i] * x[i + 1]) / diagonal[i];
    }

    return x;
}
