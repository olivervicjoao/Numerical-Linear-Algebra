#include "../include/lss.hpp"

void lss::direct::lu::factorize(matrix& A, std::vector<unsigned>& pivots)
{
}

std::vector<real> lss::direct::lu::solve(const matrix& LU, const std::vector<unsigned>& pivots, const std::vector<real>& b)
{
}

void lss::direct::lu::band::factorize(matrix& A)
{
}

std::vector<real> lss::direct::lu::band::solve(const matrix& LU, const std::vector<real>& b)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void lss::direct::cholesky::factorize(matrix& A)
{
}

std::vector<real> lss::direct::cholesky::solve(const matrix& L, const std::vector<real>& b)
{
}

void lss::direct::cholesky::band::factorize(matrix& A)
{
}

std::vector<real> lss::direct::cholesky::band::solve(const matrix& L, const std::vector<real>& b)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void lss::direct::thomas::factorize(std::vector<real>& lower, std::vector<real>& diagonal, const std::vector<real>& upper)
{
}

std::vector<real>
lss::direct::thomas::solve(const std::vector<real>& lower, const std::vector<real>& diagonal, const std::vector<real>& upper, const std::vector<real>& b)
{
}
