#ifndef LSS_HPP
#define LSS_HPP

#include "defines.hpp"
#include "matrix.hpp"

namespace lss
{
namespace direct
{

    namespace lu
    {

        std::vector<real> solve_LU(matrix& A, const std::vector<real>& b);
        // LU densa com pivoteamento parcial
        void factorize(matrix& A, std::vector<unsigned>& pivots);

        std::vector<real> solve(const matrix& LU, const std::vector<unsigned>& pivots, const std::vector<real>& b);

        namespace band
        {
            // LU em banda sem pivoteamento
            void factorize(matrix& A);

            std::vector<real> solve(const matrix& LU, const std::vector<real>& b);
        }  // namespace band
    }  // namespace lu

    namespace cholesky
    {
        // Cholesky densa: A = LL^T
        void factorize(matrix& A);

        std::vector<real> solve(const matrix& L, const std::vector<real>& b);

        namespace band
        {
            // Cholesky em banda
            void factorize(matrix& A);

            std::vector<real> solve(const matrix& L, const std::vector<real>& b);
        }  // namespace band
    }  // namespace cholesky

    namespace thomas
    {
        // Fatoração tridiagonal
        void factorize(std::vector<real>& lower, std::vector<real>& diagonal, const std::vector<real>& upper);

        // Resolução utilizando a fatoração
        std::vector<real> solve(const std::vector<real>& lower, const std::vector<real>& diagonal, const std::vector<real>& upper, const std::vector<real>& b);
    }  // namespace thomas
}  // namespace direct

}  // namespace lss

#endif
