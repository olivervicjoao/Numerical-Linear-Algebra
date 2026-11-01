#ifndef ERRORS_HPP
#define ERRORS_HPP
#include "defines.hpp"
#include "matrix.hpp"

namespace lss::validation
{
real relative_residual_inf(const matrix& A, const std::vector<real>& x, const std::vector<real>& b);
}

#endif
