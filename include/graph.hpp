#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "defines.hpp"
#include "matrix.hpp"

namespace reordering
{
using Graph = std::vector<std::vector<std::size_t>>;

// Constrói o grafo associado à matriz
Graph build_graph(const matrix& A);
}  // namespace reordering

#endif
