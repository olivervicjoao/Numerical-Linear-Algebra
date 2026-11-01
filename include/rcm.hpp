#ifndef RCM_HPP
#define RCM_HPP

#include "graph.hpp"

namespace reordering
{
// Retorna a permutação RCM
std::vector<std::size_t> reverse_cuthill_mckee(const Graph& graph);

void save_permuted_matrix_mtx(const matrix& A, const std::vector<std::size_t>& permutation, const std::filesystem::path& output_file);
}  // namespace reordering

#endif
