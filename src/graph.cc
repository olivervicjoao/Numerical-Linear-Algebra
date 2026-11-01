#include "../include/graph.hpp"

reordering::Graph reordering::build_graph(const matrix& A)
{
    const std::size_t n = A.rows();

    if (n == 0 || A.cols() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    const std::size_t p = A.bi();
    const std::size_t q = A.bs();

    Graph graph(n);

    // Percorrer as colunas da matriz
    for (std::size_t j = 0; j < n; j++)
    {
        // Limites da banda na coluna j
        const std::size_t i_min = (j >= q) ? j - q : 0;
        const std::size_t i_max = j + std::min(p, n - 1 - j);

        for (std::size_t i = i_min; i <= i_max; i++)
        {
            // Ignorar a diagonal principal
            if (i == j)
            {
                continue;
            }

            // Existe conexão entre i e j?
            if (A(i, j) != 0.0)
            {
                graph[i].push_back(j);
                graph[j].push_back(i);
            }
        }
    }

    // Remover arestas duplicadas
    for (auto& neighbors : graph)
    {
        std::sort(neighbors.begin(), neighbors.end());

        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    return graph;
}
