#include "../include/rcm.hpp"

std::vector<std::size_t> reordering::reverse_cuthill_mckee(const Graph& graph)
{
    const std::size_t n = graph.size();

    // Grau de cada vértice
    std::vector<std::size_t> degree(n);

    for (std::size_t i = 0; i < n; i++)
    {
        degree[i] = graph[i].size();

        for (const std::size_t neighbor : graph[i])
        {
            if (neighbor >= n)
            {
                throw std::invalid_argument("Invalid vertex in adjacency list!");
            }
        }
    }

    // Ordenar os candidatos a vértice inicial por grau
    std::vector<std::size_t> roots(n);

    std::iota(roots.begin(), roots.end(), 0);

    std::sort(roots.begin(),
              roots.end(),
              [&degree](std::size_t a, std::size_t b)
              {
                  if (degree[a] != degree[b])
                  {
                      return degree[a] < degree[b];
                  }

                  return a < b;
              });

    std::vector<bool> visited(n, false);

    std::vector<std::size_t> order;
    order.reserve(n);

    std::queue<std::size_t> queue;

    // Processar todas as componentes conexas
    for (const std::size_t start : roots)
    {
        if (visited[start])
        {
            continue;
        }

        visited[start] = true;
        queue.push(start);

        // Busca em largura (BFS)
        while (!queue.empty())
        {
            const std::size_t v = queue.front();
            queue.pop();

            order.push_back(v);

            // Selecionar vizinhos ainda não visitados
            std::vector<std::size_t> neighbors;

            for (const std::size_t u : graph[v])
            {
                if (!visited[u])
                {
                    neighbors.push_back(u);
                }
            }

            // Ordenar vizinhos por grau crescente
            std::sort(neighbors.begin(),
                      neighbors.end(),
                      [&degree](std::size_t a, std::size_t b)
                      {
                          if (degree[a] != degree[b])
                          {
                              return degree[a] < degree[b];
                          }

                          return a < b;
                      });

            // Inserir os vizinhos na fila
            for (const std::size_t u : neighbors)
            {
                if (!visited[u])
                {
                    visited[u] = true;
                    queue.push(u);
                }
            }
        }
    }

    // Cuthill-McKee -> Reverse Cuthill-McKee
    std::reverse(order.begin(), order.end());

    return order;
}

void reordering::save_permuted_matrix_mtx(const matrix& A, const std::vector<std::size_t>& permutation, const std::filesystem::path& output_file)
{
    const std::size_t n = A.rows();

    if (n == 0 || A.cols() != n || permutation.size() != n)
    {
        throw std::invalid_argument("Incompatible dimensions!");
    }

    // Construir a permutação inversa: original -> novo índice
    std::vector<std::size_t> new_position(n);
    std::vector<bool>        visited(n, false);

    for (std::size_t i = 0; i < n; i++)
    {
        const std::size_t old_index = permutation[i];

        if (old_index >= n || visited[old_index])
        {
            throw std::invalid_argument("Invalid permutation!");
        }

        visited[old_index]      = true;
        new_position[old_index] = i;
    }

    const std::size_t p = A.bi();
    const std::size_t q = A.bs();

    // Percorrer apenas os coeficientes potencialmente não nulos
    auto for_each_nonzero = [&](const auto& visit)
    {
        for (std::size_t j = 0; j < n; j++)
        {
            const std::size_t i_min = (j >= q) ? j - q : 0;
            const std::size_t i_max = j + std::min(p, n - 1 - j);

            for (std::size_t i = i_min; i <= i_max; i++)
            {
                const real value = A(i, j);

                if (value != 0.0)
                {
                    visit(new_position[i], new_position[j], value);
                }
            }
        }
    };

    // Contar os coeficientes não nulos
    std::size_t nnz = 0;

    for_each_nonzero([&](std::size_t, std::size_t, real) { ++nnz; });

    // Criar o arquivo Matrix Market
    std::ofstream file(output_file);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open output file: " + output_file.string());
    }

    file << "%%MatrixMarket matrix coordinate real general\n";
    file << "% Matrix reordered using Reverse Cuthill-McKee\n";
    file << n << " " << n << " " << nnz << "\n";

    file << std::scientific << std::setprecision(std::numeric_limits<real>::max_digits10);

    // Escrever A_rcm = P*A*P^T
    for_each_nonzero(
        [&](std::size_t i, std::size_t j, real value)
        {
            // Matrix Market utiliza índices começando em 1
            file << i + 1 << " " << j + 1 << " " << value << "\n";
        });

    if (!file)
    {
        throw std::runtime_error("Error writing RCM matrix!");
    }
}
