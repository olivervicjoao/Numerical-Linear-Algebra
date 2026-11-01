#include "../include/utils.hpp"

void save_vector_mtx(const std::vector<real>& x, const std::filesystem::path& output_dir, const std::string& filename)
{
    // Caminho completo do arquivo
    const std::filesystem::path output_file = output_dir / filename;

    std::ofstream file(output_file);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open output file: " + output_file.string());
    }

    // Cabeçalho Matrix Market
    file << "%%MatrixMarket matrix array real general\n";
    file << "% Solution vector\n";

    // Dimensões: n linhas e 1 coluna
    file << x.size() << " 1\n";

    // Precisão suficiente para preservar os valores de double
    file << std::scientific << std::setprecision(std::numeric_limits<real>::max_digits10);

    // Salvar os valores
    for (const real value : x)
    {
        file << value << '\n';
    }

    if (!file)
    {
        throw std::runtime_error("Error writing file: " + output_file.string());
    }
}
