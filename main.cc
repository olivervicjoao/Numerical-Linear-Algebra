#include "include/errors.hpp"
#include "include/lss.hpp"
#include "include/matrix.hpp"
#include "include/utils.hpp"

int main(int argc, char **argv)
{
    std::string path_to_matrix(argv[1]);
    std::string path_to_save(argv[2]);

    matrix m1(false);
    m1.read_from_file(path_to_matrix);

    // m1.show();

    if (std::atoi(argv[3]) == 1)
    {
        std::string command_for_python = "python3 scripts/spy_plot.py " + path_to_matrix + " -o " + path_to_save;

        int status = std::system(command_for_python.c_str());

        if (status != 0)
        {
            std::cerr << "Error generating sparsity plot!\n";
        }
    }

    std::string command_folder = "python3 scripts/create_folder.py " + path_to_save;

    int status2 = std::system(command_folder.c_str());

    if (status2 != 0)
    {
        std::cerr << "Error creating output directory!\n";
    }

    // Preserva a matriz original
    matrix A_lu_band = m1;

    // Vetor b para a matriz SPD tridiagonal 5x5
    std::vector<real> b = { 1.0, 1.0, 0.0, 1.0, 0.0 };

    // Fatoração LU em banda (in-place)
    lss::direct::lu::band::factorize(A_lu_band);

    // Resolução do sistema LUx = b
    std::vector<real> x = lss::direct::lu::band::solve(A_lu_band, b);

    // Mostrar solução
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "\nSolution = [ ";

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "[" << x[i] << "] ";
    }

    std::cout << "]\n" << std::endl;

    // Salvar solução no formato Matrix Market
    save_vector_mtx(x, path_to_save, "solution_lu_band.mtx");

    // Calcular resíduo relativo usando a matriz ORIGINAL
    real residual = lss::validation::relative_residual_inf(m1, x, b);

    std::cout << std::scientific << std::setprecision(12) << "Relative residual = " << residual << '\n';

    return 0;
}
