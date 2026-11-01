#include "include/errors.hpp"
#include "include/lss.hpp"
#include "include/matrix.hpp"
#include "include/rcm.hpp"
#include "include/utils.hpp"

int main(int argc, char **argv)
{
    std::string path_to_matrix(argv[1]);
    std::string path_to_save(argv[2]);

    matrix m1(true);
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
    matrix A_cholesky_band = m1;

    // Vetor b para a matriz SPD tridiagonal 5x5
    std::vector<real> b(m1.rows(), 1.0);

    // Fatoração e resolução usando Cholesky em banda
    std::vector<real> x = lss::direct::cholesky::solve_cholesky(A_cholesky_band, b);

    // Mostrar solução
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "\nSolution = [ ";

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "[" << x[i] << "] ";
    }

    std::cout << "]\n" << std::endl;

    // Salvar solução
    save_vector_mtx(x, path_to_save, "solution_cholesky_band.mtx");

    // Calcular resíduo relativo com a matriz ORIGINAL
    real residual = lss::validation::relative_residual_inf(m1, x, b);

    std::cout << std::scientific << std::setprecision(12) << "Relative residual = " << residual << '\n';

    return 0;
}
