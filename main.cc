#include "include/matrix.hpp"
#include <string>

int main(int argc, char **argv)
{
    std::string path_to_matrix(argv[1]);
    std::string path_to_save(argv[2]);

    matrix m1;
    m1.read_from_file(path_to_matrix);

    // m1.show();

    std::string command_for_python = "python3 scripts/spy_plot.py " + path_to_matrix + " -o " + path_to_save;

    int status = std::system(command_for_python.c_str());

    if (status != 0)
    {
        std::cerr << "Error generating sparsity plot!\n";
    }

    return 0;
}
