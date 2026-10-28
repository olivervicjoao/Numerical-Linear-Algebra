#include "include/matrix.hpp"

int main(void)
{
    matrix m1;
    m1.read_from_file("/home/joaooliveira/Projects/Numerical_Linear_Algebra/matrix/NS/5x5_teste/NS_5x5_A.mtx");

    m1.show();

    int status = std::system("python3 scripts/spy_plot.py /home/joaooliveira/Projects/Numerical_Linear_Algebra/matrix/NS/5x5_teste/NS_5x5_A.mtx -o tests/");

    if (status != 0)
    {
        std::cerr << "Error generating sparsity plot!\n";
    }

    return 0;
}
