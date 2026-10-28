#include "../include/matrix.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

matrix::matrix(bool dense)
{
    /*
     * Para quem estiver lendo isso no futuro (especialmente pra mim pq eu esqueço as coisas):
     *      Eu sei que esse construtor poderia estar melhor, eu sei que provavelmente vc vai achar que consegue otimiza-lo! Mas não!
     *      Ele funciona, e funciona rápido! Lembre-se da frase: "If It Works Dont Touch It"!
     */

    this->_dense_ = dense;

    this->_rows_            = 0;
    this->_cols_            = 0;
    this->_nnz_             = 0;
    this->_lower_bandwidth_ = 0;
    this->_upper_bandwidth_ = 0;

    this->_density_            = 0.0;
    this->_storage_size_bytes_ = 0.0;
}

double& matrix::operator()(std::size_t i, std::size_t j)
{
    if (i >= _rows_ || j >= _cols_)
        throw std::out_of_range("Matrix index out of range");

    return _data_.at(j * _rows_ + i);
}

const double& matrix::operator()(std::size_t i, std::size_t j) const
{
    if (i >= _rows_ || j >= _cols_)
        throw std::out_of_range("Matrix index out of range");

    return _data_.at(j * _rows_ + i);
}

unsigned matrix::rows() const
{
    return this->_rows_;
}

unsigned matrix::cols() const
{
    return this->_cols_;
}

void matrix::prepare_matrix_obj()
{
    /*
     * Lógica:
     *      -> Se for densa reserva todos os espaços e já coloca os zeros
     *      -> Se não for densa (banda) reserva apenas a banda
     */

    if (this->_dense_ == true)
    {
        this->_data_.resize(this->_rows_ * this->_cols_, 0.0);
    }
    else
    {
        this->_data_.resize(this->_upper_bandwidth_ + this->_lower_bandwidth_ + 1, 0.0);
        // std::cout << "TODO!" << std::endl;
    }
}

void matrix::show()
{
    std::cout << std::fixed << std::setprecision(2);

    for (std::size_t i = 0; i < _rows_; i++)
    {
        std::cout << "[ ";

        for (std::size_t j = 0; j < _cols_; j++)
        {
            std::cout << std::setw(10) << (*this)(i, j) << " ";
        }

        std::cout << " ]" << std::endl;
    }

    std::cout << std::defaultfloat;
    std::cout << std::endl;
}

void matrix::read_from_file(std::string file)
{
    /*
     * Leitura de matriz no formato Matrix Market (.mtx).
     * As triplas (i, j, a_ij) são convertidas de índices
     * iniciados em 1 para índices iniciados em 0.
     *
     * O armazenamento denso é column-major,
     * controlado pelo operator()(i, j).
     */

    std::ifstream input(file);

    if (!input.is_open())
    {
        throw std::runtime_error("Could not open Matrix Market file");
    }

    std::string line;

    // ---------------------------------------
    // 1. Leitura do cabeçalho
    // ---------------------------------------
    if (!std::getline(input, line))
    {
        throw std::runtime_error("Empty Matrix Market file");
    }

    std::istringstream header(line);

    std::string banner, object, format, field, symmetry;

    if (!(header >> banner >> object >> format >> field >> symmetry))
    {
        throw std::runtime_error("Invalid Matrix Market header");
    }

    if (banner != "%%MatrixMarket" || object != "matrix" || format != "coordinate" || field != "real" || (symmetry != "general" && symmetry != "symmetric"))
    {
        throw std::runtime_error("Unsupported Matrix Market format");
    }

    this->_symetric_ = (symmetry == "symmetric");

    // ---------------------------------------
    // 2. Leitura das dimensões
    // ---------------------------------------
    bool found_dimensions = false;

    while (std::getline(input, line))
    {
        auto first = line.find_first_not_of(" \t\r\n");

        if (first == std::string::npos || line[first] == '%')
            continue;

        found_dimensions = true;
        break;
    }

    if (!found_dimensions)
    {
        throw std::runtime_error("Missing matrix dimensions");
    }

    std::stringstream stream(line);

    unsigned file_entries;

    if (!(stream >> this->_rows_ >> this->_cols_ >> file_entries))
    {
        throw std::runtime_error("Invalid matrix dimensions");
    }

    if (this->_rows_ == 0 || this->_cols_ == 0)
    {
        throw std::runtime_error("Invalid matrix size");
    }

    if (this->_symetric_ && this->_rows_ != this->_cols_)
    {
        throw std::runtime_error("Symmetric matrix must be square");
    }

    // ---------------------------------------
    // 3. Inicialização dos atributos
    // ---------------------------------------
    this->_nnz_             = 0;
    this->_density_         = 0.0;
    this->_lower_bandwidth_ = 0;
    this->_upper_bandwidth_ = 0;

    unsigned i, j;
    real     value;

    std::vector<std::pair<std::pair<unsigned, unsigned>, real>> info;

    info.reserve(file_entries);

    // ---------------------------------------
    // 4. Leitura das triplas
    // ---------------------------------------
    while (info.size() < file_entries && std::getline(input, line))
    {
        auto first = line.find_first_not_of(" \t\r\n");

        if (first == std::string::npos || line[first] == '%')
            continue;

        stream.clear();
        stream.str(line);

        if (!(stream >> i >> j >> value))
        {
            throw std::runtime_error("Invalid matrix entry");
        }

        if (i == 0 || j == 0 || i > _rows_ || j > _cols_)
        {
            throw std::runtime_error("Invalid matrix indices");
        }

        info.emplace_back(std::make_pair(i - 1, j - 1), value);

        if (value != 0.0)
        {
            this->_nnz_++;

            if (this->_symetric_ && i != j)
                this->_nnz_++;

            if (i > j)
            {
                this->_lower_bandwidth_ = std::max(this->_lower_bandwidth_, i - j);
            }
            else if (j > i)
            {
                this->_upper_bandwidth_ = std::max(this->_upper_bandwidth_, j - i);
            }
        }
    }

    if (info.size() != file_entries)
    {
        throw std::runtime_error("Incomplete Matrix Market file");
    }

    // ---------------------------------------
    // 5. Tratamento da simetria
    // ---------------------------------------
    if (this->_symetric_)
    {
        unsigned b = std::max(this->_lower_bandwidth_, this->_upper_bandwidth_);

        this->_lower_bandwidth_ = b;
        this->_upper_bandwidth_ = b;
    }

    // ---------------------------------------
    // 6. Cálculo da densidade
    // ---------------------------------------
    this->_density_ = static_cast<real>(this->_nnz_) / (static_cast<real>(this->_rows_) * this->_cols_);

    // ---------------------------------------
    // 7. Preparação do armazenamento
    // ---------------------------------------
    this->prepare_matrix_obj();

    // ---------------------------------------
    // 8. Transferência para a matriz
    // ---------------------------------------
    for (const auto& entry : info)
    {
        unsigned i     = entry.first.first;
        unsigned j     = entry.first.second;
        real     value = entry.second;

        (*this)(i, j) = value;

        if (this->_symetric_ && i != j)
        {
            (*this)(j, i) = value;
        }
    }

    this->_storage_size_bytes_ = this->_data_.capacity() * sizeof(real);

    // ---------------------------------------
    // 9. Estimativa do número de condição
    // ---------------------------------------

    if (this->_rows_ != this->_cols_)
    {
        throw std::runtime_error("Condition number requires a square matrix");
    }

    lapack_int n = static_cast<lapack_int>(this->_rows_);

    // 9.1. Norma 1 da matriz original
    real anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, this->_data_.data(), n);

    // 9.2. Cópia para não modificar a matriz original
    std::vector<real> data_copy = this->_data_;

    std::vector<lapack_int> ipiv(n);

    // 9.3. Fatoração LU
    lapack_int lapack_info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, data_copy.data(), n, ipiv.data());

    if (lapack_info < 0)
    {
        throw std::runtime_error("LAPACKE_dgetrf failed");
    }

    // 9.4. Estimativa de kappa_1(A)
    real k1 = std::numeric_limits<real>::infinity();

    if (lapack_info == 0)
    {
        real rcond = 0.0;

        lapack_info = LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', n, data_copy.data(), n, anorm, &rcond);

        if (lapack_info != 0)
        {
            throw std::runtime_error("LAPACKE_dgecon failed");
        }

        if (rcond > 0.0)
        {
            k1 = 1.0 / rcond;
        }
    }

    // ---------------------------------------
    // 10. Informações da matriz
    // ---------------------------------------
    std::cout << std::endl;
    std::cout << "Rows: " << this->_rows_ << std::endl;
    std::cout << "Columns: " << this->_cols_ << std::endl;
    std::cout << "NNZ: " << this->_nnz_ << std::endl;
    std::cout << "Density: " << this->_density_ << std::endl;

    std::cout << "Symmetric: " << std::boolalpha << this->_symetric_ << std::endl;

    std::cout << "Lower bandwidth: " << this->_lower_bandwidth_ << std::endl;

    std::cout << "Upper bandwidth: " << this->_upper_bandwidth_ << std::endl;

    std::cout << "Storage in memory: " << std::scientific << std::setprecision(3) << this->_storage_size_bytes_ / (1024.0 * 1024.0) << " MiB"
              << std::defaultfloat << std::endl;

    std::cout << "Condition number (k1): " << k1 << std::endl;

    std::cout << std::endl;
}
