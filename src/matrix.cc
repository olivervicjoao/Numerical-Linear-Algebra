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

real matrix::get(std::size_t i, std::size_t j) const
{
    // Verificar se os índices são válidos
    if (i >= this->_rows_ || j >= this->_cols_)
    {
        throw std::out_of_range("Matrix index out of range!");
    }

    // Matriz densa
    if (this->_dense_)
    {
        return (*this)(i, j);
    }

    // Matriz em banda: verificar se o elemento
    // está fora da banda inferior
    if (i > j && (i - j) > this->_lower_bandwidth_)
    {
        return 0.0;
    }

    // Verificar se está fora da banda superior
    if (j > i && (j - i) > this->_upper_bandwidth_)
    {
        return 0.0;
    }

    // Elemento dentro da banda
    return (*this)(i, j);
}

real& matrix::operator()(std::size_t i, std::size_t j)
{
    if (i >= _rows_ || j >= _cols_)
        throw std::out_of_range("Matrix index out of range");

    if (_dense_)
        return _data_.at(j * static_cast<std::size_t>(_rows_) + i);

    // Fora da banda, A(i,j) = 0 estruturalmente; nao pode ser escrito.
    if ((i >= j && i - j > _lower_bandwidth_) || (j > i && j - i > _upper_bandwidth_))
    {
        throw std::out_of_range("Cannot modify an element outside the band");
    }

    const std::size_t stride   = static_cast<std::size_t>(_lower_bandwidth_) + _upper_bandwidth_ + 1;
    const std::size_t band_row = (i >= j) ? static_cast<std::size_t>(_upper_bandwidth_) + (i - j) : static_cast<std::size_t>(_upper_bandwidth_) - (j - i);

    return _data_.at(j * stride + band_row);
}

const real& matrix::operator()(std::size_t i, std::size_t j) const
{
    if (i >= _rows_ || j >= _cols_)
        throw std::out_of_range("Matrix index out of range");

    if (_dense_)
        return _data_.at(j * static_cast<std::size_t>(_rows_) + i);

    if ((i >= j && i - j > _lower_bandwidth_) || (j > i && j - i > _upper_bandwidth_))
    {
        // Referencia constante segura: o zero nao pode ser alterado.
        static const real zero = 0.0;
        return zero;
    }

    const std::size_t stride   = static_cast<std::size_t>(_lower_bandwidth_) + _upper_bandwidth_ + 1;
    const std::size_t band_row = (i >= j) ? static_cast<std::size_t>(_upper_bandwidth_) + (i - j) : static_cast<std::size_t>(_upper_bandwidth_) - (j - i);

    return _data_.at(j * stride + band_row);
}

unsigned matrix::rows() const
{
    return this->_rows_;
}

unsigned matrix::cols() const
{
    return this->_cols_;
}

unsigned matrix::bi() const
{
    return this->_lower_bandwidth_;
}

unsigned matrix::bs() const
{
    return this->_upper_bandwidth_;
}

void matrix::prepare_matrix_obj()
{
    /*
     * Lógica:
     *      -> Se for densa reserva todos os espaços e já coloca os zeros
     *      -> Se não for densa (banda) reserva apenas a banda
     */
    std::size_t count = 0;

    if (this->_dense_)
    {
        count = static_cast<std::size_t>(this->_rows_) * this->_cols_;
    }
    else
    {
        const std::size_t stride = static_cast<std::size_t>(this->_lower_bandwidth_) + this->_upper_bandwidth_ + 1;
        count                    = stride * this->_cols_;
    }

    std::vector<real> fresh_data(count, 0.0);
    this->_data_.swap(fresh_data);

    this->_storage_size_bytes_ = this->_data_.capacity() * sizeof(real);
}

void matrix::show()
{
    // Preserva os ajustes de formato de std::cout.
    const auto old_flags     = std::cout.flags();
    const auto old_precision = std::cout.precision();

    std::cout << std::fixed << std::setprecision(2);

    // Forca o operator() const: leituras fora da banda retornam zero.
    const matrix& view = *this;

    for (std::size_t i = 0; i < this->_rows_; ++i)
    {
        std::cout << "[ ";
        for (std::size_t j = 0; j < this->_cols_; ++j)
            std::cout << std::setw(10) << view(i, j) << " ";
        std::cout << "]\n";
    }
    std::cout << std::endl;

    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
}

void matrix::verify_symmetry()
{
    if (this->_rows_ != this->_cols_)
    {
        this->_symetric_ = false;
        return;
    }

    const real    tolerance = 1e-12;
    const matrix& view      = *this;  // Leitura segura fora da banda.

    this->_symetric_ = true;

    for (std::size_t i = 0; i < this->_rows_; ++i)
    {
        for (std::size_t j = i + 1; j < this->_cols_; ++j)
        {
            const real a = view(i, j);
            const real b = view(j, i);

            if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tolerance * std::max({ 1.0, std::abs(a), std::abs(b) }))
            {
                this->_symetric_ = false;
                return;
            }
        }
    }
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
        throw std::runtime_error("Could not open Matrix Market file: " + file);

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("Empty Matrix Market file");

    // 1. Cabecalho Matrix Market.
    std::istringstream header(line);
    std::string        banner, object, format, field, symmetry;
    if (!(header >> banner >> object >> format >> field >> symmetry))
        throw std::runtime_error("Invalid Matrix Market header");

    if (banner != "%%MatrixMarket" || object != "matrix" || format != "coordinate" || field != "real" || (symmetry != "general" && symmetry != "symmetric"))
    {
        throw std::runtime_error("Unsupported Matrix Market format");
    }
    const bool declared_symmetric = (symmetry == "symmetric");

    // 2. Ignora comentarios e localiza a linha das dimensoes.
    bool found_dimensions = false;
    while (std::getline(input, line))
    {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '%')
            continue;
        found_dimensions = true;
        break;
    }
    if (!found_dimensions)
        throw std::runtime_error("Missing Matrix Market dimensions");

    std::istringstream dimensions(line);
    unsigned           file_entries = 0;
    unsigned           rows         = 0;
    unsigned           cols         = 0;
    if (!(dimensions >> rows >> cols >> file_entries) || rows == 0 || cols == 0)
        throw std::runtime_error("Invalid Matrix Market dimensions");
    if (declared_symmetric && rows != cols)
        throw std::runtime_error("A symmetric matrix must be square");

    // 3. Leitura temporaria de triplas (indices convertidos para base zero).
    using entry_type = std::pair<std::pair<unsigned, unsigned>, real>;
    std::vector<entry_type> info;
    info.reserve(file_entries);

    while (info.size() < file_entries && std::getline(input, line))
    {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '%')
            continue;

        std::istringstream entry_stream(line);
        unsigned           i = 0, j = 0;
        real               value = 0.0;
        if (!(entry_stream >> i >> j >> value))
            throw std::runtime_error("Invalid Matrix Market entry");
        if (i == 0 || j == 0 || i > rows || j > cols)
            throw std::runtime_error("Matrix Market index out of range");
        if (!std::isfinite(value))
            throw std::runtime_error("Non-finite value in Matrix Market file");

        info.emplace_back(std::make_pair(i - 1, j - 1), value);
    }

    if (info.size() != file_entries)
        throw std::runtime_error("Incomplete Matrix Market file");

    // 4. Consolida coordenadas repetidas somando os valores.
    // Isso evita contar entradas repetidas como varios nao zeros e
    // trata cancelamentos (por exemplo, 1.0 + (-1.0) = 0.0).
    std::map<std::pair<unsigned, unsigned>, real> entries;
    for (const auto& entry : info)
        entries[entry.first] += entry.second;

    // 5. Caracteriza a matriz completa antes de escolher a alocacao.
    this->_rows_            = rows;
    this->_cols_            = cols;
    this->_nnz_             = 0;
    this->_density_         = 0.0;
    this->_lower_bandwidth_ = 0;
    this->_upper_bandwidth_ = 0;
    this->_symetric_        = declared_symmetric;

    for (const auto& entry : entries)
    {
        const unsigned i     = entry.first.first;
        const unsigned j     = entry.first.second;
        const real     value = entry.second;

        if (value == 0.0)
            continue;

        this->_nnz_ += (declared_symmetric && i != j) ? 2u : 1u;

        if (i > j)
            this->_lower_bandwidth_ = std::max(this->_lower_bandwidth_, i - j);
        else if (j > i)
            this->_upper_bandwidth_ = std::max(this->_upper_bandwidth_, j - i);
    }

    if (declared_symmetric)
    {
        const unsigned b        = std::max(this->_lower_bandwidth_, this->_upper_bandwidth_);
        this->_lower_bandwidth_ = b;
        this->_upper_bandwidth_ = b;
    }

    this->_density_ = static_cast<real>(this->_nnz_) / (static_cast<real>(rows) * cols);

    // 6. Alocacao densa ou em banda (sempre por colunas).
    this->prepare_matrix_obj();

    // 7. Transfere os dados, expandindo a simetria se declarada.
    for (const auto& entry : entries)
    {
        const unsigned i     = entry.first.first;
        const unsigned j     = entry.first.second;
        const real     value = entry.second;

        if (value == 0.0)
            continue;

        (*this)(i, j) = value;
        if (declared_symmetric && i != j)
            (*this)(j, i) = value;
    }

    // 8. Verificacao numerica de simetria (tambem para arquivos 'general').
    this->verify_symmetry();

    // 9. Numero de condicao: LAPACKE opera sobre armazenamento DENSO.
    // Se a matriz atual estiver em banda, monta-se somente uma copia
    // temporaria densa para estimar kappa_1, sem mudar _data_.
    if (rows != cols)
        throw std::runtime_error("Condition number requires a square matrix");
    if (static_cast<std::size_t>(rows) > static_cast<std::size_t>(std::numeric_limits<lapack_int>::max()))
    {
        throw std::runtime_error("Matrix dimension exceeds LAPACKE limits");
    }

    const lapack_int  n = static_cast<lapack_int>(rows);
    std::vector<real> lapack_data;

    if (this->_dense_)
    {
        lapack_data = this->_data_;
    }
    else
    {
        lapack_data.assign(static_cast<std::size_t>(rows) * cols, 0.0);
        const matrix& view = *this;
        for (std::size_t j = 0; j < cols; ++j)
        {
            const std::size_t first_row = (j > _upper_bandwidth_) ? (j - _upper_bandwidth_) : 0;
            const std::size_t last_row  = std::min<std::size_t>(rows - 1, j + _lower_bandwidth_);

            for (std::size_t i = first_row; i <= last_row; ++i)
                lapack_data[j * rows + i] = view(i, j);
        }
    }

    const real anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, lapack_data.data(), n);

    std::vector<lapack_int> ipiv(rows);
    lapack_int              lapack_info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, lapack_data.data(), n, ipiv.data());

    if (lapack_info < 0)
        throw std::runtime_error("LAPACKE_dgetrf: invalid argument");

    real k1 = std::numeric_limits<real>::infinity();
    if (lapack_info == 0)
    {
        real rcond  = 0.0;
        lapack_info = LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', n, lapack_data.data(), n, anorm, &rcond);

        if (lapack_info != 0)
            throw std::runtime_error("LAPACKE_dgecon failed");
        if (rcond > 0.0)
            k1 = 1.0 / rcond;
    }

    // 10. Saida informativa (o formato anterior foi mantido).
    const auto old_flags     = std::cout.flags();
    const auto old_precision = std::cout.precision();

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
    std::cout << "Condition number (k1): " << std::scientific << std::setprecision(3) << k1 << std::defaultfloat << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(old_flags);
    std::cout.precision(old_precision);
}
