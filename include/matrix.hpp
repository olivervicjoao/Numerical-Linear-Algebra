#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "defines.hpp"
#include <string>
#include <vector>

class matrix
{
public:
    matrix(bool dense = true);

    void read_from_file(std::string file);

    real get(unsigned i, unsigned j) const;

    unsigned rows() const;
    unsigned cols() const;

    double&       operator()(std::size_t i, std::size_t j);
    const double& operator()(std::size_t i, std::size_t j) const;

    void reverse_Cuthill_McKee();

    void show();

    unsigned bi() const;
    unsigned bs() const;

private:
    unsigned _rows_;
    unsigned _cols_;
    unsigned _nnz_;
    unsigned _lower_bandwidth_;
    unsigned _upper_bandwidth_;

    real _density_;
    real _storage_size_bytes_;

    bool _symetric_;
    bool _dense_;

    void prepare_matrix_obj();

    void verify_symmetry();

    std::vector<real> _data_;
};

#endif
