#pragma once

#ifndef BASICMAT_H
#define BASICMAT_H

#include <vector>
#include <stdexcept>
#include <cmath> // For mathematical functions

namespace hps
{

    class BasicMat {
    public:
        // Constructors
        BasicMat();
        BasicMat(int rows, int cols);

        // Destructor
        ~BasicMat();

        // Access element
        float& at(int row, int col);
        const float& at(int row, int col) const;

        // Getters
        int rows() const { return rows_; }
        int cols() const { return cols_; }

        // Data pointer
        float* data() { return data_; }
        const float* data() const { return data_; }

        // Matrix multiplication
        BasicMat operator*(const BasicMat& other) const;

    private:
        int rows_;
        int cols_;
        float* data_;
    };

    // Inline definitions

    inline BasicMat::BasicMat() : rows_(0), cols_(0), data_(nullptr) {}

    inline BasicMat::BasicMat(int rows, int cols)
        : rows_(rows), cols_(cols) {
        data_ = new float[rows * cols]();
    }

    inline BasicMat::~BasicMat() {
        delete[] data_;
    }

    inline float& BasicMat::at(int row, int col) {
        if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        return data_[row * cols_ + col];
    }

    inline const float& BasicMat::at(int row, int col) const {
        if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        return data_[row * cols_ + col];
    }

    inline BasicMat BasicMat::operator*(const BasicMat& other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }

        BasicMat result(rows_, other.cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < other.cols_; ++j) {
                for (int k = 0; k < cols_; ++k) {
                    result.at(i, j) += at(i, k) * other.at(k, j);
                }
            }
        }
        return result;
    }

} // namespace MyNamespace

#endif // BASICMAT_H

