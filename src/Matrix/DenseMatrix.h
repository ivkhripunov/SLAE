//
// Created by ivankhripunov on 02.03.23.
//

#ifndef SLAE_DENSEMATRIX_H
#define SLAE_DENSEMATRIX_H

#include <vector>
#include <tuple>

template<typename Type>
class DenseMatrix {
private:
    std::vector<Type> matrix_;
    std::size_t height_, width_;

public:

    DenseMatrix(const std::size_t &height, const std::size_t &width) : height_(height), width_(width),
                                                                       matrix_(std::vector<Type>(height * width)) {};

    DenseMatrix(const std::size_t &height, const std::size_t &width,
                std::vector<std::tuple<std::size_t, std::size_t, Type>> &data) : height_(height), width_(width) {
        matrix_.resize(height * width);

        for (const auto &triplet: data)
            matrix_[width * std::get<0>(triplet) + std::get<1>(triplet)] = std::get<2>(triplet);
    }

    Type &operator()(const std::size_t &i, const std::size_t &j) {

        return matrix_[i * width_ + j];
    };


    [[nodiscard]] const Type &operator()(const std::size_t &i, const std::size_t &j) const {

        return matrix_[i * width_ + j];
    }

    [[nodiscard]] const std::size_t &get_geight() const {
        return height_;
    }

    [[nodiscard]] const std::size_t &get_width() const {
        return width_;
    }

};

#endif //SLAE_DENSEMATRIX_H
