//
// Created by ivankhripunov on 17.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <vector>
#include <tuple>
#include "DenseMatrix.h"


template<typename Type>
class CSR {
private:
    std::vector<Type> values;
    std::vector<Type> column_index;
    std::vector<Type> line_index;
    std::size_t height_, width_;

public:

    CSR(const std::initializer_list<Type> &new_values,
        const std::initializer_list<Type> &new_column_index,
        const std::initializer_list<Type> &new_row_count,
        const std::size_t &height,
        const std::size_t &width) :
            values(new_values),
            column_index(new_column_index),
            line_index(new_row_count),
            height_(height),
            width_(width) {};


    CSR(const std::vector<std::tuple<std::size_t, std::size_t, Type>> &triplets, int height, int width) :
            height_(height), width_(width) {

        std::vector<std::size_t> line_counts(height_, 0);

        for (const auto &triple: triplets) line_counts[std::get<0>(triple)]++;

        line_index.resize(height_ + 1, 0);

        for (int i = 1; i <= height_; ++i) line_index[i] = line_index[i - 1] + line_counts[i - 1];

        values.resize(triplets.size());
        column_index.resize(triplets.size());
        std::vector<std::size_t> current_line(height_, 0);
        for (const auto &triplet: triplets) {
            double val = std::get<2>(triplet);
            int i = std::get<0>(triplet);
            int j = std::get<1>(triplet);
            values[line_index[i] + current_line[i]] = val;
            column_index[line_index[i] + current_line[i]] = j;
            current_line[i]++;
        }
    }


    [[nodiscard]]Type operator()(const std::size_t &i, const std::size_t &j) const {
        {

            for (std::size_t k = line_index[i]; k < line_index[i + 1]; ++k) if (column_index[k] == j) return values[k];

            return static_cast<Type>(0);
        }

    }


    [[nodiscard]] std::vector<Type> operator*(const std::vector<Type> &vector) const {
        std::vector<Type> result(height_, 0.);

        for (size_t j = 0; j < line_index.size() - 1; ++j) {
            for (size_t i = line_index[j]; i < line_index[j + 1]; ++i) {
                result[j] += values[i] * vector[column_index[i]];
            }
        }

        return result;
    }

    [[nodiscard]] std::size_t get_width() const {
        return width_;
    }

};

#endif //SLAE_CSR_H
