//
// Created by ivankhripunov on 17.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <vector>
#include <algorithm>

template<typename Type>
class Element {

private:
    std::size_t i_, j_;
    Type value_;

public:
    Element(const std::size_t &new_i,
            const std::size_t &new_j,
            const Type &new_value) :
            i_(new_i),
            j_(new_j),
            value_(new_value) {};

    Element(std::initializer_list<Type> &initializer) :
            i_(initializer[0]),
            j_(initializer[1]),
            value_(initializer[2]) {};

    [[nodiscard]] std::size_t get_i() const {
        return i_;
    }

    [[nodiscard]] std::size_t get_j() const {
        return j_;
    }

    Type get_value() const {
        return value_;
    }

    bool operator<(Element<Type> const &second) const {
        return this->i_ < second.i_ or (this->i_ == second.i_ and this->j_ < second.j_);
    };
};

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

    CSR(std::vector<Element<Type>> vector_of_elements,
        const std::size_t &height,
        const std::size_t &width) : height_(height), width_(width) {

        sort(vector_of_elements.begin(), vector_of_elements.end());
        values.resize(vector_of_elements.size());
        column_index.resize(vector_of_elements.size());
        line_index.resize(height_ + 1);
        line_index[0] = 0;

        for (size_t i = 0; i < vector_of_elements.size(); ++i) {
            values[i] = vector_of_elements[i].get_value();
            column_index[i] = vector_of_elements[i].get_j();
            line_index[vector_of_elements[i].get_i() + 1] += 1;
        }

        for (size_t i = 0; i < line_index.size() - 1; ++i) {
            line_index[i + 1] += line_index[i];
        }
    }

    Type operator()(const std::size_t &i, const std::size_t &j) const {
        {

            for (std::size_t k = line_index[i]; k < line_index[i + 1]; ++k) if (column_index[k] == j) return values[k];

            return 0;
        }

    }

    std::vector<Type> operator*(const std::vector<Type> &vector) {
        std::vector<Type> result(vector.size());

    }

};

#endif //SLAE_CSR_H
