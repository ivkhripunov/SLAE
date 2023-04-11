//
// Created by ivankhripunov on 17.02.23.
//

#ifndef SLAE_CSR_H
#define SLAE_CSR_H

#include <vector>
#include <tuple>
#include <fstream>
#include "DenseMatrix.h"

#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"


template<typename Type>
struct Triplet {
    std::size_t i;
    std::size_t j;
    Type value_;

    bool operator<(const Triplet<Type> &triplet) const {
        if (i != triplet.i) return i < triplet.i;
        else return j < triplet.j;
    }

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

    CSR(const std::vector<Triplet<Type>> &triplets, const std::size_t &height, const std::size_t &width) : height_(height), width_(width) {

        std::vector<Triplet<Type>> triplets_vector = triplets;

        std::sort(triplets_vector.begin(), triplets_vector.end());

        line_index.reserve(height_ + 1);
        values.reserve(triplets_vector.size());
        column_index.reserve(triplets_vector.size());
        line_index.push_back(0);

        std::size_t element_counter = 0;
        std::size_t line_counter = 0;
        for (const Triplet<Type> &triplet: triplets_vector) {
            while (line_counter < triplet.i) {
                line_index.push_back(element_counter);
                line_counter++;
            }

            if (triplet.i == line_counter) {
                values.push_back(triplet.value_);
                column_index.push_back(triplet.j);
                element_counter++;
            }
        }
        line_index.push_back(element_counter);
    }


    CSR(const std::initializer_list<Triplet<Type>> &triplets, const std::size_t &height, const std::size_t &width) : height_(height), width_(width) {

        std::vector<Triplet<Type>> triplets_vector(triplets);

        std::sort(triplets_vector.begin(), triplets_vector.end());

        line_index.reserve(height_ + 1);
        values.reserve(triplets_vector.size());
        column_index.reserve(triplets_vector.size());
        line_index.push_back(0);

        std::size_t element_counter = 0;
        std::size_t line_counter = 0;
        for (const Triplet<Type> &triplet: triplets_vector) {
            while (line_counter < triplet.i) {
                line_index.push_back(element_counter);
                line_counter++;
            }

            if (triplet.i == line_counter) {
                values.push_back(triplet.value_);
                column_index.push_back(triplet.j);
                element_counter++;
            }
        }
        line_index.push_back(element_counter);
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

    std::vector<Type>
    Jacobi(const std::vector<Type> &b, const std::vector<Type> &initial_guess,
           const Type &tolerance);

    std::vector<Type>
    GaussSeidel(const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance);


    std::vector<Type>
    SOR(const std::vector<Type> &b, const std::vector<Type> &initial_guess,
        const Type &tolerance, const Type &w);

    std::vector<Type>
    SSOR(const std::vector<Type> &b, const std::vector<Type> &initial_guess,
        const Type &tolerance, const Type &w);


};

template<typename Type>
std::vector<Type>
CSR<Type>::Jacobi(const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance) {

    std::vector<Type> tmp(b.size()), result = initial_guess;
    Type diagonal_element;

    while (second_norm((*this) * result - b) > tolerance) {

        for (std::size_t k = 0; k < b.size(); ++k) {
            tmp[k] = b[k];
            for (size_t i = line_index[k]; i < line_index[k + 1]; ++i) {
                if (column_index[i] == k) diagonal_element = values[i];
                else tmp[k] -= values[i] * result[column_index[i]];
            }

            tmp[k] /= diagonal_element;
        }
        result = tmp;

    }

    return result;

}

template<typename Type>
std::vector<Type>
CSR<Type>::GaussSeidel(const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance) {

    std::vector<Type> tmp = initial_guess, result(b.size());
    Type diagonal_element;

    while (second_norm((*this) * result - b) > tolerance) {

        for (std::size_t k = 0; k < b.size(); ++k) {
            result[k] = b[k];

            for (size_t i = line_index[k]; i < line_index[k + 1]; ++i) {
                if (column_index[i] < k) result[k] -= values[i] * result[column_index[i]];
                else if (column_index[i] == k) diagonal_element = values[i];
                else result[k] -= values[i] * tmp[column_index[i]];
            }

            result[k] /= diagonal_element;
        }

        tmp = result;

    }

    return result;
}

template<typename Type>
std::vector<Type>
CSR<Type>::SOR(const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance,
               const Type &w) {

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/SOR.txt");
    std::size_t counter = 0;

    std::vector<Type> result = initial_guess;
    Type diagonal_element;

    while (second_norm((*this) * result - b) > tolerance) {

        fout << counter << " " << second_norm((*this) * result - b) << std::endl;
        counter++;

        for (std::size_t i = 0; i < b.size(); ++i) {
            Type sum = b[i];

            for (size_t j = line_index[i]; j < line_index[i + 1]; ++j)
                if (column_index[j] != i) sum -= values[j] * result[column_index[j]];
                else diagonal_element = values[j];

            result[i] = (1 - w) * result[i] + w * sum / diagonal_element;

        }

    }

    return result;
}

template<typename Type>
std::vector<Type>
CSR<Type>::SSOR(const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance,
               const Type &w) {

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/SSOR.txt");
    std::size_t counter = 0;

    std::vector<Type> result = initial_guess;
    Type diagonal_element;

    while (second_norm((*this) * result - b) > tolerance) {

        fout << counter << " " << second_norm((*this) * result - b) << std::endl;
        counter++;

        for (std::size_t i = 0; i < b.size(); ++i) {
            Type sum = b[i];

            for (size_t j = line_index[i]; j < line_index[i + 1]; ++j) {
                if (column_index[j] != i) sum -= values[j] * result[column_index[j]];
                else diagonal_element = values[j];
            }

            result[i] = (1 - w) * result[i] + w * sum / diagonal_element;

        }

        fout << counter << " " << second_norm((*this) * result - b) << std::endl;
        counter++;

        for (int i = b.size() - 1; i > 0; i--) {
            Type sum = b[i];

            for (size_t j = line_index[i]; j < line_index[i + 1]; ++j) {
                if (column_index[j] != i) sum -= values[j] * result[column_index[j]];
                else diagonal_element = values[j];
            }

            result[i] = (1 - w) * result[i] + w * sum / diagonal_element;

        }

    }

    return result;
}

#endif //SLAE_CSR_H
