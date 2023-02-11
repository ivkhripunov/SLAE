//
// Created by ivankhripunov on 10.02.23.
//

#ifndef SLAE_TRIDIAGONALMATRIX_H
#define SLAE_TRIDIAGONALMATRIX_H

#include <vector>
#include <iostream>
#include <initializer_list>

template<typename Type>
class Trio {
public:
    Type first_element;
    Type second_element;
    Type third_element;

    Trio(Type first_value, Type second_value, Type third_value) : first_element(first_value),
                                                                  second_element(second_value),
                                                                  third_element(third_value) {}

    [[nodiscard]] Type operator[](std::size_t i) const {
        switch (i) {
            case 0:
                return first_element;
            case 1:
                return second_element;
            case 2:
                return third_element;
            default:
                return;
        }

    }

};

template<typename Type>
class TridiagonalMatrix {
private:
    std::vector<Trio<Type>> matrix;
    size_t size_;

public:
    TridiagonalMatrix(const std::vector<Type> &lower_diagonal,
                      const std::vector<Type> &main_diagonal,

                      const std::vector<Type> &upper_diagonal) {

        size_ = main_diagonal.size();

        matrix.reserve(size_);

        for (std::size_t i = 0; i < size_; ++i) {
            Trio tmp(lower_diagonal[i], main_diagonal[i], upper_diagonal[i]);
            matrix.push_back(tmp);
        }
    }

    TridiagonalMatrix(std::initializer_list<Type> &init_lower_diagonal,
                      std::initializer_list<Type> &init_main_diagonal,
                      std::initializer_list<Type> &init_upper_diagonal) {

        size_ = init_main_diagonal.size();

        matrix.reserve(size_);

        for (std::size_t i = 0; i < size_; ++i) {
            Trio tmp(*(init_lower_diagonal.begin() + i),
                     *(init_main_diagonal.begin() + i),
                     *(init_upper_diagonal.begin() + i));
            matrix[i] = tmp;
        }

    }

    [[nodiscard]] size_t size() const {
        return size_;
    }

    [[nodiscard]] Trio<Type> &operator[](std::size_t i) {
        return matrix[i];
    }

    [[nodiscard]] Trio<Type> operator[](std::size_t i) const {
        return matrix[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const TridiagonalMatrix<Type> &Matrix) {
        for (size_t i = 0; i < Matrix.size_; ++i) {
            std::cout << Matrix.matrix[i].first_element << " "
                      << Matrix.matrix[i].second_element << " "
                      << Matrix.matrix[i].third_element << std::endl;
        }

        return os;
    }
};

#endif //SLAE_TRIDIAGONALMATRIX_H
