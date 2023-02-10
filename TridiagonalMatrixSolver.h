//
// Created by ivankhripunov on 10.02.23.
//

#ifndef SLAE_TRIDIAGONALMATRIXSOLVER_H
#define SLAE_TRIDIAGONALMATRIXSOLVER_H

#include "TridiagonalMatrix.h"

template<typename Type>
std::vector<Type> solve_tridiagonal(const TridiagonalMatrix<Type> &matrix, const std::vector<Type> &right_hand_column) {

    size_t size_ = matrix.size();
    std::vector<Type> result(size_);

    if (matrix.size() != right_hand_column.size()) {
        std::cout << "different sizes";
    } else {
        std::vector<Type> p_values(size_);
        std::vector<Type> q_values(size_);

        calc_coeffs(matrix, right_hand_column, p_values, q_values);

        result[size_ - 1] = q_values[size_ - 1];
        for (int i = size_ - 2; i >= 0; --i) {
            result[i] = p_values[i] * result[i + 1] + q_values[i];
        }
    }
    return result;
}

template<typename Type>
void calc_coeffs(const TridiagonalMatrix<Type> &Matrix, const std::vector<Type> &right_hand_column,
                 std::vector<Type> &p_coeffs, std::vector<Type> &q_coeffs) {

    p_coeffs[0] = -Matrix[0].third_element / Matrix[0].second_element;
    q_coeffs[0] = right_hand_column[0] / Matrix[0].second_element;

    for (size_t i = 0; i < Matrix.size() - 1; ++i) {
        p_coeffs[i + 1] =
                -Matrix[i + 1].third_element / (Matrix[i + 1].first_element * p_coeffs[i] + Matrix[i + 1].second_element);

        q_coeffs[i + 1] =
                (right_hand_column[i + 1] - Matrix[i + 1].first_element * q_coeffs[i]) /
                (Matrix[i + 1].first_element * p_coeffs[i] + Matrix[i + 1].second_element);

    }
}

#endif //SLAE_TRIDIAGONALMATRIXSOLVER_H
