//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_SIMPLEITERATION_H
#define SLAE_SIMPLEITERATION_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include <cmath>

template<typename Type>
std::vector<Type>
SimpleIteration(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance,
                const Type &step) {

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (inf_norm(r) > tolerance) {
        result = result - r * step;
        r = A * result - b;
    }

    return result;

}

template<typename Type>
std::vector<Type> cacl_polynom_roots(const std::size_t &R, const Type &min_eigen_value, const Type &max_eigen_value) {

    std::size_t roots_count = std::pow(2, R);

    std::vector<Type> roots(roots_count);

    const Type sin_a = std::sin(M_PI / roots_count);
    const Type cos_a = std::cos(M_PI / roots_count);
    Type sin_b = std::sin(M_PI / (2 * roots_count));
    roots[0] = std::cos(M_PI / (2 * roots_count));

    for (std::size_t i = 1; i < roots_count; ++i) {
        roots[i] = roots[i - 1] * cos_a - sin_a * sin_b;
        sin_b = roots[i - 1] * sin_a + sin_b * cos_a;
        roots[i - 1] = (max_eigen_value + min_eigen_value) / 2 + (max_eigen_value - min_eigen_value) / 2 * roots[i - 1];
    }

    roots[roots_count - 1] =
            (max_eigen_value + min_eigen_value) / 2 + (max_eigen_value - min_eigen_value) / 2 * roots[roots_count - 1];


    return roots;
}

std::vector<std::size_t> calc_indexes_sequence(const std::size_t &R) {

    std::size_t roots_count = std::pow(2, R);
    std::vector<std::size_t> sequence(roots_count);

    for (std::size_t i = 2; i <= R; ++i) {
        int step = static_cast<int>(pow(2, R - i));
        for (int j = 0; j < roots_count; j += 2 * step) {
            sequence[j + step] = static_cast<size_t>(pow(2, i)) - sequence[j] - 1;
        }
    }

    return sequence;
}

template<typename Type>
std::vector<Type>
MPI_ChebyshevAcceleration(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                          const std::size_t &R, const Type &min_eigen_value,
                          const Type &max_eigen_value) {

    std::vector<Type> step_array(cacl_polynom_roots(R, min_eigen_value, max_eigen_value));
    std::vector<std::size_t> sequence(calc_indexes_sequence(R));

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    for (const std::size_t &i: sequence) {

        result = result - r / step_array[i];
        r = A * result - b;

    }

    return result;
}

#endif //SLAE_SIMPLEITERATION_H
