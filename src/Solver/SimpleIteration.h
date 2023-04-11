//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_SIMPLEITERATION_H
#define SLAE_SIMPLEITERATION_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include <cmath>
#include <fstream>


template<typename Type>
std::vector<Type>
SimpleIteration(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance,
                const Type &step) {

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/log/MPI.txt");
    std::size_t counter = 0;

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (second_norm(r) > tolerance) {

        fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
        counter++;

        result = result - r * step;
        r = A * result - b;
    }

    fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;

    fout.close();

    return result;

}


std::vector<std::size_t> calc_indexes_sequence(const std::size_t &roots_count) {

    std::vector<std::size_t> sequence(roots_count);
    sequence[0] = 1;

    for (std::size_t i = 1; i < roots_count; i *= 2) {
        for (int j = 0; j < i; ++j) {
            sequence[roots_count * (j * 2 + 1) / (2 * i)] = 2 * i + 1 - sequence[roots_count * 2 * j / (2 * i)];
        }
    }


    return sequence;
}

template<typename Type>
std::vector<Type>
cacl_sequenced_polynom_roots(const std::size_t &R, const Type &min_eigen_value, const Type &max_eigen_value) {

    std::size_t roots_count = std::pow(2, R);
    std::vector<std::size_t> sequence(calc_indexes_sequence(roots_count));

    std::vector<Type> roots(roots_count);

    const Type sin_a = std::sin(M_PI / static_cast<Type>(roots_count));
    const Type cos_a = std::cos(M_PI / static_cast<Type>(roots_count));
    Type sin_b = std::sin(M_PI / (2 * static_cast<Type>(roots_count)));
    roots[0] = std::cos(M_PI / (2 * static_cast<Type>(roots_count)));

    for (std::size_t i = 1; i < roots_count; ++i) {
        roots[i] = roots[i - 1] * cos_a - sin_a * sin_b;
        sin_b = roots[i - 1] * sin_a + sin_b * cos_a;
    }

    std::vector<Type> sequenced_roots = roots;
    for (std::size_t i = 0; i < roots_count; ++i) {
        sequenced_roots[i] = 2 / ((max_eigen_value + min_eigen_value) +
                                  (max_eigen_value - min_eigen_value) * roots[sequence[i] - 1]);
    }


    return sequenced_roots;
}

template<typename Type>
std::vector<Type>
MPI_ChebyshevAccelerationLog(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                             const std::size_t &R, const Type &min_eigen_value,
                             const Type &max_eigen_value, const Type &accuracy) {

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/log/MPI_CHEB.txt");
    std::size_t counter = 0;

    std::vector<Type> step_array = cacl_sequenced_polynom_roots(R, min_eigen_value, max_eigen_value);

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (second_norm(r) > accuracy) {

        fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
        counter++;

        result = result - r * step_array[counter % 8];
        r = A * result - b;

        counter++;
    }

    fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;

    return result;
}


template<typename Type>
std::vector<Type>
FastestGradientDescent(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                   const Type &tolerance) {

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/log/FGD.txt");
    std::size_t counter = 0;

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (second_norm(r) > tolerance) {

        fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
        counter++;

        Type step = r * r / (r * (A * r));

        result = result - r * step;
        r = A * result - b;
    }

    fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;

    fout.close();

    return result;

}
#endif //SLAE_SIMPLEITERATION_H
