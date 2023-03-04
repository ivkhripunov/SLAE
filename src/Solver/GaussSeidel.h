//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_GAUSSSEIDEL_H
#define SLAE_GAUSSSEIDEL_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"

template<typename T>
T max(T a, T b) {
    return a > b ? a : b;
}


template<typename Type>
std::vector<Type>
GaussSeidel(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
            const Type &tolerance) {

    std::vector<Type> result = initial_guess;

    while (third_norm(A * result - b) > tolerance) {
        std::cout << third_norm(A * result - b) << std::endl;
        for (long long i = 0; i < b.size(); ++i) {
            Type sum = static_cast<Type>(0);

            for (long long j = 0; j < i - 1; ++j) sum += A(i, j) * result[j];


            for (long long j = i; j < b.size(); ++j) sum += A(i, j) * result[j];

            result[i] = (1 / A(i, i)) * (b[i] - sum);

        }
    }
    return result;
}
#endif //SLAE_GAUSSSEIDEL_H
