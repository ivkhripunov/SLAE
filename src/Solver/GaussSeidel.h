//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_GAUSSSEIDEL_H
#define SLAE_GAUSSSEIDEL_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include <fstream>

template<typename T>
T max(T a, T b) {
    return a > b ? a : b;
}


template<typename Type>
std::vector<Type>
GaussSeidel(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
            const Type &tolerance) {

    std::ofstream fout;
    fout.open("/home/ivankhripunov/CLionProjects/SLAE/tests/gaussseidel.txt");
    unsigned long long counter = 0;

    std::vector<Type> result = initial_guess;

    while (counter < 1000) {

        for (long long i = 0; i < b.size(); ++i) {
            Type sum = static_cast<Type>(0);

            for (long long j = 0; j < i; ++j) sum += A(i, j) * result[j];


            for (long long j = i + 1; j < b.size(); ++j) sum += A(i, j) * result[j];

            result[i] = (1 / A(i, i)) * (b[i] - sum);

        }

        fout << counter << " " << third_norm(A * result - b) << std::endl;

        counter++;
    }
    return result;
}
#endif //SLAE_GAUSSSEIDEL_H
