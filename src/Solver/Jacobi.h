//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_JACOBI_H
#define SLAE_JACOBI_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include <fstream>

template<typename Type>
std::vector<Type>
Jacobi(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance) {

    std::ofstream fout;
    fout.open("/home/ivankhripunov/CLionProjects/SLAE/tests/jacobi.txt");
    unsigned long long counter = 0;

    std::vector<Type> result = initial_guess;
    std::vector<Type> tmp = initial_guess;

    while (counter < 1000) {
        tmp = result;
        for (int i = 0; i < b.size(); ++i) {
            Type sum = static_cast<Type>(0);

            for (std::size_t j = 0; j < A.get_width(); ++j) {
                if (i != j) sum += A(i, j) * tmp[j];
            }

            result[i] = (b[i] - sum) / A(i, i);

        }

        fout << counter << " " << third_norm(A * result - b) << std::endl;

        counter++;
    }
    return result;

}

#endif //SLAE_JACOBI_H
