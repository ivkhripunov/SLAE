//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_JACOBI_H
#define SLAE_JACOBI_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"

template<typename Type>
std::vector<Type>
Jacobi(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess, const Type &tolerance) {

    std::vector<Type> result = initial_guess;
    std::vector<Type> tmp = initial_guess;

    while (third_norm(A * result - b) > tolerance) {
        tmp = result;
        for (int i = 0; i < b.size(); ++i) {
            Type sum = static_cast<Type>(0);
            int skip = 0;
            int count = 0;
            for (int k = skip; k < skip + count; ++k) {
                if (i != A.cols_[k]) sum += A.values_[k] * tmp[i];
                else continue;
            }
            result[i] = (b[i] - sum) / A(i, i);
        }
    }
    return result;

}

#endif //SLAE_JACOBI_H
