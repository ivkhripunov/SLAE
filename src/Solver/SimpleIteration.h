//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_SIMPLEITERATION_H
#define SLAE_SIMPLEITERATION_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"

template<typename Type>
std::vector<Type>
SimpleIteration(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance,
                const Type &step) {

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (third_norm(r) > tolerance) {
        result = result - r * step;
        r = A * result - b;
    }
    return result;
}

#endif //SLAE_SIMPLEITERATION_H
