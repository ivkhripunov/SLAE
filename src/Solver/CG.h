
#ifndef SLAE_CG_H
#define SLAE_CG_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include <cmath>

template<typename Type>
std::vector<Type>
CG(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance) {

    std::vector<Type> result = initial_guess;
    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> d = r;
    Type d_r = d * r;
    Type alpha = d_r / (d * (A * d));

    while (second_norm(r) > tolerance) {

        result = result - alpha * d;

        r = A * result - b;

        d_r = d * r;

        if (d_r <= 1e-15) return CG(A, b, result, tolerance);

        d = r + ((r * r) / d_r) * d;
        alpha = (d * r) / (d * (A * d));
    }

    return result;
}

#endif //SLAE_CG_H
