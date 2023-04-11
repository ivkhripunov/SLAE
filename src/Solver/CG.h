
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

    std::ofstream fout("/home/ivankhripunov/CLionProjects/SLAE/tests/log/CG.txt");
    std::size_t counter = 0;

    std::vector<Type> result = initial_guess;
    std::vector<Type> r = A * result - b, old_r;
    std::vector<Type> d = r;
    Type alpha;

    while (second_norm(r) > tolerance) {

        fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
        counter++;

        alpha = d * r / (d * (A * d));

        result = result - alpha * d;

        old_r = r;
        r = A * result - b;

        d = r + ((r * r) / (d * old_r)) * d;
    }

    fout << " " << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;

    return result;
}

#endif //SLAE_CG_H
