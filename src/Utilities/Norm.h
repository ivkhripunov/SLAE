//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_NORM_H
#define SLAE_NORM_H

#include <vector>
#include <algorithm>
#include <cmath>

template<typename Type>
Type second_norm(const std::vector<Type> &vector) {

    Type norm = 0;

    for (const auto &elem : vector) norm += elem * elem;
    return pow(norm, 0.5);
    /*return std::abs(*std::max_element(vector.begin(), vector.end(), [](const Type &a, const Type &b) {
        return std::abs(a) < std::abs(b);
    }));*/
}

template<typename Type>
Type inf_norm(const std::vector<Type> &vector) {

    return std::abs(*std::max_element(vector.begin(), vector.end(), [](const Type &a, const Type &b) {
        return std::abs(a) < std::abs(b);
    }));
}

#endif //SLAE_NORM_H
