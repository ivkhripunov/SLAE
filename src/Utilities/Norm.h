//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_NORM_H
#define SLAE_NORM_H

#include <vector>

template<typename Type>
Type inf_norm(const std::vector<Type> &vector) {

    return std::abs(*std::max_element(vector.begin(), vector.end(), [](const Type &a, const Type &b) {
        return std::abs(a) < std::abs(b);
    }));
}

#endif //SLAE_NORM_H
