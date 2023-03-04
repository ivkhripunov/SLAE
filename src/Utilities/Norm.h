//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_NORM_H
#define SLAE_NORM_H

#include <vector>

template<typename Type>
Type third_norm(const std::vector<Type> &vector) {
    Type result = static_cast<Type>(0);

    for (const Type &element: vector) result += element * element;

    return result;
}

#endif //SLAE_NORM_H
