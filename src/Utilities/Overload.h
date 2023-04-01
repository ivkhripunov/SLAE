//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_OVERLOAD_H
#define SLAE_OVERLOAD_H

#include <vector>

template<typename Type>
std::vector<Type> operator-(const std::vector<Type>& a, const std::vector<Type>& b) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator+(const std::vector<Type>& a, const std::vector<Type>& b) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const std::vector<Type>& a, const Type &c) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] * c;
    }
    return result;
}

template<typename Type>
std::vector<Type> operator/(const std::vector<Type>& a, const Type &c) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] / c;
    }
    return result;
}


#endif //SLAE_OVERLOAD_H
