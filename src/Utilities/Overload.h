//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_OVERLOAD_H
#define SLAE_OVERLOAD_H

#include <ostream>
#include <vector>

template<typename Type>
std::vector<Type> operator-(const std::vector<Type> &a, const std::vector<Type> &b) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator+(const std::vector<Type> &a, const std::vector<Type> &b) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const std::vector<Type> &a, const Type &c) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] * c;
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const Type &c, const std::vector<Type> &a) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] * c;
    }
    return result;
}

template<typename Type>
Type operator*(const std::vector<Type> &a, const std::vector<Type> &b) {
    Type result = static_cast<Type>(0);

    for (std::size_t i = 0; i < a.size(); ++i) result += a[i] * b[i];

    return result;
}

template<typename Type>
std::vector<Type> operator/(const std::vector<Type> &a, const Type &c) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] / c;
    }
    return result;
}
template<typename Type>
std::ostream& operator<<(std::ostream& os, const std::vector<Type>& vec) {
    os << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}



#endif //SLAE_OVERLOAD_H
