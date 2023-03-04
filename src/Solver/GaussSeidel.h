//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_GAUSSSEIDEL_H
#define SLAE_GAUSSSEIDEL_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"

template<typename T>
T max(T a, T b) {
    return a > b ? a : b;
}


template<typename Type>
std::vector<Type>
GaussSeidel(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
            const Type &tolerance) {

    std::vector<Type> result = initial_guess;

    while (third_norm(A * result - b) > tolerance) {
        for (long long i = 0; i < b.size(); ++i) {
            Type sum = static_cast<Type>(0);

            for (long long j = 0; j < i - 1; ++j) sum += A(i, j) * result[j];


            for (long long j = i; j < b.size(); ++j) sum += A(i, j) * result[j];

            result[i] = (1 / A(i, i)) * (b[i] - sum);

        }
    }
    return result;
}

/*
template <typename T>
std::vector<T> GaussSeidel(const CSR<T> &A, const std::vector<T> &b, const std::vector<T> &x0, const T tolerance)
{
    std::vector<T> xk = x0;
    T sum_1 = 0, sum_2 = 0;
    while(third_norm(A * xk - b) > tolerance){
        for(long long k = 0; k < x0.size(); k ++){
            for(long long i = 0; i < k - 1; i ++)
                sum_1 += A(k, i) * xk[i];
            for(long long i = k + 1; i < x0.size(); i ++)
                sum_2 += A(k, i) * xk[i];
            xk[k] = (1 / A(k, k)) * (b[k] - sum_1 - sum_2);
            sum_1 = 0;
            sum_2 = 0;
        }
    }
    return xk;
}*/



#endif //SLAE_GAUSSSEIDEL_H
