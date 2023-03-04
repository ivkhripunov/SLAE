//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_HOUSEHOLDER_H
#define SLAE_HOUSEHOLDER_H

#include "DenseMatrix.h"

template<typename Type>
std::tuple<DenseMatrix<Type>, DenseMatrix<Type>> QR_Householder(const DenseMatrix<Type> &A) {
    int width = A.get_width();
    int height = A.get_height();



}

#endif //SLAE_HOUSEHOLDER_H
