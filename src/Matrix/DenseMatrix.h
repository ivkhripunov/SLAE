//
// Created by ivankhripunov on 02.03.23.
//

#ifndef SLAE_DENSEMATRIX_H
#define SLAE_DENSEMATRIX_H

#include <vector>

template<typename Type>
class DenseMatrix {
private:
    std::vector<Type> matrix_;
    std::size_t height_, width_;

public:

    DenseMatrix(const std::size_t &h, const std::size_t &w) : height_(h), width_(w),
                                                             matrix_(std::vector<Type>(h * w)) {};


};

#endif //SLAE_DENSEMATRIX_H
