//
// Created by ivankhripunov on 18.02.23.
//

#include <gtest/gtest.h>
#include "../src/Matrix/CSR.h"

TEST(CSR_MATRIX, ELEMENT_ACCESS) {
    CSR<double> matrix({1, 2, 3, 4, 1, 11}, {0, 1, 3, 2, 1, 3}, {0, 3, 4, 6}, 3, 4);

    ASSERT_DOUBLE_EQ(matrix(2, 3), 11);
    ASSERT_DOUBLE_EQ(matrix(1, 1), 0);
    ASSERT_DOUBLE_EQ(matrix(0, 0), 1);
}

TEST(CSR_MATRIX, CONSTRUCTOR) {
    std::vector<Element<double>> vector = {
            {0, 1, 2},
            {0, 0, 1},
            {0, 3, 3},
            {1, 2, 4},
            {2, 1, 1},
            {2, 3, 11}};

    CSR matrix(vector, 3, 4);
}