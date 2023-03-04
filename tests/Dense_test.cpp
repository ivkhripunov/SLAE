//
// Created by ivankhripunov on 04.03.23.
//

#include <gtest/gtest.h>
#include "../src/Matrix/DenseMatrix.h"

TEST(DENSE, CONSTRUCTOR) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    DenseMatrix A(2, 3, data);

    ASSERT_DOUBLE_EQ(A(0, 0), 1.0);
    ASSERT_DOUBLE_EQ(A(0, 1), 2.0);
    ASSERT_DOUBLE_EQ(A(0, 2), 3.0);
    ASSERT_DOUBLE_EQ(A(1, 0), 4.0);
    ASSERT_DOUBLE_EQ(A(1, 1), 5.0);
    ASSERT_DOUBLE_EQ(A(1, 2), 6.0);
}

TEST(DENSE, MULTIPLICATION) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    DenseMatrix A(2, 3, data);

    std::vector<double> x = {1.0, 2.0, 3.0};


    std::vector<double> y = A * x;

    ASSERT_EQ(y.size(), 2);

    ASSERT_DOUBLE_EQ(y[0], 14.0);
    ASSERT_DOUBLE_EQ(y[1], 32.0);
}

TEST(DENSE, ADD) {

    std::vector<double> data1 = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> data2 = {5.0, 6.0, 7.0, 8.0};
    DenseMatrix A(2, 2, data1);
    DenseMatrix B(2, 2, data2);

    DenseMatrix C = A + B;


    ASSERT_DOUBLE_EQ(C(0, 0), 6.0);
    ASSERT_DOUBLE_EQ(C(0, 1), 8.0);
    ASSERT_DOUBLE_EQ(C(1, 0), 10.0);
    ASSERT_DOUBLE_EQ(C(1, 1), 12.0);
}