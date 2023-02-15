//
// Created by ivankhripunov on 10.02.23.
//

#include <gtest/gtest.h>
#include "../src/Solver/TridiagonalMatrixSolver.h"

TEST(TRIDIOGONAL_TRIO, CONSTRUCTOR) {
    Trio<double> trio_1(1.0, -999.4, 25);

    ASSERT_DOUBLE_EQ(1.0, trio_1[0]);
    ASSERT_DOUBLE_EQ(-999.4, trio_1[1]);
    ASSERT_DOUBLE_EQ(25, trio_1[2]);
}


TEST(TRIDIAGONAL, CONSTRUCTOR) {
    std::vector<double> a = {0, 1, 1, 1};
    std::vector<double> b = {2, 10, -5, 4};
    std::vector<double> c = {1, -5, 2, 0};

    TridiagonalMatrix<double> matrix_1(a, b, c);

    for (size_t i = 0; i < matrix_1.size(); ++i) {
        ASSERT_DOUBLE_EQ(a[i], matrix_1[i][0]);
        ASSERT_DOUBLE_EQ(b[i], matrix_1[i][1]);
        ASSERT_DOUBLE_EQ(c[i], matrix_1[i][2]);
    }
}

TEST(TRIDIAGONAL, SOLVE_COEFFS) {
    TridiagonalMatrix<double> matrix_2({0, 1, 1, 1},
                                       {2, 10, -5, 4},
                                       {1, -5, 2, 0});

    std::vector<double> p_true = {-0.5, 10. / 19, 38. / 85, 0.};
    std::vector<double> q_true = {-5. / 2, -31. / 19, 729. / 85, -8.};

    std::vector<double> b = {-5, -18, -40, -27};
    std::vector<double> p(b.size());
    std::vector<double> q(b.size());

    calc_coeffs(matrix_2, b, p, q);

    for (size_t i = 0; i < b.size(); ++i) {
        ASSERT_DOUBLE_EQ(p[i], p_true[i]);
        ASSERT_DOUBLE_EQ(q[i], q_true[i]);
    }
}

TEST(TRIDIAGONAL, SOLVE_1) {
    TridiagonalMatrix<double> matrix_3({0, 1, 1, 1},
                                       {2, 10, -5, 4},
                                       {1, -5, 2, 0});

    std::vector<double> b = {-5, -18, -40, -27};
    std::vector<double> result_true = {-3, 1, 5, -8};

    std::vector<double> result = solve_tridiagonal(matrix_3, b);

    for (size_t i = 0; i < matrix_3.size(); ++i) {
        ASSERT_DOUBLE_EQ(result[i], result_true[i]);
    }
}

TEST(TRIDIAGONAL, SOLVE_2) {
    TridiagonalMatrix<double> matrix_4({0, 7, 5.3, -22.3},
                                       {25, -46, 98.5, 27},
                                       {5, -6.2, 4.2, 0});

    std::vector<double> b = {24.5, 345.1, -18.2, 13.5};
    std::vector<double> result_true = {2.411698578850735,
                                       -7.158492894253677,
                                       0.172994062520050,
                                       0.642880281266560};

    std::vector<double> result = solve_tridiagonal(matrix_4, b);

    for (size_t i = 0; i < matrix_4.size(); ++i) {
        ASSERT_NEAR(result[i], result_true[i], 1e-15);
    }
}
