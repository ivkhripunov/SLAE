//
// Created by ivankhripunov on 04.03.23.
//

#include <gtest/gtest.h>
#include "../src/Solver/SimpleIteration.h"
#include "../src/Solver/Jacobi.h"
#include "../src/Solver/GaussSeidel.h"
/*
TEST(SIMPLEITERATION, SOLVE_1) {
    CSR<double> matrix({1, 2, 3, 4, 1, 11}, {0, 1, 3, 2, 1, 3}, {0, 3, 4, 6}, 3, 4);

    CSR<double> A = {{2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0}, {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9}, 3,
                     3};

    std::vector<double> b = {1.0, 0.0, 1.0};

    std::vector<double> initial = {0, 0, 0};

    std::vector<double> expected = {1, 1, 1};

    std::vector<double> result = SimpleIteration(A, b, initial, 1e-3, 1e-4);

    for (std::size_t i = 0; i < expected.size(); ++i) ASSERT_NEAR(expected[i], result[i], 5e-2);
}

TEST(JACOBI, SOLVE_1) {
    CSR<double> matrix({1, 2, 3, 4, 1, 11}, {0, 1, 3, 2, 1, 3}, {0, 3, 4, 6}, 3, 4);

    CSR<double> A = {{2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0}, {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9}, 3,
                     3};

    std::vector<double> b = {1.0, 0.0, 1.0};

    std::vector<double> initial = {0, 0, 0};

    std::vector<double> expected = {1, 1, 1};

    std::vector<double> result = Jacobi(A, b, initial, 1e-10);

    for (std::size_t i = 0; i < expected.size(); ++i) ASSERT_NEAR(expected[i], result[i], 1e-5);
}

TEST(GAUSSSEIDEL, SOLVE_1) {
    CSR<double> matrix({1, 2, 3, 4, 1, 11}, {0, 1, 3, 2, 1, 3}, {0, 3, 4, 6}, 3, 4);

    CSR<double> A = {{2.0,  -1.0, 0.0, -1.0, 2.0,  -1.0, 0.0,  -1.0, 2.0}, {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9}, 3, 3};

    std::vector<double> b = {1.0, 0.0, 1.0};

    std::vector<double> initial = {0, 0, 0};

    std::vector<double> expected = {1, 1, 1};

    std::vector<double> result = GaussSeidel(A, b, initial, 1e-1);

    for (std::size_t i = 0; i < expected.size(); ++i) ASSERT_NEAR(expected[i], result[i], 1e-2);
}*/

TEST(MPI, TASK_3) {
    CSR<double> A = {{10, 1, 0, 1, 7, 0, 0, 0.1, 1}, {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9}, 3, 3};
    std::vector<double> x0 = {0, 0, 0};
    std::vector<double> b = {20, 30, 1};

    std::vector<double> result = b;

    result = SimpleIteration(A, b, x0, 1e-12, 0.1935);

    //for (const auto &element: result) std::cout << element << " ";
}