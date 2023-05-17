#include "../src/Matrix/CSR.h"
#include <gtest/gtest.h>
#include "../src/Solver/SimpleIteration.h"
#include "../src/Solver/CG.h"

TEST(Samostoyayelnaya, Task_1) {

    std::vector<Triplet<double>> triplets;

    double a_ = 11, b_ = 25;
    double lambda_max = 2 * (b_ + 2 * a_ * cos(M_PI / 18));
    double lambda_min = 2 * (b_ - 2 * a_ * cos(M_PI / 18));

    for (std::size_t i = 0; i < 289; ++i) {
        triplets.push_back({i, i, 2 * b_});
    }

    for (std::size_t i = 0; i < 288; ++i) {
        triplets.push_back({i + 1, i, a_});
        triplets.push_back({i, i + 1, a_});
    }

    for (std::size_t i = 0; i < 272; ++i) {
        triplets.push_back({i + 17, i, a_});
        triplets.push_back({i, i + 17, a_});
    }

    CSR<double> A{triplets, 289, 289};

    //std::cout << lambda_min << " " << lambda_max;

    std::vector<double> b(289, 1);

    std::vector<double> initial(289, 0);

    //std::vector<double> result_mpi_1 = SimpleIteration(A, b, initial, 1e-13, 1 / lambda_max);
    //std::vector<double> result_mpi_2 = SimpleIteration(A, b, initial, 1e-13, 2 / (lambda_max + lambda_min));
    std::vector<double> result_mpi_3 = MPI_ChebyshevAccelerationLog(A, b, initial, 3, lambda_min, lambda_max, 1e-13);
    //std::vector<double> result_mpi_4 = A.SSOR(b, initial, 1e-13, 1.05);
    //std::vector<double> result_mpi_5 = A.SOR(b, initial, 1e-13, 1.05);

}

TEST(Samostoyayelnaya, Task_1_1) {

    std::vector<Triplet<double>> triplets;

    double a_ = 11, b_ = 25;
    double lambda_max = 2 * (b_ + 2 * a_ * cos(M_PI / 18));
    double lambda_min = 2 * (b_ - 2 * a_ * cos(M_PI / 18));

    for (std::size_t i = 0; i < 289; ++i) {
        triplets.push_back({i, i, 2 * b_});
    }

    for (std::size_t i = 0; i < 288; ++i) {
        triplets.push_back({i + 1, i, a_});
        triplets.push_back({i, i + 1, a_});
    }

    for (std::size_t i = 0; i < 272; ++i) {
        triplets.push_back({i + 17, i, a_});
        triplets.push_back({i, i + 17, a_});
    }

    CSR<double> A{triplets, 289, 289};

    //std::cout << lambda_min << " " << lambda_max;

    std::vector<double> b(289, 1);

    std::vector<double> initial(289, 0);

    for (int i = -1000; i <= 1000; ++i) {
        double step = 5e-3;
        std::vector<double> result_mpi_3 = MPI_ChebyshevAccelerationLog(A, b, initial, 3, lambda_min - i * step,
                                                                        lambda_max + i * step, 1e-13);
    }

}

TEST(Samostoyayelnaya, Task_2) {

    double lambda_max = 27;
    double lambda_min = 18;

    CSR<double> A{{{0, 0, 18}, {1, 1, 21}, {2, 2, 24}, {3, 3, 27}}, 4, 4};

    std::vector<double> b(4, 1);

    std::vector<double> initial(4, 0);

    std::vector<double> result_mpi_1 = SimpleIteration(A, b, initial, 1e-13, 0.9 * 2 / lambda_max);
    std::vector<double> result_mpi_2 = SimpleIteration(A, b, initial, 1e-13, 2 / (lambda_max + lambda_min));
    std::vector<double> result_mpi_3 = FastestGradientDescent(A, b, initial, 1e-13);
    std::vector<double> result_mpi_4 = MPI_ChebyshevAccelerationLog(A, b, initial, 3, lambda_min, lambda_max, 1e-13);
    std::vector<double> result_mpi_5 = CG(A, b, initial, 1e-13);
}

TEST(SLAE_EXAM, Task_4) {
    CSR<double> A{{{0, 0, 10},
                  {1, 0, 3},
                  {2, 0, 6},
                  {0, 1, 3},
                  {1, 1, 5},
                  {2, 1, 1},
                  {0, 2, 6},
                  {1, 2, 1},
                  {2, 2, 8}}, 3, 3};

    std::cout << cacl_sequenced_polynom_roots(3, 2.251314271566045, 15.882409871887806);
}

TEST(GRADIENT, HeavyBall) {
    CSR<double> A{{{0, 0, 18}, {1, 1, 21}, {2, 2, 24}, {3, 3, 27}}, 4, 4};

    std::vector<double> b(4, 1);

    std::vector<double> initial(4, 1);

    std::vector<double> result_heavy_ball = HeavyBall(A, b, initial, 1e-13);

    std::vector<double> reference_result = FastestGradientDescent(A, b, initial, 1e-13);

    ASSERT_NEAR(reference_result[0], result_heavy_ball[0], 1e-13);
    ASSERT_NEAR(reference_result[1], result_heavy_ball[1], 1e-13);
    ASSERT_NEAR(reference_result[2], result_heavy_ball[2], 1e-13);
    ASSERT_NEAR(reference_result[3], result_heavy_ball[3], 1e-13);
}