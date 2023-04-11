#include "../src/Matrix/CSR.h"
#include <gtest/gtest.h>
#include "../src/Solver/SimpleIteration.h"

TEST(Samostoyayelnaya, Task_1) {

    std::vector<DOK<double>> triplets;

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

    //std::vector<double> result_mpi_1 = SimpleIterationLog(A, b, initial, 1e-13, 1 / lambda_max);
    //std::vector<double> result_mpi_2 = SimpleIterationLog(A, b, initial, 1e-13, 2 / (lambda_max + lambda_min));
    //std::vector<double> result_mpi_3 = MPI_ChebyshevAccelerationLog(A, b, initial, 3, lambda_min, lambda_max, 1e-13);
    std::vector<double> result_mpi_4 = A.SSOR(b, initial, 1e-13, 1.05);
    std::vector<double> result_mpi_5 = A.SOR(b, initial, 1e-13, 1.05);

}