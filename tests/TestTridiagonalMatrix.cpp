//
// Created by ivankhripunov on 10.02.23.
//

#include "../src/Solvers/TridiagonalMatrixSolver.h"


int main() {
    std::vector<double> a = {666, 666, 666, 666, 666};
    std::vector<double> b = {666, 666, 666, 666, 666};
    std::vector<double> c = {11, 12, 13, 14, 0};
    std::vector<double> d = {-5, -18, -40, -27};

    TridiagonalMatrix<double> matrix_1(a, b, c);
    TridiagonalMatrix<double> matrix_2({0, 1, 1, 1}, {2, 10, -5, 4}, {1, -5, 2, 0});

    std::cout << matrix_1;

    calc_coeffs(matrix_2, d, a, b);
    std::vector<double> result = solve_tridiagonal(matrix_2, d);

    for (const auto x: result) std::cout << x << " ";
    std::cout << std::endl;

    for (const auto x: a) std::cout << x << " ";
    std::cout << std::endl;
    for (const auto x: b) std::cout << x << " ";
}
