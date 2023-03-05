//
// Created by ivankhripunov on 04.03.23.
//

#ifndef SLAE_SIMPLEITERATION_H
#define SLAE_SIMPLEITERATION_H

#include "../Matrix/CSR.h"
#include "../Utilities/Norm.h"
#include "../Utilities/Overload.h"
#include<fstream>

template<typename Type>
std::vector<Type>
SimpleIteration(const CSR<Type> &A, const std::vector<Type> &b, const std::vector<Type> &initial_guess,
                const Type &tolerance,
                const Type &step) {

    std::ofstream fout;
    fout.open("/home/ivankhripunov/CLionProjects/SLAE/tests/mpi.txt");

    long long int  counter = 0;

    std::vector<Type> r = A * initial_guess - b;
    std::vector<Type> result = initial_guess;

    while (counter < 1000) {
        result = result - r * step;
        r = A * result - b;

        counter ++;

        fout << counter << " " << third_norm(r) << std::endl;


    }

    fout.close();

    return result;


}

#endif //SLAE_SIMPLEITERATION_H
