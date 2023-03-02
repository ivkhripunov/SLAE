cmake_minimum_required(VERSION 3.24)
project(SLAE)


set(CMAKE_CXX_STANDARD 17)
file(GLOB_RECURSE source *.hpp *.h *.cpp)
add_library(SLAE INTERFACE ${source} Matrix Solver Matrix/DenseMatrix.h)