# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
#include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
#add_subdirectory(Google_tests/lib)
add_subdirectory(eigen)

add_executable(threediagonal_test threediagonal_test.cpp) # здесь название исполняемого файла и файл в котором лежит тест
add_test(NAME threediagonal_test COMMAND threediagonal_test.cpp)
target_link_libraries(threediagonal_test SLAE gtest gtest_main)  # Здесь линковка, тот же исполняемый файл, название проекта и эти же либы

add_executable(CSR_test Samostoyatelnaya.cpp)
#add_test(NAME CSR_test COMMAND CSR_test.cpp)
target_link_libraries(CSR_test SLAE gtest gtest_main)

add_executable(Dense_test Dense_test.cpp)
add_test(NAME Dense_test COMMAND Dense_test.cpp)
target_link_libraries(Dense_test SLAE gtest gtest_main)

#add_executable(SimpleIteration_test Iteration_test.cpp)
#add_test(NAME SimpleIteration_test COMMAND Iteration_test.cpp)
#target_link_libraries(SimpleIteration_test SLAE gtest gtest_main)

#add_executable(Samostoyatelnaya Samostoyatelnaya.cpp)
#add_test(NAME Samostoyatelnaya COMMAND Samostoyatelnaya.cpp)
#target_link_libraries(Samostoyatelnaya SLAE gtest gtest_main)
