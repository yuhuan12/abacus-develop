cmake_minimum_required(VERSION 3.20)

project(CG_test)
project(TN_test)

set(CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED True)

enable_testing()

add_subdirectory(../ ./opt_binary)

add_executable(CG_test CG_test.cpp test_tools.cpp)
add_executable(TN_test TN_test.cpp test_tools.cpp)

target_link_libraries(CG_test PUBLIC opt gtest gtest_main pthread)
target_link_libraries(TN_test PUBLIC opt gtest gtest_main pthread)

include(GoogleTest)
gtest_discover_tests(CG_test)
gtest_discover_tests(TN_test)