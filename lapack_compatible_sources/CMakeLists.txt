cmake_minimum_required(VERSION 3.10)

project(hqrrp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# find_package(MKL CONFIG REQUIRED)
# message(STATUS "MKL's imported targets are\n-- \t${MKL_IMPORTED_TARGETS}") 
# message(STATUS "The value of MKL::mkl_intel_ilp64 is \n-- \t ${MKL::mkl_intel_ilp64}")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(OMP_LIB "OpenMP::OpenMP_CXX")
endif()

find_package(blaspp REQUIRED)
find_package(lapackpp REQUIRED)

add_library(hqrrp NoFLA_HQRRP_WY_blk_var4.cpp)
target_include_directories(hqrrp PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(hqrrp PUBLIC blaspp lapackpp)

add_executable(simple_test simple_test.cpp)
target_include_directories(simple_test PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(simple_test PUBLIC hqrrp lapackpp blaspp)
