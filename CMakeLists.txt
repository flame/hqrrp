cmake_minimum_required(VERSION 3.10)

project(hqrrp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# from HQRRP-build:
#   cmake -DCMAKE_INSTALL_PREFIX=`pwd`/../HQRRP-install -Dblaspp_DIR=/opt/mklpp/lib/blaspp/ -Dlapackpp_DIR=/opt/mklpp/lib/lapackpp/ -DCMAKE_BINARY_DIR=`pwd` -DCMAKE_BUILD_TYPE=Release ../lapack_compatible_sources/

# find_package(MKL CONFIG REQUIRED)
# message(STATUS "MKL's imported targets are\n-- \t${MKL_IMPORTED_TARGETS}") 
# message(STATUS "The value of MKL::mkl_intel_ilp64 is \n-- \t ${MKL::mkl_intel_ilp64}")

enable_testing()


find_package(blaspp REQUIRED)
find_package(lapackpp REQUIRED)

set(
    hqrrp_cxx
    ${PROJECT_SOURCE_DIR}/src/hqrrp.cpp
    ${PROJECT_SOURCE_DIR}/src/util.cpp
)
add_library(
    hqrrp ${hqrrp_cxx}
)
target_include_directories(
    hqrrp PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${Random123_DIR}
)
target_link_libraries(hqrrp PUBLIC blaspp lapackpp)


set(
    simple_test_cxx
    ${PROJECT_SOURCE_DIR}/src/simple_test.cpp
)
add_executable(
    simple_test ${simple_test_cxx}
)
target_include_directories(
    simple_test PUBLIC ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(
    simple_test PUBLIC hqrrp
)

add_subdirectory(test)


# I needed to build blaspp, lapackpp, and this library with
#   LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:/usr/bin
#   LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:/usr/bin
# Building blaspp should use -Dblas=mkl and -Dblas_int=int64.
#
#   ... seems this library can still compile and run when blaspp and lapackpp
#   where compiled with LD_LIBRARY_PATH=/home/riley/intel/oneapi/mkl/2022.0.2/lib/intel64
#   and LIBRARY_PATH=$LD_LIBRARY_PATH
#
# *** OR ***
#      
#   Consider when all are compiled with
#       LD_LIBRARY_PATH=/home/riley/intel/oneapi/mkl/2022.0.2/lib/intel64
#       LIBRARY_PATH=$LD_LIBRARY_PATH
#   yeah, that works too.
#
# It might be the case that the value of LD_LIBRARY_PATH doesn't affect this
# specific project (so long as it doesn't link explicitly to MKL).
