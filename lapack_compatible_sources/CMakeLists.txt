cmake_minimum_required(VERSION 3.10)

project(hqrrp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# from HQRRP-build:
#   cmake -DCMAKE_INSTALL_PREFIX=`pwd`/../HQRRP-install -Dblaspp_DIR=/opt/mklpp/lib/blaspp/ -Dlapackpp_DIR=/opt/mklpp/lib/lapackpp/ -DCMAKE_BINARY_DIR=`pwd` -DCMAKE_BUILD_TYPE=Release ../lapack_compatible_sources/

# find_package(MKL CONFIG REQUIRED)
# message(STATUS "MKL's imported targets are\n-- \t${MKL_IMPORTED_TARGETS}") 
# message(STATUS "The value of MKL::mkl_intel_ilp64 is \n-- \t ${MKL::mkl_intel_ilp64}")

find_package(blaspp REQUIRED)
find_package(lapackpp REQUIRED)

add_library(hqrrp NoFLA_HQRRP_WY_blk_var4.cpp)
target_include_directories(hqrrp PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(hqrrp PUBLIC blaspp lapackpp)

add_executable(simple_test simple_test.cpp)
target_include_directories(simple_test PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(simple_test PUBLIC hqrrp)

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
