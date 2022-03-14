cmake_minimum_required(VERSION 3.10)

project(simple_test)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message(STATUS "Checking for hqrrp ... ")
find_package(hqrrp REQUIRED)
message(STATUS "Done checking for hqrrp. ...")

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
    simple_test PUBLIC hqrrp blaspp
)
