message(STATUS "Checking for regression tests ...")
find_package(GTest)
set(tmp FALSE)
if (GTest_FOUND)
    set(tmp TRUE)

    add_executable(hqrrp_tests test_hqrrp.cpp)
    target_link_libraries(hqrrp_tests hqrrp GTest::GTest GTest::Main)
    gtest_discover_tests(hqrrp_tests)

endif()
message(STATUS "Checking for regression tests ... ${tmp}")
