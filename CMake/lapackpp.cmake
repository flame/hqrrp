message(STATUS "Checking for lapackpp ... ")
find_package(lapackpp REQUIRED)
message(STATUS "Checking for lapackpp ... ${lapackpp_VERSION}")

add_library(HQRRP_lapackpp INTERFACE)

target_link_libraries(HQRRP_lapackpp INTERFACE lapackpp)

install(TARGETS HQRRP_lapackpp EXPORT HQRRP_lapackpp)

install(EXPORT HQRRP_lapackpp
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
    EXPORT_LINK_INTERFACE_LIBRARIES    
)