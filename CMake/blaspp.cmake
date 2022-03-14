message(STATUS "Checking for blaspp ... ")
find_package(blaspp REQUIRED)
message(STATUS "Checking for blaspp ... ${blaspp_VERSION}")

add_library(HQRRP_blaspp INTERFACE)

target_link_libraries(HQRRP_blaspp INTERFACE blaspp)

install(TARGETS HQRRP_blaspp EXPORT HQRRP_blaspp)

install(EXPORT HQRRP_blaspp
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
    EXPORT_LINK_INTERFACE_LIBRARIES    
)