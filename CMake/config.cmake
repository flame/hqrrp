configure_file(CMake/HQRRPConfig.cmake.in
    ${CMAKE_INSTALL_LIBDIR}/cmake/HQRRPConfig.cmake @ONLY)

install(FILES
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/cmake/HQRRPConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
