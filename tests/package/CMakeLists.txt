# Consumer project used by tests/package/check_package.cmake:
# it must find an installed aerobus through find_package(Aerobus).
cmake_minimum_required(VERSION 3.14...3.31)
project(aerobus_package_consumer LANGUAGES CXX)

find_package(Aerobus REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE aerobus::aerobus)
set_target_properties(consumer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin")
