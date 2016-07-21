cmake_minimum_required(VERSION 3.0)
project(galactica_bindings)

add_subdirectory(pybind11)

set(
    BINDING_SOURCES
    "../algorithms.cpp"
    "../geometry.cpp"
    "../image.cpp"
    "../module.cpp"
    "../operations.cpp"
)

include_directories("../../include" "../include")
link_directories("../../bin")

set(BINDING_NAME "py_galactica")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wfatal-errors -Wextra -O3")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

pybind11_add_module(${BINDING_NAME} ${BINDING_SOURCES})

set(
    LIB_NAMES
    "pthread"
    "galactica"
)

target_link_libraries(${BINDING_NAME} ${LIB_NAMES})
