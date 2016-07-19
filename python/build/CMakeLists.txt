cmake_minimum_required(VERSION 3.0)
project(galactica_bindings)

add_subdirectory(pybind11)

set(
    BINDING_SOURCES
    "../module.cpp"
    "../image.cpp"
)

include_directories("../../include")

set(BINDING_NAME "tomography")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

pybind11_add_module(${BINDING_NAME} ${BINDING_SOURCES})
