cmake_minimum_required(VERSION 3.0)
project(galactica_bindings)

add_subdirectory(pybind11)

set(
    BINDING_SOURCES
    "../src/module.cpp"
)

include_directories("../../include" "../../ext/fmt" "../include")
link_directories("../../bin")

find_package(CUDA 7.5 QUIET)

set(
    LIB_NAMES
    "pthread"
    "galactica"
)

if(CUDA_FOUND)
    # if we found cuda, define -DUSE_CUDA and additionally link with
    # galactica-cuda
    set(LIB_NAMES ${LIB_NAMES} "galactica_cuda" "cudart")
    add_definitions(-DUSE_CUDA)
    include_directories("../../cuda/" "${CUDA_TOOLKIT_ROOT_DIR}/include")
    link_directories("${CUDA_TOOLKIT_ROOT_DIR}/lib64")
endif(CUDA_FOUND)

set(BINDING_NAME "py_galactica")
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wfatal-errors -Wextra -O3")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

pybind11_add_module(${BINDING_NAME} ${BINDING_SOURCES})


target_link_libraries(${BINDING_NAME} ${LIB_NAMES})
