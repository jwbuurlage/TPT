cmake_minimum_required(VERSION 3.0)
project(tomos_bindings)

add_subdirectory(pybind11)

set(
    BINDING_SOURCES
    "../src/module.cpp"
)

set(BINDING_NAME "py_tomos")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../python/tomos)

pybind11_add_module(${BINDING_NAME} ${BINDING_SOURCES})

target_include_directories(${BINDING_NAME} PUBLIC ../include)
target_link_libraries(${BINDING_NAME} PUBLIC tomos)
