cmake_minimum_required(VERSION 3.0.0)

project(tomo)

include_directories(
    "../include"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wfatal-errors -Wextra -O3 -g")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(LIB_NAME "galactica")

set(
    LIB_SOURCES
    "../src/phantoms.cpp"
)

add_library( ${LIB_NAME} STATIC ${LIB_SOURCES} )

set(
    EXAMPLE_DIR
    "../examples/"
)

set(
    EXAMPLE_SOURCES
    "reconstruction.cpp"
)

set(
    LIB_NAMES
    "pthread"
    "galactica"
)

foreach(source_file ${EXAMPLE_SOURCES})
    string(REPLACE ".cpp" "" source_name ${source_file})
    add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
    target_link_libraries( ${source_name} ${LIB_NAMES} )
endforeach(source_file)

# add_subdirectory(../test/build ../test/build)
