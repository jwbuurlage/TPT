cmake_minimum_required(VERSION 3.0.0)


project(tomo)

link_directories("~/.local/lib")

include_directories(
    "../include"
    "../ext/hemi/"
    "../ext/fmt/"
    "../ext/date/"
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
    "../ext/fmt/fmt/format.cc"
)

add_library( ${LIB_NAME} SHARED ${LIB_SOURCES} )

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
    "boost_program_options"
)

foreach(source_file ${EXAMPLE_SOURCES})
    string(REPLACE ".cpp" "" source_name ${source_file})
    add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
    target_link_libraries( ${source_name} ${LIB_NAMES} )
endforeach(source_file)


find_package(CUDA 7.5 QUIET)

if(CUDA_FOUND)
    set(CUDA_64_BIT_DEVICE_CODE ON)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)

    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 --expt-extended-lambda")

    set(
        CUDA_LIB_SOURCES
        "../cuda/forward_project.cu"
        "../cuda/sart.cu"
        "../cuda/sirt.cu"
    )

    set(
        CUDA_EXAMPLE_SOURCES
        "../examples/cuda_reconstruction.cpp"
    )

    set(CUDA_LIB_NAME "galactica_cuda")

    cuda_add_library(${CUDA_LIB_NAME} SHARED ${CUDA_LIB_SOURCES})
    cuda_add_executable("cuda_reconstruction" ${CUDA_EXAMPLE_SOURCES})
    target_link_libraries("cuda_reconstruction" ${LIB_NAMES} ${CUDA_LIB_NAME})

endif(CUDA_FOUND)

# add_subdirectory(../test/build ../test/build)
