cmake_minimum_required(VERSION 3.0.0)
project(tomo)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

link_directories("~/.local/lib")

include_directories(
    "../include"
    "../ext/hemi/"
    "../ext/fmt/"
    "../ext/date/"
    "../ext/bulk/include"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wno-missing-braces -Wfatal-errors -Werror -pedantic -O3 -g")
# ignore warnings from external libs
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ../ext/fmt -isystem ~/.local/include -isystem ../ext/bulk/include -isystem /usr/local/cuda/include")

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
    "geometry_benchmark.cpp"
)

set(
    LIB_NAMES
    "pthread"
    "galactica"
    "boost_program_options"
)

set(CMAKE_CXX_COMPILER "/ufs/buurlage/.local/bin/clang++")
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

find_package(MPI QUIET)

if(MPI_CXX_FOUND)
    set(
        MPI_DIR
        "../distributed/"
        )

    set(
        MPI_SOURCES
        "sirt_test.cpp"
        )

    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})

    # make distributed example
    foreach(source_file ${MPI_SOURCES})
        string(REPLACE ".cpp" "" source_name ${source_file})
        add_executable(${source_name} ${MPI_DIR}${source_file})
        target_link_libraries(${source_name} ${LIB_NAMES} ${MPI_LIBRARIES})
    endforeach(source_file)
else(MPI_CXX_FOUND)
    message("No MPI installation found.")
endif(MPI_CXX_FOUND)

add_subdirectory(../test/build ../test/build)
