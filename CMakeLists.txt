cmake_minimum_required(VERSION 3.0.0)
project(tomos)

include(FeatureSummary)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(FOR_RELEASE "Build with all optimizations turned on")
option(WITH_MPI "Build MPI specific parts")
option(WITH_ZMQ "Build ZeroMQ support for visualization servers")

feature_summary(WHAT ALL)

link_directories("~/.local/lib")

include_directories(
    "include"
    "ext/fmt/include/"
    "ext/cpptoml/include/"
    "ext/bulk/include"
    "../glm"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Wno-missing-braces -Wfatal-errors -Werror -pedantic -O3 -march=native")
if(FOR_RELEASE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
else(FOR_RELEASE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(FOR_RELEASE)

# ignore warnings from external libs
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ../ext/fmt/include -DFMT_HEADER_ONLY -isystem ../ext/cpptoml/include -isystem ~/.local/include -isystem ../ext/bulk/include -isystem /usr/local/cuda/include")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(LIB_NAME "tomos")

set(
    LIB_SOURCES
    "src/phantoms.cpp"
)

add_library( ${LIB_NAME} SHARED ${LIB_SOURCES} )

set(
    EXAMPLE_DIR
    "examples/"
)

set(
    EXAMPLE_SOURCES
    "minimal_example.cpp"
    "generate_mtx.cpp"
    "reconstruction.cpp"
    "read_tiff.cpp"
    "dose_deposition.cpp"
    "geometry_benchmark.cpp"
    "partitioning_tests.cpp"
)

set(
    LIB_NAMES
    "pthread"
    "tomos"
    "boost_program_options"
    "stdc++fs"
    "opencv_highgui"
)

foreach(source_file ${EXAMPLE_SOURCES})
    string(REPLACE ".cpp" "" source_name ${source_file})
    add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
    target_link_libraries( ${source_name} ${LIB_NAMES} )
endforeach(source_file)

if(WITH_ZMQ)
    link_directories("ext/libzmq/build/lib")
    include_directories(
        "ext/tomopackets/include"
        "ext/cppzmq"
        "ext/libzmq/include"
    )
    add_executable("external_plotter" "examples/external_plotter.cpp")
    target_link_libraries("external_plotter" ${LIB_NAMES})

    add_executable("read_geometry_test" "examples/read_geometry_test.cpp")
    target_link_libraries("read_geometry_test" ${LIB_NAMES})

    #    add_executable("read_legoscan_test" "examples/read_legoscan_test.cpp")
    #    target_link_libraries("read_legoscan_test" ${LIB_NAMES})

    add_executable("reserialize" "examples/reserialize.cpp")
    target_link_libraries("reserialize" ${LIB_NAMES})

    add_executable("partition" "examples/partition.cpp")
    target_link_libraries("partition" ${LIB_NAMES})

    add_executable("comvol" "examples/generate_comvol_table.cpp")
    target_link_libraries("comvol" ${LIB_NAMES})

    add_executable("preview" "examples/preview.cpp")
    target_link_libraries("preview" ${LIB_NAMES})
endif(WITH_ZMQ)

#find_package(CUDA 7.5 QUIET)
#
#if(CUDA_FOUND)
#    set(CUDA_64_BIT_DEVICE_CODE ON)
#    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#
#    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --gpu-architecture=sm_52 -std=c++11 --expt-extended-lambda -ccbin /opt/sw/gcc-5.4.0/bin")
#
#    set(
#        CUDA_LIB_SOURCES
#        "../cuda/forward_project.cu"
#        "../cuda/sart.cu"
#        "../cuda/sirt.cu"
#    )
#
#    set(
#        CUDA_EXAMPLE_SOURCES
#        "../examples/cuda_reconstruction.cpp"
#    )
#
#    set(CUDA_LIB_NAME "tomos_cuda")
#
#    cuda_add_library(${CUDA_LIB_NAME} SHARED ${CUDA_LIB_SOURCES})
#    cuda_add_executable("cuda_reconstruction" ${CUDA_EXAMPLE_SOURCES})
#    target_link_libraries("cuda_reconstruction" ${LIB_NAMES} ${CUDA_LIB_NAME})
#
#endif(CUDA_FOUND)

# add_subdirectory(../test/build ../test/build)

find_package(MPI REQUIRED)
find_package(OpenCV REQUIRED)

if(WITH_MPI AND WITH_ZMQ)
if(MPI_CXX_FOUND)
    set(
        MPI_DIR
        "distributed/"
        )

    set(
        MPI_SOURCES
        "simple_distributed.cpp"
        "run_measurement.cpp"
        "count_messages.cpp"
        )

message("${MPI_CXX_COMPILE_FLAGS}")
message("${MPI_CXX_LINK_FLAGS}")
message("${MPI_CXX_INCLUDE_PATH}")
message("${MPI_CXX_LIBRARIES}")
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    include_directories(${MPI_CXX_INCLUDE_PATH})
	include_directories( ${OpenCV_INCLUDE_DIRS} )


    # make distributed example
    foreach(source_file ${MPI_SOURCES})
        string(REPLACE ".cpp" "" source_name ${source_file})
        add_executable(mpi_${source_name} ${MPI_DIR}${source_file})
        target_link_libraries(mpi_${source_name} ${LIB_NAMES} ${MPI_CXX_LIBRARIES} ${OpenCV_LIBS})
    endforeach(source_file)
else(MPI_CXX_FOUND)
    message("No MPI installation found.")
endif(MPI_CXX_FOUND)
endif(WITH_MPI AND WITH_ZMQ)
