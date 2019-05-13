cmake_minimum_required(VERSION 3.0.0)
project(tpt)

include(FeatureSummary)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(FOR_RELEASE "Build with all optimizations turned on")
option(DISTRIBUTED "Build distributed examples based on MPI and ZMQ")
option(PYTHON_BINDINGS "Build python bindings")
option(TPT_LIB_ONLY "Do not build examples and tools")

feature_summary(WHAT ALL)

set(LIB_NAME "tpt")
add_library(${LIB_NAME} INTERFACE)

# Add compiler options for TPT applications
target_compile_options(${LIB_NAME} INTERFACE
    -std=c++17
    -Wall
    -Wextra
    -Wno-missing-braces
    -Wfatal-errors
    -Werror
    -pedantic
    -O3
    -march=native)

# Add include directories
target_include_directories(${LIB_NAME} INTERFACE "include")

# Ignore warnings from external libs
target_compile_definitions(${LIB_NAME} INTERFACE
     -DFMT_HEADER_ONLY)
target_include_directories(${LIB_NAME} SYSTEM INTERFACE
    "ext/fmt/include/"
    "ext/cpptoml/include/"
    "ext/bulk/include"
)

# Optional optimizations
if(FOR_RELEASE)
    # Add NDEBUG define
    target_compile_options(${LIB_NAME} INTERFACE -DNDEBUG)
else(FOR_RELEASE)
    # Add debug flag
    target_compile_options(${LIB_NAME} INTERFACE -g)
endif(FOR_RELEASE)

# Set output for binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# Add GLM as dependency
find_package(glm REQUIRED)

# Add Boost as dependency
find_package (Boost 1.61 REQUIRED)

# Add CLI11 as dependency
add_subdirectory("ext/CLI11")

# Link against external libraries
target_link_libraries(
    ${LIB_NAME} INTERFACE
    pthread
    stdc++fs
    glm
    Boost::boost
    CLI11::CLI11
)

# -----------------------------------------------------
# Basic (non-MPI) Examples

if (NOT TPT_LIB_ONLY)
  set(EXAMPLE_DIR "examples/")

  set(
      EXAMPLE_SOURCES
      "generate_mtx.cpp"
      "minimal_example.cpp"
  )

  foreach(source_file ${EXAMPLE_SOURCES})
      string(REPLACE ".cpp" "" source_name ${source_file})
      add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
      target_link_libraries(${source_name} ${LIB_NAME})
  endforeach(source_file)
endif()

# -----------------------------------------------------
# Basic tools

# -----------------------------------------------------
# MPI specific code

if(DISTRIBUTED AND (NOT TPT_LIB_ONLY))
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    # Find ZeroMQ
    find_package(ZeroMQ QUIET)

    if (ZeroMQ_FOUND)
        add_library(zmq INTERFACE)
        target_include_directories(zmq INTERFACE ${ZeroMQ_INCLUDE_DIR})
        target_link_libraries(zmq INTERFACE ${ZeroMQ_LIBRARY})
    else()
        message("'zmq' not installed on the system, building from source...")

        execute_process(COMMAND git submodule update --init -- ext/libzmq
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

        set(ZMQ_BUILD_TESTS OFF CACHE BOOL "disable tests" FORCE)
        set(WITH_PERF_TOOL OFF CACHE BOOL "disable perf-tools" FORCE)
        add_subdirectory(${CMAKE_SOURCE_DIR}/ext/libzmq)
        set(ZMQ_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/libzmq/include)

        # ZeroMQ names their target libzmq, which is inconsistent => create a ghost dependency
        add_library(zmq INTERFACE)
        target_link_libraries(zmq INTERFACE libzmq)
    endif()

    target_include_directories(${LIB_NAME} INTERFACE
        "ext/tomopackets/include"
        "ext/cppzmq"
        "ext/libzmq/include"
    )

    target_link_libraries(${LIB_NAME} INTERFACE "zmq")

    set(
        TOOLS_DIR
        "tools/"
        )

    set(
        TOOLS_SOURCES
        "grcb.cpp"
        "stats.cpp"
        "overlap.cpp"
        "partition.cpp"
        "run_measurement.cpp"
        "count_messages.cpp"
        )

    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    include_directories(${MPI_CXX_INCLUDE_PATH})

    # make distributed example
    foreach(source_file ${TOOLS_SOURCES})
        string(REPLACE ".cpp" "" source_name ${source_file})
        add_executable(tpt_${source_name} ${TOOLS_DIR}${source_file})
        target_link_libraries(tpt_${source_name} ${MPI_CXX_LIBRARIES} ${LIB_NAME})
    endforeach(source_file)
else(MPI_CXX_FOUND)
    message("No MPI installation found.")
endif(MPI_CXX_FOUND)
endif()

if(PYTHON_BINDINGS)
add_subdirectory(python)
endif(PYTHON_BINDINGS)
