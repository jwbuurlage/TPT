cmake_minimum_required(VERSION 3.0.0)

project(galactica_test)

include_directories(
    "../../include"
    "../../ext/catch/include"
)

# link_directories(
#     "../ext/submodule/lib"
# )

set(TEST_EXEC_NAME "galactica_test")

set(
    TEST_SOURCES
    "../test.cpp"
    "../initialization.cpp"
    "../math.cpp"
    "../geometry.cpp"
    "../../ext/fmt/fmt/format.cc"
)

set(
    LIB_NAMES
    "pthread"
)

## multiple binaries:
# foreach(source_file ${EXAMPLE_SOURCES})
#     string(REPLACE ".cpp" "" source_name ${source_file})
#     add_executable(${source_name} ${EXAMPLE_DIR}${source_file})
#     target_link_libraries( ${source_name} ${LIB_NAMES} )
# endforeach(source_file)

add_executable(${TEST_EXEC_NAME} ${TEST_SOURCES})
target_link_libraries( ${TEST_EXEC_NAME} ${LIB_NAMES} )
