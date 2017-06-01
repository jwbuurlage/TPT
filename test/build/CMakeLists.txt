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

# reignore warnings from cpptoml, there is probably a nicer way. TODO
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ../../../ext/cpptoml/include")

add_executable(${TEST_EXEC_NAME} ${TEST_SOURCES})
target_link_libraries( ${TEST_EXEC_NAME} ${LIB_NAMES} )
