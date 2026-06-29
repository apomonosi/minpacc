cmake_minimum_required(VERSION 3.14)
project(minpacc VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS    "Build test suite"       ON)

add_subdirectory(src)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
