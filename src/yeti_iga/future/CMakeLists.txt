cmake_minimum_required(VERSION 3.15)
project(BSplineLib VERSION 1.0 LANGUAGES CXX)

# Options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Default build mode : Release
# Mode Release par défaut
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")

# Find pybind11 (from system or pip package) and OpenMP
find_package(pybind11 REQUIRED)
find_package(OpenMP REQUIRED)


pybind11_add_module(bspline
    src/bspline/BSpline.cpp
    src/bspline/BSplineTensor.cpp
    src/bspline/Patch.cpp
    src/bspline/IGAAssembler.cpp
    src/python/bindings.cpp)

target_include_directories(bspline PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(bspline PUBLIC OpenMP::OpenMP_CXX)

# Flags Debug : pas d'optimisation, debug symbols
target_compile_options(bspline PRIVATE
    $<$<CONFIG:Debug>:
        -O0
        -g
    >
)

# Flags Release : optimisation complète
target_compile_options(bspline PRIVATE
    $<$<CONFIG:Release>:
        -O3
        -march=native
        -mtune=native
    >
)