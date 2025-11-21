cmake_minimum_required(VERSION 3.15)
project(BSplineLib VERSION 1.0 LANGUAGES CXX)

# Options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find pybind11 (from system or pip package)
find_package(pybind11 REQUIRED)

pybind11_add_module(bspline
    src/bspline/BSpline.cpp
    src/bspline/BSplineTensor.cpp
    src/python/bindings.cpp)

target_include_directories(bspline PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)