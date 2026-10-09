# Copyright 2019-2020 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the terms
# of the GNU Lesser General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Yeti. If not, see <https://www.gnu.org/licenses/>



#######
# Main CMakeLists.txt file

cmake_minimum_required(VERSION 3.21.0)

# Test if using scikit-build or not
if("${SKBUILD_PROJECT_NAME}" STREQUAL "")
    message(STATUS "Building without scikit-build")
    set(LEGACY TRUE)
else()
    message(STATUS "Building with scikit-build")
    set(LEGACY FALSE)
endif()


# Project definition, use languages Fortran and C
project(${SKBUILD_PROJECT_NAME} LANGUAGES Fortran C)

# gfortran < 6 doesn't support deferred-length allocatable CHARACTER components
# in derived types (src/yeti_iga/main/m_nurbspatch.f90's elt_type_patch/
# tensor_patch), failing with a cryptic "Error: Deferred-length character
# component ... is not yet supported" deep inside the Fortran compile --
# confirmed to fail with gfortran 4.8.5 and succeed with 9.3.1. On HPC/RHEL-style
# systems the system-default `f95`/`gfortran` is often this old, even when a
# newer compiler (e.g. via devtoolset) is already selected for C/CXX: CMake's
# Fortran compiler search is independent of CC/CXX, so it must be pointed at
# the newer compiler explicitly if they differ.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0)
    message(FATAL_ERROR
        "gfortran ${CMAKE_Fortran_COMPILER_VERSION} (${CMAKE_Fortran_COMPILER}) is too old "
        "to build Yeti: deferred-length allocatable CHARACTER components in derived "
        "types (src/yeti_iga/main/m_nurbspatch.f90) require gfortran >= 6. This often "
        "happens on HPC/RHEL-style systems where the system-default Fortran compiler is "
        "old even though a newer C/C++ toolchain (e.g. devtoolset) is already in use. "
        "Point CMake at a newer gfortran explicitly, e.g.:\n"
        "  CMAKE_ARGS=\"-DCMAKE_Fortran_COMPILER=/path/to/newer/gfortran\" pipx run build\n"
        "or export FC=/path/to/newer/gfortran before building.")
endif()

# CMake Variable printer for debug
include(CMakePrintHelpers)

cmake_print_variables(LEGACY)

# Set version number

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 1)

set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})


# Path to additional CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Project needs Python
# set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3
             COMPONENTS Interpreter NumPy # Development
             REQUIRED)

# F2PY headers
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)


find_package(Sphinx)

# Fortran compilation flags
SET(CMAKE_Fortran_FLAGS "-O2 -fPIC")

# Define location of Fortran module files
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
# Define location of libraries
SET(LIBRARY_OUTPUT_PATH            ${CMAKE_BINARY_DIR}/lib)
# Define output path for Python modules
if(LEGACY)
    SET(PYTHON_MODULES_OUTPUT_PATH     ${CMAKE_BINARY_DIR}/lib/python/yeti_iga)
    SET(PYTHON_LIB_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/python)
    file(MAKE_DIRECTORY ${PYTHON_MODULES_OUTPUT_PATH})
else()
    SET(PYTHON_MODULES_OUTPUT_PATH     ${SKBUILD_PROJECT_NAME})
endif()


# Specific f2py argument for Intel compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(FCOMPILER "--fcompiler=intelem")
endif()

subdirs(src)
add_subdirectory(src/yeti_iga/future)

if(LEGACY)
    enable_testing()
    # Copy benchs to build directory
    add_custom_target(benchs_init DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/benchs)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/benchs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    subdirs(benchs)

    # Copy pyproject.toml file to build/doc directory to handle package version
    add_custom_target(doc_init DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/docs/pyproject.toml)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/docs)
    subdirs(docs)
endif()