# Copyright 2019-2023 Arnaud Duval
# Copyright 2020 Thibaut Hirschler

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

# CMakeLists.txt file for /src directory

# Grab Lapack and openMP
find_package(
    LAPACK
    REQUIRED
    # BLA_STATIC ON
)

find_package(OpenMP REQUIRED)

# define a macro to build f2py module
macro(CreateF2PYModule)
    set(options)
    set(oneValueArgs NAME DESTINATION)
    set(multiValueArgs FILES LINK_LIBRARIES)

    CMAKE_PARSE_ARGUMENTS(MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.${Python3_SOABI}.so)

    add_custom_target(${MODULE_NAME} ALL
        DEPENDS ${generated_module_file}
    )

    # manage dependencies and crete link command
    set(link_string )
    foreach(dep ${MODULE_LINK_LIBRARIES})
        set(link_string ${link_string} -l${dep})
    endforeach()

    add_custom_command(
        OUTPUT ${generated_module_file}
        COMMAND "${Python3_EXECUTABLE}" -m numpy.f2py
        -m ${MODULE_NAME} #--quiet
        -c --opt=-O2
        -I${CMAKE_Fortran_MODULE_DIRECTORY}
        -L${LIBRARY_OUTPUT_PATH}
        ${link_string}
        -llapack -lblas
        ${FCOMPILER}
        --f90flags='-fopenmp' -lgomp
        ${MODULE_FILES}
        DEPENDS ${MODULE_FILES}
        DEPENDS ${MODULE_LINK_LIBRARIES}
        COMMENT "[F2PY] Building Fortran to Python interface module ${MODULE_NAME}"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND  ${CMAKE_COMMAND} -E copy ${generated_module_file} ${MODULE_DESTINATION}/${MODULE_NAME}.${Python3_SOABI}.so
    )

endmacro()


# Fortran object library for use with scikit-build
add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python3::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

# Macro for f2py build using scikit-build
macro(CreateF2PYModule2)
    set(options)
    set(oneValueArgs NAME DESTINATION)
    # set(multiValueArgs FILES LINK_LIBRARIES)          # Not used for themoment, module will only link with base library
    string(JOIN " " FILES_WITH_SPACES ${MODULE_FILES})

    # Debug print
    cmake_parse_arguments(MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    cmake_print_variables(MODULE_NAME)
    cmake_print_variables(MODULE_FILES)
    cmake_print_variables(FILES_WITH_SPACES)

    add_custom_command(
        OUTPUT "${MODULE_NAME}module.c" "${MODULE_NAME}-f2pywrappers.f"
        DEPENDS ${MODULE_FILES}
        VERBATIM
        COMMAND "${Python3_EXECUTABLE}" -m numpy.f2py
                "${FILES_WITH_SPACES}" -m "${MODULE_NAME}" --lower)

    python_add_library(
        ${MODULE_NAME} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}module.c"  #Generated
        "${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}-f2pywrappers.f"
        "${MODULE_FILES}" WITH_SOABI        # Fortran sources
    )


    target_link_libraries(${MODULE_NAME}
                          PRIVATE fortranobject
                          PRIVATE base)
    target_include_directories(${MODULE_NAME} PRIVATE "${Python_NumPy_INCLUDE_DIRS}")

    install(TARGETS ${MODULE_NAME} DESTINATION ${MODULE_DESTINATION})
endmacro(CreateF2PYModule2)



# Base Fortran modules
add_library(base        parameters.f90
                        igaparmetrization.f90
                        embeddedmapping.f90
                        main/material_lib.f90
                        main/material_lib_HO.f90
                        main/dersbasisfuns.f90
                        main/shap.f90
                        main/shap_HO.f90
                        main/shapPress.f
                        main/operateurs.f
                        main/GaussLegendre.f90
                        main/UELMAT.f90
                        main/UELMAT_HO.f90
                        main/stiffmatrix.F90
                        main/stiffmatrix_HO.f90
                        main/UMASSMAT.f
                        main/evaluatenurbsfcts.f
                        main/evaluateStress.f
                        main/point_inversion.f90
                        main/projection.f90
                        main/isElemOnFace.f90
                        main/shell/nurbsbasisfuns.f
                        main/shell/USFMEM_shell.f
                        main/shell/USFBND_shell.f
                        main/shell/curvilinearCoordinates.f
                        main/shell/shapPress_shell.f
                        main/shell/UELMAT_shell.f
                        main/shell/ComputeMemStrain_shell.f
                        main/shell/ComputeBndStrain_shell.f
                        main/shell/evaluateStress_shell.f
                        main/embedded/mapEntities2parametricSpaces.f
                        main/embedded/UELMAT_embeddedshell.f
                        main/embedded/UELMAT_embeddedvol.f90
                        )

# Create Python modules massmtrx, DOF, reconstructionSOL, stiffmtrx_elemstorage

CreateF2PYModule(NAME massmtrx
                 FILES  ${CMAKE_CURRENT_SOURCE_DIR}/main/compute_CMASSMatrix.f
                 LINK_LIBRARIES base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

CreateF2PYModule(NAME DOF
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/main/indDOF.f
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

CreateF2PYModule(NAME reconstructionSOL
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/main/reconstruction.f
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

CreateF2PYModule(NAME stiffmtrx_elemstorage
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/main/build_system_elemStorage.f90
                 LINK_LIBRARIES base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})
CreateF2PYModule(NAME utils
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/utils/GaussPts.f90
                 LINK_LIBRARIES base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

# Create a Python module for basic functions, already contained in compiled library "base"
CreateF2PYModule(NAME tools
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/main/dersbasisfuns.f90
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})



# Define subdirectories
subdirs(postprocessing optim fitting coupling iga_wq_mf)# solver) #preprocessing

# solver and preprocessing directories contain pure python script that should be copied into build directory
add_custom_target(copy_solver ALL)
add_custom_target(copy_preprocessing ALL)

file(GLOB_RECURSE
    solver
    LIST_DIRECTORIES true
    ${CMAKE_CURRENT_SOURCE_DIR}/solver)

file(GLOB_RECURSE
    preprocessing
    LIST_DIRECTORIES true
    ${CMAKE_CURRENT_SOURCE_DIR}/preprocessing)

add_custom_command(TARGET copy_solver PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/solver ${PYTHON_MODULES_OUTPUT_PATH}/solver
                   DEPENDS solver)

add_custom_command(TARGET copy_preprocessing PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/preprocessing ${PYTHON_MODULES_OUTPUT_PATH}/preprocessing
                   DEPENDS preprocessing)

