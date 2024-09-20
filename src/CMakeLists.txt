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
        COMMAND ${F2PY_EXECUTABLE}
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
                        main/stiffmatrix.f90
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

add_library(nlg         main/UGEOMAT.f
                        main/UELEXTVEC.f90
                        main/UELINTVEC.f90
                        main/UELTANMAT.f90
                        main/stiffgeommatrix.f90
                        main/geommatrix.f90
                        main/geomvector.f90

                        )
target_link_libraries(nlg base)


# Create Python modules massmtrx, DOF, reconstructionSOL, stiffmtrx_elemstorage, stiffmtrx_elemstorageNLG
CreateF2PYModule(NAME tanmtrx_elemStorageNLG
                 FILES  ${CMAKE_CURRENT_SOURCE_DIR}/main/tanmtrx_elemStorage.f90
                 LINK_LIBRARIES base nlg
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

CreateF2PYModule(NAME intvec_elemStorageNLG
                 FILES  ${CMAKE_CURRENT_SOURCE_DIR}/main/intvec_elemStorage.f90
                 LINK_LIBRARIES base nlg
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

CreateF2PYModule(NAME extvec_elemStorageNLG
                 FILES  ${CMAKE_CURRENT_SOURCE_DIR}/main/extvec_elemStorage.f90
                 LINK_LIBRARIES base nlg
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH})

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

