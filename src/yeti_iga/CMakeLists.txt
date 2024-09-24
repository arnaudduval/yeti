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

if(LEGACY)
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
endif()
