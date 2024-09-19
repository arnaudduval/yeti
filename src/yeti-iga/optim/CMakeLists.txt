# Copyright 2019-2020 Arnaud Duval
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

# CMakeLists.txt file for /src/optim directory

# Optim library
set(optim_src_file ${CMAKE_CURRENT_SOURCE_DIR}/gradUELMAT_solid.f90
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradUELMAT_shell.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradUELMAT_embeddedshell.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradCPLGMAT_disp.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradCPLGMAT_rot.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradUELMASSMAT_solid.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/gradUELMAT_embeddedsolid.f90
                     )


add_library(optim ${optim_src_file})
target_link_libraries(optim base)

CreateF2PYModule(NAME volume
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/volume.f90
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradvolume.f90
                 LINK_LIBRARIES base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

# Warning : following code create an import error in Python if dependencies are in order 'base optim coupling'
CreateF2PYModule(NAME vibration
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradvibration_fullAN.f
                 LINK_LIBRARIES optim base coupling
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

#CreateF2PYModule(NAME gradvolume
#                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradvolume.f
#                 LINK_LIBRARIES base
#                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

CreateF2PYModule(NAME gradcompliance
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradcompliance_semiAN.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/gradcompliance_fullAN.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/gradcoupling.f
                 LINK_LIBRARIES optim base coupling
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

CreateF2PYModule(NAME gradcoupling
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradcoupling.f
                 LINK_LIBRARIES coupling optim base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

#CreateF2PYModule(NAME graddisplacement
#                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/graddisplacement_fullAN.f
#                       ${CMAKE_CURRENT_SOURCE_DIR}/adjointRHS_disp.f
#                 LINK_LIBRARIES  optim base
#                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

CreateF2PYModule(NAME gradadjointwork
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/gradtotalwork_fullAN.f90
                 LINK_LIBRARIES optim base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

CreateF2PYModule(NAME stress
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/aggregate_stress.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/adjointRHS_aggregStress.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/partialDERV_aggregStress.f
                 LINK_LIBRARIES optim base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

CreateF2PYModule(NAME disp
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/aggregate_disp.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/adjointRHS_disp.f
                 LINK_LIBRARIES optim base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/optim)

add_custom_target(copy_optim ALL)

add_custom_command(TARGET copy_optim PRE_BUILD
                   COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PYTHON_MODULES_OUTPUT_PATH}/optim/__init__.py
                   DEPENDS optim)


