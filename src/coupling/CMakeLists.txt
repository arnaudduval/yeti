# Copyright 2019-2021 Arnaud Duval

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

# CMakeLists.txt file for /src/coupling directory

# Coupling library
set(coupling_src_file ${CMAKE_CURRENT_SOURCE_DIR}/GPsInfoAtInterface.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/evaluatelgrgefcts.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/UCPGDISPMATRIX.f
                     ${CMAKE_CURRENT_SOURCE_DIR}/UCPGROTMATRIX.f
                     )

                     
add_library(coupling ${coupling_src_file})
target_link_libraries(coupling base)

# 'cplgmatrix' Python module

CreateF2PYModule(NAME cplgmatrix
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/build_cplgmatrix.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/build_cplgmatrixU5.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/build_cplgmatrixU5_wdatabase.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/build_dirichletBCmatrix.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/cplg_U5_wEmbedded.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/weakInterfaceInfos.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/eval_coupling_gpts.f90
                 LINK_LIBRARIES coupling base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/coupling)

add_custom_target(copy_coupling ALL)
                 
add_custom_command(TARGET copy_coupling PRE_BUILD
                   COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PYTHON_MODULES_OUTPUT_PATH}/coupling/__init__.py
                   DEPENDS coupling)
