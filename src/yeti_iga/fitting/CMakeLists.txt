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

# CMakeLists.txt file for /src/fitting directory

CreateF2PYModule(NAME interpolate
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/interpolationMatrix.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/grevillePhysicalCoords.f
                 LINK_LIBRARIES base
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/fitting)

# CreateF2PYModule(NAME evaluate
#                  FILES ${CMAKE_CURRENT_SOURCE_DIR}/grevillePhysicalCoords.f
#                        ${CMAKE_CURRENT_SOURCE_DIR}/grevilleAbscissae.f
#                  LINK_LIBRARIES base
#                  DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/fitting)

add_custom_target(copy_fitting ALL)

add_custom_command(TARGET copy_fitting PRE_BUILD
                   COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PYTHON_MODULES_OUTPUT_PATH}/fitting/__init__.py
                   DEPENDS interpolate
                   DEPENDS base
                #    DEPENDS evaluate
                )
