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

# CMakeLists.txt file for /src/postprocessing directory

CreateF2PYModule(NAME postproc
                 FILES ${CMAKE_CURRENT_SOURCE_DIR}/postproc_basisfuncts.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/postproc_vtu.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/generate_polydata_vtk.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/writepatch.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/evaluate_disp.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/evaluate_stress.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/postproc_cplgInterface.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_solid.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_solid_HO.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_shell.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_embdedshell.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_embdedsolid.f
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_gram_matrix.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/compute_svars_solid_rhs.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/generate_proj_vtu.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/generate_coupling_vtu.f90
                       #${CMAKE_CURRENT_SOURCE_DIR}/postproc_cplg.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/postproc_faces.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/postproc_curve_2D.f90
                       ${CMAKE_CURRENT_SOURCE_DIR}/postproc_vtu_bezier.f90
                 LINK_LIBRARIES base coupling
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/postprocessing)

add_custom_target(copy_postprocessing ALL)

add_custom_command(TARGET copy_postprocessing PRE_BUILD
                   COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PYTHON_MODULES_OUTPUT_PATH}/postprocessing/__init__.py
                   DEPENDS postproc)
