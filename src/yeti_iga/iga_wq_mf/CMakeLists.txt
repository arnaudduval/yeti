# Copyright 2019-2021 Arnaud Duval
# Copyright 2022 Joaquin Cornejo 

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

# CMakeLists.txt file for /src/iga_wq_mf directory

# Add library
add_library(iga_wq_mf 
            "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/algebra.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/bspline.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/iga_wq_methods.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/tensor_algebra.f90"
)

target_link_libraries(iga_wq_mf base OpenMP::OpenMP_Fortran ${LAPACK_lapack_LIBRARY})  

CreateF2PYModule(NAME basis_weights
                 FILES  "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/iga_wq_basis_weights.f90"
                 LINK_LIBRARIES iga_wq_mf base 
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/iga_wq_mf
)

CreateF2PYModule(NAME assembly
                 FILES  "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/iga_assembly.f90"
                        "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/wq_assembly.f90"
                 LINK_LIBRARIES iga_wq_mf base 
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/iga_wq_mf
)

CreateF2PYModule(NAME solver
                 FILES  "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/mf_iga_solver.f90"
                        "${CMAKE_CURRENT_SOURCE_DIR}/fsrc/mf_wq_solver.f90"  
                 LINK_LIBRARIES iga_wq_mf base 
                 DESTINATION ${PYTHON_MODULES_OUTPUT_PATH}/iga_wq_mf
)

add_custom_target(copy_iga_wq_mf ALL)
                 
add_custom_command(TARGET copy_iga_wq_mf PRE_BUILD
                   COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${PYTHON_MODULES_OUTPUT_PATH}/iga_wq_mf/__init__.py
                   DEPENDS iga_wq_mf)
