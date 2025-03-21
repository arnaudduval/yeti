# Copyright 2020-2023 Arnaud Duval

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

# Add validation tests

set(BENCH_PATH ${CMAKE_BINARY_DIR}/benchs)

set(RUN_LONG_TESTS OFF CACHE BOOL "Run long duration tests")

# Compute eigenvalues, 3D solid
add_test(NAME test_platevib COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/test_platevib.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)

# Optimize with eigenvalues objective function
# Finitie differences
if(RUN_LONG_TESTS)
    add_test(NAME OPTplatevib_FD COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/OPTplatevib_FD.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)
endif(RUN_LONG_TESTS)
# Semi analytic
if(RUN_LONG_TESTS)
    add_test(NAME OPTplatevib_semi_AN COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/OPTplatevib_semi_AN.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)
endif(RUN_LONG_TESTS)
# Full analytic
add_test(NAME test_grad_an_platevib COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/test_grad_an_platevib.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)

# Optimize compliance with volume kept constant, full analytic gradients, 2D solid structure
add_test(NAME test_grad_plate_with_hole COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/plateWithHole/test_grad_plate_with_hole.py WORKING_DIRECTORY ${BENCH_PATH}/plateWithHole)

# Optimize compliance with volume kept constant, full analytic gradients, Kirchhoff-Love shell
add_test(NAME test_grad_square_shell_roof COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/squareShellRoof/test_grad_square_shell_roof.py WORKING_DIRECTORY ${BENCH_PATH}/squareShellRoof)

# Optimize compliance with volume kept constant, full analytic gradients, 3D solid
if(RUN_LONG_TESTS)
    add_test(NAME OPT3-volumetricBeam COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT3-volumetricBeam/OPT3-volumetricBeam.py WORKING_DIRECTORY ${BENCH_PATH}/OPT3-volumetricBeam)
endif(RUN_LONG_TESTS)

# Optimize max displacement with volume kept constant, full analytic gradients, 2D solid structure
add_test(NAME test_grad_cantilever_beam COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/cantileverBeam/test_grad_cantilever_beam.py WORKING_DIRECTORY ${BENCH_PATH}/cantileverBeam)

# Optimize max displacement with volume kept constant, full analytic gradients, Kirchhoff-love shell
add_test(NAME test_grad_square_shell_roof_disp COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/squareShellRoofDisp/test_grad_square_shell_roof_disp.py WORKING_DIRECTORY ${BENCH_PATH}/squareShellRoofDisp)

# Optimize max Von Mises stress, full analytic gradients, 2D solid structure
add_test(NAME test_grad_tensile_specimen COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/tensileSpecimen/test_grad_tensile_specimen.py WORKING_DIRECTORY ${BENCH_PATH}/tensileSpecimen)

# Optimize bending moment with volume kept constant, full analytic gradients, Kirchhoff-Love shell
add_test(NAME test_grad_catenary COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/catenary/test_grad_catenary.py WORKING_DIRECTORY ${BENCH_PATH}/catenary)

# optimize eigenfrequencies with volume kept constant, full analytic gradients, 3D solid
add_test(NAME test_grad_elephant_trunk COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/ElephantTrunk/test_grad_elephant_trunk.py WORKING_DIRECTORY ${BENCH_PATH}/ElephantTrunk)

# Static analysis with distributed pressure
add_test(NAME test_dload_beam_triangular_1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/test_dload_beam_triangular_1.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME test_dload_beam_triangular_2 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/test_dload_beam_triangular_2.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME test_dload_beam_uniform COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/test_dload_beam_dist.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME test_dload_beam_uniform_field COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/test_dload_beam_dist_field.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME test_dload_beam_unused_distrib COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/test_dload_beam_unused_distrib.py WORKING_DIRECTORY ${BENCH_PATH}/dload)

# Combination of different loads with extra parameters, given in different orders
add_test(NAME test_load_combination COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/test_beam_dload_combination.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)
add_test(NAME test_load_combination_order COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/test_beam_dload_triangle_centrif.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)
add_test(NAME test_load_combination_order_embedded COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/test_embeddedbeam_load_combination.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)

# Static shell analysis where a single element is loaded with pressure
add_test(NAME test_shell_1_element_load_1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/shell_1_element_load/test_shell_load_1elt.py WORKING_DIRECTORY ${BENCH_PATH}/shell_1_element_load)
add_test(NAME test_shell_1_element_load_2 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/shell_1_element_load/test_shell_load_9elts.py WORKING_DIRECTORY ${BENCH_PATH}/shell_1_element_load)

# Static analyis of a 2D palne stress plate with a round hole and least square projection postprocessing
add_test(NAME test_leastsquare_projection_1patch COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/1_patch/test_leastsquare_projection_1_patch.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/1_patch)
add_test(NAME test_leastsquare_projection_2patches COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/2_patches/test_leastsquare_projection_2_patches.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/2_patches)
add_test(NAME test_leastsquare_projection_2patchescplg COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/2_patches_coupling/test_leastsquare_projection_2_patchescplg.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/2_patches_coupling)

# Coupling static analysis and post processing of interface patches
add_test(NAME test_coupling_postprocessing COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/postproc_interfaces/test_coupling_postprocessing.py WORKING_DIRECTORY ${BENCH_PATH}/postproc_interfaces)

# Embedded solid static analysis - geometrically nonlinear mapping
add_test(NAME test_embedded_solid_nonlinear_mapping COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/EMBDSOL1-patch_test_NL_map/test_embedded_solid_nonlinear_mapping.py WORKING_DIRECTORY ${BENCH_PATH}/EMBDSOL1-patch_test_NL_map)

# Embedded solid static analysis - geometrically nonlinear embedded solid
add_test(NAME test_nonlinear_embedded_solid COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/EMBDSOL2-patch_test_NL_embded/test_nonlinear_embedded_solid.py WORKING_DIRECTORY ${BENCH_PATH}/EMBDSOL2-patch_test_NL_embded)

# Couling between solid patches with projected Gauss points
add_test(NAME test_coupling_u5_2d_compatible
COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/coupling_U5/test_coupling_u5_2d_compatible.py
         WORKING_DIRECTORY ${BENCH_PATH}/coupling_U5)
add_test(NAME test_coupling_u5
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/coupling_U5/test_coupling_u5.py
         WORKING_DIRECTORY ${BENCH_PATH}/coupling_U5)
# Comparison between U5 and U4 method
add_test(NAME test_compare_u4_u5_2d
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/coupling_U5/test_compare_u4_u5_2d.py
         WORKING_DIRECTORY ${BENCH_PATH}/coupling_U5)

# plate with hole with Abaqus parameters
add_test(NAME test_ROM1_plateWithHole COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/ROM1-plateWithHole/test_romplatewithhole.py WORKING_DIRECTORY ${BENCH_PATH}/ROM1-plateWithHole)

# Postprocessing of an IGA model to plot faces numbers in a VTU file
add_test(NAME test_post_faces COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/post_face/test_post_faces.py WORKING_DIRECTORY ${BENCH_PATH}/post_face)

# Benchs for centrifugal body force
add_test(NAME test_centrif_u1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1/test_centrif_u1.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1)
add_test(NAME test_centrif_u10 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U10/test_centrif_u10.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U10)
add_test(NAME test_centrif_u1_c0 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_C0/test_centrif_u1_c0.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_C0)
add_test(NAME test_centrif_u1_u1_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_U1_1load/test_centrif_u1_u1_1load.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_U1_1load)
add_test(NAME test_centrif_u1_u1_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_U1_2loads/test_centrif_u1_u1_2loads.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_U1_2loads)
add_test(NAME test_centrif_u1_cpl_u1_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_1load/test_centrif_u1_cpl_u1_1load.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_1load)
add_test(NAME test_centrif_u1_cpl_u1_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_2loads/test_centrif_u1_cpl_u1_2loads.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_2loads)
add_test(NAME test_centrif_u1_cpl_u10_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_1load/test_centrif_u1_cpl_u10_1load.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_1load)
add_test(NAME test_centrif_u1_cpl_u10_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_2loads/test_centrif_u1_cpl_u10_2loads.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_2loads)

# Benchs for compliance optimisation with centrifugal force
add_test(NAME test_grad_centrif_coupling_u1_c0 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_C0/test_grad_centrif_coupling_u1_c0.py WORKING_DIRECTORY ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_C0)
add_test(NAME test_grad_centrif_u1_cplg_u1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_cplg_U1/test_grad_centrif_u1_cplg_u1.py WORKING_DIRECTORY ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_cplg_U1)

# Solid element sensibilities
# 1 - centrifugal body force loading
add_test(NAME test_grad_solid_unit_cube_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_solid/test_grad_solid_unit_cube_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_solid
)

# Embedded solid element sensibilities
# 1 - volume gradient
add_test(NAME test_gradvol_embdsol_cube_nonli_hull
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradvol_embdsol_cube_nonli_hull.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
add_test(NAME test_gradvol_embdsol_cube_nonli_embd
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradvol_embdsol_cube_nonli_embd.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
add_test(NAME test_gradvol_embdsol_cube_nonli_hull_embd
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradvol_embdsol_cube_nonli_hull_embd.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)

# 2 - compliance gradient
# linear, pressure loading
add_test(NAME test_gradcomp_embdsol_cube_lin_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_lin_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear, centrifugal body force loading
add_test(NAME test_gradcomp_embdsol_cube_lin_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_lin_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear mapping ot embedded entity, non linear hull, pressure loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_hull_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_hull_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear mapping of embedded entity, non linear hull, centrifugal body force loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_hull_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_hull_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, linear mapping of hull, pressure loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_embd_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_embd_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, linear mapping of hull, centrifugal body force loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_embd_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_embd_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, non linear mapping of hull, pressure loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_hull_embd_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_hull_embd_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, non linear mapping of hull, centrifugal body force loading
add_test(NAME test_gradcomp_embdsol_cube_nonli_hull_embd_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/test_gradcomp_embdsol_cube_nonli_hull_embd_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)

# VTU Bezier output (output on 2 and 8 patchs may be redundant)
add_test(NAME test_vtu_output_2patchs
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/test_vtu_output_2patchs.py
         WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
)
add_test(NAME test_vtu_output_8patchs
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/test_vtu_output_8patchs.py
         WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
)
# add_test(NAME test_vtu_output_2patchs_coupling
#          COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/test_vtu_output_2patchs_coupling.py
#          WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
# )


# Benchs for sensibilities computation on embedded solid element
#add_test(NAME grad_embdsol_unit_square COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/grad_embdsol_unit_square.py WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol)

set_tests_properties(
                     test_platevib
                     test_grad_an_platevib
                     test_grad_plate_with_hole
                     test_grad_square_shell_roof
                     test_grad_cantilever_beam
                     test_grad_square_shell_roof_disp
                     test_grad_tensile_specimen
                     test_grad_catenary
                     test_grad_elephant_trunk
                     test_dload_beam_triangular_1
                     test_dload_beam_triangular_2
                     test_dload_beam_uniform
                     test_dload_beam_uniform_field
                     test_dload_beam_unused_distrib
                     test_load_combination
                     test_load_combination_order
                     test_load_combination_order_embedded
                     test_shell_1_element_load_1
                     test_shell_1_element_load_2
                     test_leastsquare_projection_1patch
                     test_leastsquare_projection_2patches
                     test_leastsquare_projection_2patchescplg
                     test_coupling_postprocessing
                     test_embedded_solid_nonlinear_mapping
                     test_nonlinear_embedded_solid
                     test_coupling_u5_2d_compatible
                     test_coupling_u5
                     test_compare_u4_u5_2d
                     test_ROM1_plateWithHole
                     test_post_faces
                     test_centrif_u1
                     test_centrif_u10
                     test_centrif_u1_c0
                     test_centrif_u1_u1_1load
                     test_centrif_u1_u1_2loads
                     test_centrif_u1_cpl_u1_1load
                     test_centrif_u1_cpl_u1_2loads
                     test_centrif_u1_cpl_u10_1load
                     test_centrif_u1_cpl_u10_2loads
                     test_grad_centrif_coupling_u1_c0
                     test_grad_centrif_u1_cplg_u1
                     test_grad_solid_unit_cube_centrif
                     test_gradvol_embdsol_cube_nonli_hull
                     test_gradvol_embdsol_cube_nonli_embd
                     test_gradvol_embdsol_cube_nonli_hull_embd
                     test_gradcomp_embdsol_cube_lin_pressure
                     test_gradcomp_embdsol_cube_lin_centrif
                     test_gradcomp_embdsol_cube_nonli_hull_pressure
                     test_gradcomp_embdsol_cube_nonli_hull_centrif
                     test_gradcomp_embdsol_cube_nonli_embd_pressure
                     test_gradcomp_embdsol_cube_nonli_embd_centrif
                     test_gradcomp_embdsol_cube_nonli_hull_embd_pressure
                     test_gradcomp_embdsol_cube_nonli_hull_embd_centrif
                     test_vtu_output_2patchs
                     test_vtu_output_8patchs
                    #  test_vtu_output_2patchs_coupling
                     PROPERTIES ENVIRONMENT PYTHONPATH=${PYTHON_LIB_OUTPUT_PATH}:$ENV{PYTHONPATH})

# Properties for long tests
if(RUN_LONG_TESTS)
    set_tests_properties(OPTplatevib_semi_AN
                         OPTplatevib_FD
                         OPT3-volumetricBeam
                         PROPERTIES ENVIRONMENT PYTHONPATH=${PYTHON_LIB_OUTPUT_PATH}:$ENV{PYTHONPATH}
                         )
endif(RUN_LONG_TESTS)



