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
add_test(NAME platevib COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/platevib.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)

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
add_test(NAME OPTplatevib_AN COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/platevib/OPTplatevib_AN.py WORKING_DIRECTORY ${BENCH_PATH}/platevib)

# Optimize compliance with volume kept constant, full analytic gradients, 2D solid structure
add_test(NAME OPT1-plateWithHole COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT1-plateWithHole/OPT1-plateWithHole.py WORKING_DIRECTORY ${BENCH_PATH}/OPT1-plateWithHole)

# Optimize compliance with volume kept constant, full analytic gradients, Kirchhoff-Love shell
add_test(NAME OPT2-squareShellRoof COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT2-squareShellRoof/OPT2-squareShellRoof.py WORKING_DIRECTORY ${BENCH_PATH}/OPT2-squareShellRoof)

# Optimize compliance with volume kept constant, full analytic gradients, 3D solid
if(RUN_LONG_TESTS)
    add_test(NAME OPT3-volumetricBeam COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT3-volumetricBeam/OPT3-volumetricBeam.py WORKING_DIRECTORY ${BENCH_PATH}/OPT3-volumetricBeam)
endif(RUN_LONG_TESTS)

# Optimize max displacement with volume kept constant, full analytic gradients, 2D solid structure
add_test(NAME OPT4-cantileverBeam COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT4-cantileverBeam/OPT4-cantileverBeam.py WORKING_DIRECTORY ${BENCH_PATH}/OPT4-cantileverBeam)

# Optimize max displacement with volume kept constant, full analytic gradients, Kirchhoff-love shell
add_test(NAME OPT5-squareShellRoofDisp COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT5-squareShellRoofDisp/OPT5-squareShellRoofDisp.py WORKING_DIRECTORY ${BENCH_PATH}/OPT5-squareShellRoofDisp)

# Optimize max Von Mises stress, full analytic gradients, 2D solid structure
add_test(NAME OPT6-tensileSpecimen COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT6-tensileSpecimen/OPT6-tensileSpecimen.py WORKING_DIRECTORY ${BENCH_PATH}/OPT6-tensileSpecimen)

# Optimize bending moment with volume kept constant, full analytic gradients, Kirchhoff-Love shell
add_test(NAME OPT7-catenary COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT7-catenary/OPT7-catenary.py WORKING_DIRECTORY ${BENCH_PATH}/OPT7-catenary)

# optimize eigenfrequencies with volume kept constant, full analytic gradients, 3D solid
add_test(NAME OPT8-ElephantTrunk COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/OPT8-ElephantTrunk/OPT8-ElephantTrunk.py WORKING_DIRECTORY ${BENCH_PATH}/OPT8-ElephantTrunk)

# Static analysis with distributed pressure
add_test(NAME dload-triangle-1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/dload_beam_triangular_1.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME dload-triangle-2 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/dload_beam_triangular_2.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME dload-uniform COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/dload_beam_dist.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME dload-uniform-field COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/dload_beam_dist_field.py WORKING_DIRECTORY ${BENCH_PATH}/dload)
add_test(NAME dload_unused_distrib COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/dload/dload_beam_unused_distrib.py WORKING_DIRECTORY ${BENCH_PATH}/dload)

# Combination of different loads with extra parameters, given in different orders
add_test(NAME load_combination COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/beam_dload_combination.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)
add_test(NAME load_combination_order COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/beam_dload_triangle_centrif.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)
add_test(NAME load_combination_order_embedded COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/load_combination/embeddedbeam_load_combination.py WORKING_DIRECTORY ${BENCH_PATH}/load_combination)

# Static shell analysis where a single element is loaded with pressure
add_test(NAME shell_1_element_load-1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/shell_1_element_load/linear_analysis_1elt.py WORKING_DIRECTORY ${BENCH_PATH}/shell_1_element_load)
add_test(NAME shell_1_element_load-2 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/shell_1_element_load/linear_analysis_9elts.py WORKING_DIRECTORY ${BENCH_PATH}/shell_1_element_load)

# Static analyis of a 2D palne stress plate with a round hole and least square projection postprocessing
add_test(NAME leastsquare-projection-1patch COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/1_patch/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/1_patch)
add_test(NAME leastsquare-projection-2patches COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/2_patches/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/2_patches)
add_test(NAME leastsquare-projection-2patchescplg COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/leastsquare_projection/2_patches_coupling/coupling_analysis_mono.py WORKING_DIRECTORY ${BENCH_PATH}/leastsquare_projection/2_patches_coupling)

# Coupling static analysis and post processing of interface patches
add_test(NAME coupling-postprocessing COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/postproc_interfaces/coupling_analysis_pp_interfaces.py WORKING_DIRECTORY ${BENCH_PATH}/postproc_interfaces)

# Embedded solid static analysis - geometrically nonlinear mapping
add_test(NAME embeddedsolid1-NLmap COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/EMBDSOL1-patch_test_NL_map/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/EMBDSOL1-patch_test_NL_map)

# Embedded solid static analysis - geometrically nonlinear embedded solid
add_test(NAME embeddedsolid2-NLembded COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/EMBDSOL2-patch_test_NL_embded/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/EMBDSOL2-patch_test_NL_embded)

# Couling between solid patches with projected Gauss points
add_test(NAME coupling_U5_2D_compatible
COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/coupling_U5/couplingU5_2D_compatible.py
         WORKING_DIRECTORY ${BENCH_PATH}/coupling_U5)
add_test(NAME coupling_U5
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/coupling_U5/couplingU5.py
         WORKING_DIRECTORY ${BENCH_PATH}/coupling_U5)

# plate with hole with Abaqus parameters
add_test(NAME ROM1-plateWithHole COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/ROM1-plateWithHole/romplatewithhole.py WORKING_DIRECTORY ${BENCH_PATH}/ROM1-plateWithHole)

# Postprocessing of an IGA model to plot faces numbers in a VTU file
add_test(NAME post_faces COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/post_face/post_faces.py WORKING_DIRECTORY ${BENCH_PATH}/post_face)

# Benchs for centrifugal body force
add_test(NAME centrif_U1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1)
add_test(NAME centrif_U10 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U10/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U10)
add_test(NAME centrif_U1_C0 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_C0/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_C0)
add_test(NAME centrif_U1_U1_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_U1_1load/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_U1_1load)
add_test(NAME centrif_U1_U1_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_U1_2loads/linear_analysis.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_U1_2loads)
add_test(NAME centrif_U1_cpl_U1_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_1load/couplingU5.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_1load)
add_test(NAME centrif_U1_cpl_U1_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_2loads/couplingU5.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U1_2loads)
add_test(NAME centrif_U1_cpl_U10_1load COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_1load/couplingU5.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_1load)
add_test(NAME centrif_U1_cpl_U10_2loads COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_2loads/couplingU5.py WORKING_DIRECTORY ${BENCH_PATH}/centrif/centrif_U1_cpl_U10_2loads)

# Benchs for compliance optimisation with centrifugal force
add_test(NAME opt_centrif_U1_C0 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_C0/opt_centrif.py WORKING_DIRECTORY ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_C0)
add_test(NAME opt_centrif_U1_cplg_U1 COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_cplg_U1/OPTcentrif.py WORKING_DIRECTORY ${BENCH_PATH}/opt-centrif-coupling/opt_centrif_U1_cplg_U1)

# Solid element sensibilities
# 1 - centrifugal body force loading
add_test(NAME grad_solid_unit_cube_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_solid/grad_solid_unit_cube_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_solid
)

# Embedded solid element sensibilities
# 1 - volume gradient
add_test(NAME gradvol_embdsol_cube_nonli_hull
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradvol_embdsol_cube_nonli_hull.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
add_test(NAME gradvol_embdsol_cube_nonli_embd
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradvol_embdsol_cube_nonli_embd.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
add_test(NAME gradvol_embdsol_cube_nonli_hull_embd
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradvol_embdsol_cube_nonli_hull_embd.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)

# 2 - compliance gradient
# linear, pressure loading
add_test(NAME gradcomp_embdsol_cube_lin_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_lin_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear, centrifugal body force loading
add_test(NAME gradcomp_embdsol_cube_lin_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_lin_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear mapping ot embedded entity, non linear hull, pressure loading
add_test(NAME gradcomp_embdsol_cube_nonli_hull_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_hull_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# linear mapping of embedded entity, non linear hull, centrifugal body force loading
add_test(NAME gradcomp_embdsol_cube_nonli_hull_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_hull_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, linear mapping of hull, pressure loading
add_test(NAME gradcomp_embdsol_cube_nonli_embd_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_embd_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, linear mapping of hull, centrifugal body force loading
add_test(NAME gradcomp_embdsol_cube_nonli_embd_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_embd_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, non linear mapping of hull, pressure loading
add_test(NAME gradcomp_embdsol_cube_nonli_hull_embd_pressure
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_hull_embd_pressure.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)
# nonlinear embedded entity, non linear mapping of hull, centrifugal body force loading
add_test(NAME gradcomp_embdsol_cube_nonli_hull_embd_centrif
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/gradcomp_embdsol_cube_nonli_hull_embd_centrif.py
         WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol
)

# VTU Bezier output (output on 2 and 8 patchs may be redundant)
add_test(NAME vtu_output_2patchs
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/linear_analysis_2patchs.py
         WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
)
add_test(NAME vtu_output_8patchs
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/linear_analysis_8patchs.py
         WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
)
add_test(NAME vtu_output_2patchs_coupling
         COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/vtk_bezier_output/linear_analysis_halftube_2patchs_coupling_U5.py
         WORKING_DIRECTORY ${BENCH_PATH}/vtk_bezier_output
)


# Benchs for sensibilities computation on embedded solid element
#add_test(NAME grad_embdsol_unit_square COMMAND ${Python3_EXECUTABLE} ${BENCH_PATH}/grad_embdsol/grad_embdsol_unit_square.py WORKING_DIRECTORY ${BENCH_PATH}/grad_embdsol)

set_tests_properties(
                     platevib
                     OPTplatevib_AN
                     OPT1-plateWithHole
                     OPT2-squareShellRoof
                     OPT4-cantileverBeam
                     OPT5-squareShellRoofDisp
                     OPT6-tensileSpecimen
                     OPT7-catenary
                     OPT8-ElephantTrunk
                     dload-triangle-1
                     dload-triangle-2
                     dload-uniform
                     dload-uniform-field
                     dload_unused_distrib
                     load_combination
                     load_combination_order
                     load_combination_order_embedded
                     shell_1_element_load-1
                     shell_1_element_load-2
                     leastsquare-projection-1patch
                     leastsquare-projection-2patches
                     leastsquare-projection-2patchescplg
                     coupling-postprocessing
                     embeddedsolid1-NLmap
                     embeddedsolid2-NLembded
                     coupling_U5_2D_compatible
                     coupling_U5
                     ROM1-plateWithHole
                     post_faces
                     centrif_U1
                     centrif_U10
                     centrif_U1_C0
                     centrif_U1_U1_1load
                     centrif_U1_U1_2loads
                     centrif_U1_cpl_U1_1load
                     centrif_U1_cpl_U1_2loads
                     centrif_U1_cpl_U10_1load
                     centrif_U1_cpl_U10_2loads
                     opt_centrif_U1_C0
                     opt_centrif_U1_cplg_U1
                     grad_solid_unit_cube_centrif
                     gradvol_embdsol_cube_nonli_hull
                     gradvol_embdsol_cube_nonli_embd
                     gradvol_embdsol_cube_nonli_hull_embd
                     gradcomp_embdsol_cube_lin_pressure
                     gradcomp_embdsol_cube_lin_centrif
                     gradcomp_embdsol_cube_nonli_hull_pressure
                     gradcomp_embdsol_cube_nonli_hull_centrif
                     gradcomp_embdsol_cube_nonli_embd_pressure
                     gradcomp_embdsol_cube_nonli_embd_centrif
                     gradcomp_embdsol_cube_nonli_hull_embd_pressure
                     gradcomp_embdsol_cube_nonli_hull_embd_centrif
                     vtu_output_2patchs
                     vtu_output_8patchs
                     vtu_output_2patchs_coupling
                     PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib/python:$ENV{PYTHONPATH})

# Properties for long tests
if(RUN_LONG_TESTS)
    set_tests_properties(OPTplatevib_semi_AN
                         OPTplatevib_FD
                         OPT3-volumetricBeam
                         PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib/python:$ENV{PYTHONPATH}
                         )
endif(RUN_LONG_TESTS)



