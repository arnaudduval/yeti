!! Copyright 2021 Arnaud Duval

!! This file is part of Yeti.
!!
!! Yeti is free software: you can redistribute it and/or modify it under the terms
!! of the GNU Lesser General Public License as published by the Free Software
!! Foundation, either version 3 of the License, or (at your option) any later version.
!!
!! Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!! PURPOSE. See the GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License along
!! with Yeti. If not, see <https://www.gnu.org/licenses/>


!! Get parametric coordinates from physical ones on the surface of a solid
!! First try to inverse point with projection, then use dichotomy search

!! Parameters
!! point : physical coordinates of point
!! iface : index of face
!! coords3D : control points coordinates
!! nb_cp : number of control points
!! isHull : boolean indicating if projection is made on a hull
!! xi : evaluated parametric coordinates
!! info : return code :    0 = standard exit
subroutine point_inversion_surface(point, iface, coords3D, nb_cp, is_hull, xi, info)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: point
    dimension point(3)
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp
    logical :: is_hull

    !! Output variables
    double precision, intent(out) :: xi
    dimension xi(3)
    integer, intent(out) :: info

    !! local variables
    integer :: projection_info
    integer :: maxiter
    double precision :: maxstep
    double precision :: u0, v0
    double precision :: S, d, dist
    dimension S(3), d(3), dist(4)
    integer :: icp, iquad
    double precision :: R
    dimension R(nnode_patch)
    double precision :: u_dicho, v_dicho, umin, umax, umed, vmin, vmax, vmed
    integer :: iter
    !!
    double precision :: lim_u, lim_v
    dimension lim_u(2), lim_v(2)
    double precision :: step_u, step_v
    integer :: num_u, num_v

    info = 0
    maxiter = 100
    maxstep = 0.1D0
    num_u = 15
    num_v = 15

    !! 1. Use default start values: u0 = 0.5, v0 = 0.5
    u0 = 0.5
    v0 = 0.5
    call projection_surface(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter,  &
        &  u0, v0, xi, projection_info)

    if (projection_info /= 0) then
        !! 2. Use sampling to define starting parameters values
        lim_u = (/ 0.0d0, 1.0d0 /)
        lim_v = (/ 0.0d0, 1.0d0 /)
        call compute_starting_point(point, iface, coords3D, nb_cp, is_hull, num_u, num_v, &
            &   lim_u, lim_v, u0, v0)
        call projection_surface(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter,  &
            &  u0, v0, xi, projection_info)

        if (projection_info /= 0) then
            !! 3. Use sampling to define starting parameters values, with more sampling points
            lim_u = (/ 0.0d0, 1.0d0 /)
            lim_v = (/ 0.0d0, 1.0d0 /)
            num_u = num_u * 2
            num_v = num_v * 2
            call compute_starting_point(point, iface, coords3D, nb_cp, is_hull, num_u, num_v, &
                &   lim_u, lim_v, u0, v0)
            call projection_surface(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter,  &
                &  u0, v0, xi, projection_info)

            if (projection_info /= 0) then
                !! 4. Use sampling again, centering on previously computed parameters values
                ! Compute previous step value
                step_u = (lim_u(2) - lim_u(1)) / (1.0d0 * num_u)
                step_v = (lim_v(2) - lim_v(1)) / (1.0d0 * num_v)
                ! Update limits & ensure it stays between 0 and 1
                lim_u(1) = u0 - step_u
                if (lim_u(1) < 0.0d0) lim_u(1) = 0.0d0
                lim_u(2) = u0 + step_u
                if (lim_u(2) > 1.0d0) lim_u(2) = 1.0d0
                lim_v(1) = v0 - step_v
                if (lim_v(1) < 0.0d0) lim_v(1) = 0.0d0
                lim_v(2) = v0 + step_v
                if (lim_v(2) > 1.0d0) lim_v(2) = 1.0d0
                ! Compute params. & run proj.
                call compute_starting_point(point, iface, coords3D, nb_cp, is_hull, num_u, num_v, &
                    &   lim_u, lim_v, u0, v0)
                call projection_surface(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter,  &
                    &  u0, v0, xi, projection_info)
                if (projection_info /= 0) then
                    !! 5. Reduce step & increase max. nb. of iterations
                    do
                        call projection_surface(point, iface, coords3D, nb_cp, is_hull, &
                            &   maxstep, maxiter, u0, v0, xi, projection_info)
                        if (projection_info == 0) return
                        maxstep = maxstep / 10.D0
                        maxiter = maxiter * 10
                        if (maxiter > 10000) then
                            !! Just giving up
                            exit
                        end if
                    end do
                end if
            end if
        end if
    end if

    info = projection_info

end subroutine point_inversion_surface


!! Compute starting parameter values
subroutine compute_starting_point(point, iface, coords3D, nb_cp, is_hull, num_u, num_v, &
        &   lim_u, lim_v, u, v)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: point
    dimension point(3)
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp
    logical :: is_hull
    integer, intent(in) :: num_u, num_v
    double precision, intent(in) :: lim_u, lim_v
    dimension lim_u(2), lim_v(2)

    !! Output variable
    double precision, intent(out) :: u, v


    !! Local variables
    integer :: i_u, i_v, icp
    double precision :: param_u, param_v, uvw, N, eval_pt, coords, &
        &   diff, dist, idx, step_u, step_v
    dimension uvw(3), N(nnode_patch), eval_pt(3), coords(3, nnode_patch), &
        &   diff(3), dist(num_u + 2, num_v + 2), idx(2)


    !! Initialisations
    dist(:, :) = zero
    step_u = (lim_u(2) - lim_u(1)) / (1.0d0 * num_u)
    step_v = (lim_v(2) - lim_v(1)) / (1.0d0 * num_v)

    !! Evaluate equally spaced values & compute distance to point
    do i_v = 1, num_v + 2
        param_v = lim_v(1) + (i_v - 1) * 1.0d0 * step_v
        do i_u = 1, num_u + 2
            param_u = lim_u(1) + (i_u - 1) * 1.0d0 * step_u
            if (is_hull) then    !! Hull object
                call evaluate_surf_and_compute_diff_hull(point, iface, coords3D, &
                    &   nb_cp, param_u, param_v, &
                    &   diff)
            else    !! Embedded entity or classical solid
                call evaluate_surf_and_compute_diff(point, iface, coords3D, &
                    &   nb_cp, param_u, param_v, &
                    &   diff)
            end if
            !! Compute distance & add to list
            dist(i_u, i_v) = sqrt(dot_product(diff, diff))
        end do
    end do

    !! Locate minimal distance
    idx = minloc(dist)

    !! Compute starting parameters values
    u = lim_u(1) + (idx(1) - 1) * 1.0d0 * step_u
    v = lim_v(1) + (idx(2) - 1) * 1.0d0 * step_v

end subroutine compute_starting_point


!! Evaluate equally spaced values & compute distance to point
subroutine evaluate_surf_and_compute_diff(point, iface, coords3D, nb_cp, param_u, param_v, &
        &   diff)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: point
    dimension point(3)
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp
    double precision :: param_u, param_v

    !! Output variables
    double precision, intent(out) :: diff
    dimension diff(3)

    !! Local variables
    integer :: icp
    double precision :: uvw
    dimension uvw(3)
    double precision :: coords
    dimension coords(3, nnode_patch)
    double precision :: R
    dimension R(nnode_patch)
    double precision :: eval_pt
    dimension eval_pt(3)


    !! Define solid parameters
    call point_on_solid_face(param_u, param_v, iface, uvw)
    !! Get knot span and CP coordinates
    call updateElementNumber(uvw)
    do icp = 1, nnode_patch
        coords(:, icp) = COORDS3D(:3, IEN_patch(icp, current_elem))
    end do
    !! Evaluate functions
    call evalnurbs_noder(uvw, R)
    !! Evaluate surface
    eval_pt(:) = zero
    do icp = 1, nnode_patch
        eval_pt(:) = eval_pt(:) + R(icp) * coords(:, icp)
    end do
    !! Compute distance to point
    diff = eval_pt(:) - point(:)

end subroutine

!! Evaluate equally spaced values & compute distance to point (hull)
subroutine evaluate_surf_and_compute_diff_hull(point, iface, coords3D, nb_cp, param_u, param_v, &
        &   diff)
    use parameters
    use nurbspatch
    use embeddedmapping

    implicit none

    !! Input variables
    double precision, intent(in) :: point
    dimension point(3)
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp
    double precision :: param_u, param_v

    !! Output variables
    double precision, intent(out) :: diff
    dimension diff(3)

    !! Local variables
    integer :: icp
    double precision :: uvw
    dimension uvw(3)
    double precision :: coords
    dimension coords(3, nnode_map)
    double precision :: R
    dimension R(nnode_map)
    double precision :: eval_pt
    dimension eval_pt(3)


    !! Define solid parameters
    call point_on_solid_face_map(param_u, param_v, iface, uvw)
    !! Get knot span and CP coordinates
    call updateMapElementNumber(uvw)
    do icp = 1, nnode_map
        coords(:, icp) = COORDS3D(:3, IEN_map(icp, current_map_elem))
    end do
    !! Evaluate functions
    call evalnurbs_mapping_noder(uvw, R)
    !! Evaluate surface
    eval_pt(:) = zero
    do icp = 1, nnode_map
        eval_pt(:) = eval_pt(:) + R(icp) * coords(:, icp)
    end do
    !! Compute distance to point
    diff = eval_pt(:) - point(:)

end subroutine


!! Point inversion on a plane curve
!! Parameters
!! point : physical coordinates of point
!! iface : index of face
!! coords3D : control points coordinates
!! nb_cp : number of control points
!! isHull : boolean indicating if projection is made on a hull
!! xi : evaluated parametric coordinates
!! info : return code :    0 = standard exit
subroutine point_inversion_plane_curve(point, iface, coords3D, nb_cp, is_hull, xi, info)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: point
    dimension point(3)
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp
    logical :: is_hull

    !! Output variables
    double precision, intent(out) :: xi
    dimension xi(2)
    integer, intent(out) :: info

    !! Local variables
    integer :: maxiter
    integer :: maxstep
    double precision :: u0
    integer :: projection_info


    maxiter = 100
    maxstep = 0.1D0

    u0 = 0.5

    call projection_curve(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter,  &
        &  u0, xi, projection_info)

    info = projection_info

end subroutine point_inversion_plane_curve
