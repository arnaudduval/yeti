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

!! The algorithm is taken from Piegl, Les. "The NURBS Book". Springer-Verlag:
!! Berlin 1995; p. 232-234.


!! Project a point at the surface of a solid
!! Parameters
!! point : coordinates of point to be projected
!! iface : index of solid face (yeti convention) on which projection is performed
!! coords3D : coordinates of control points (1st index : direction, 2nd index : pt index)
!! nb_cp : number of control points
!! is_hull : logical indicating if solid is a hulle of not
!! maxstep : max step allowed for parametric coordinates at each Newton iteration
!! maxiter : max number of iterations for Newton methode
!! u0 : initial value for 1st parametric coordinate
!! v0 : initial value for 2nd parametric coordinate
!! xi : return value : parametric coordinates of the projected point
!! info : return code :    0 = standard exit (point coincidence or zero cosine)
!! 1 = maximum number of iterations reached
!! 2 = exit because parameters do not change significantly
!! 3 = exit for any other reason
subroutine projection_surface(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter, u0, v0, &
        &   xi, info)
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
    double precision, intent(in) :: maxstep
    integer, intent(in) :: maxiter
    double precision, intent(in) :: u0, v0

    !! Output variable
    double precision, intent(out) :: xi
    dimension xi(3)
    integer, intent(out) :: info


    !! Local variables
    !! tolerance for zero consine condition
    double precision, parameter :: tol2 = 1.D-8  !! tol2 modif 21/03/2022
    !! tolerance for point coincidence condition
    double precision, parameter :: tol1 = 1.D-12
    double precision :: u, v, uprev, vprev, du, dv
    double precision :: S, Su, Sv, Suu, Svv, Suv
    dimension S(3), Su(3), Sv(3), Suu(3), Svv(3), Suv(3)
    double precision :: r, err1, err2
    dimension r(3)
    double precision :: Jmat, Kvec, Jinv, detJ
    dimension Jmat(2, 2), Kvec(2), Jinv(2, 2)
    integer :: i

    !! max correction value
    double precision :: step

    info = 3
    u = u0
    v = v0

    !! TODO : revoir les affectations pour etre plus rapide

    do i = 1, maxiter

        !! Evaluation at u,v
        if (is_hull) then
            call derivative_surface_mapping(u, v, iface, coords3D, nb_cp, S, Su, Sv, Suu, Suv, Svv)
        else
            call derivative_surface(u, v, iface, coords3D, nb_cp, S, Su, Sv, Suu, Suv, Svv)
        end if

        !! Check conditions
        ! 1. Point coincidence
        r(:) = S(:) - point(:)
        if (sqrt(dot_product(r, r)) <= tol1) then
            info = 0
            exit
        end if
        ! 2. Zero cosine
        err1 = abs(dot_product(Su, r)) / sqrt(dot_product(Su, Su) * dot_product(r, r))
        err2 = abs(dot_product(Sv, r)) / sqrt(dot_product(Sv, Sv) * dot_product(r, r))
        if ((err1 <= tol2) .and. (err2 <= tol2)) then
            info = 0
            exit
        end if

        Jmat(1, 1) = dot_product(Su, Su) + dot_product(r, Suu)
        Jmat(1, 2) = dot_product(Su, Sv) + dot_product(r, Suv)
        Jmat(2, 1) = dot_product(Su, sv) + dot_product(r, Suv)
        Jmat(2, 2) = dot_product(Sv, Sv) + dot_product(r, Svv)

        call MatrixInv(Jinv, Jmat, detJ, 2)

        Kvec(1) = - dot_product(r, Su)
        Kvec(2) = - dot_product(r, Sv)

        uprev = u
        vprev = v

        du = Jinv(1, 1) * Kvec(1) + Jinv(1, 2) * Kvec(2)
        dv = Jinv(2, 1) * Kvec(1) + Jinv(2, 2) * Kvec(2)

        !! Correction damping
        step = sqrt(du * du + dv * dv)
        if (step >= maxStep) then
            du = du * maxStep / step
            dv = dv * maxStep / step
        end if

        u = u + du
        v = v + dv

        select case (iface)
            case (1, 2)
            if ( u < minval(Ukv2_patch)) u = minval(Ukv2_patch)
            if ( u > maxval(Ukv2_patch)) u = maxval(Ukv2_patch)
            if ( v < minval(Ukv3_patch)) v = minval(Ukv3_patch)
            if ( v > maxval(Ukv3_patch)) v = maxval(Ukv3_patch)
            case (3, 4)
            if ( u < minval(Ukv1_patch)) u = minval(Ukv1_patch)
            if ( u > maxval(Ukv1_patch)) u = maxval(Ukv1_patch)
            if ( v < minval(Ukv3_patch)) v = minval(Ukv3_patch)
            if ( v > maxval(Ukv3_patch)) v = maxval(Ukv3_patch)
            case (5, 6)
            if ( u < minval(Ukv1_patch)) u = minval(Ukv1_patch)
            if ( u > maxval(Ukv1_patch)) u = maxval(Ukv1_patch)
            if ( v < minval(Ukv2_patch)) v = minval(Ukv2_patch)
            if ( v > maxval(Ukv2_patch)) v = maxval(Ukv2_patch)
        end select

        if (norm2((u - uprev) * Su + (v - vprev) * Sv) <= tol1) then
            info = 2
            exit
        end if

    end do

    if (i == maxiter) then
        write(*, *) "Warning : max number of iterations reached during point to surface projection"
        info = 1
    end if

    if (is_hull) then
        call point_on_solid_face_map(u, v, iface, xi)
    else
        call point_on_solid_face(u, v, iface, xi)
    end if

end subroutine projection_surface

!! Compute 3D parameters of a point on a given face
subroutine point_on_solid_face(u, v, iface, xi)
    use nurbspatch

    implicit none

    double precision, intent(in) :: u, v
    integer, intent(in) :: iface
    double precision, intent(out) :: xi
    dimension xi(3)

    select case (iface)
        case (1)
        xi(1) = minval(ukv1_patch)
        xi(2) = u
        xi(3) = v
        case (2)
        xi(1) = maxval(ukv1_patch)
        xi(2) = u
        xi(3) = v
        case (3)
        xi(1) = u
        xi(2) = minval(ukv2_patch)
        xi(3) = v
        case (4)
        xi(1) = u
        xi(2) = maxval(ukv2_patch)
        xi(3) = v
        case (5)
        xi(1) = u
        xi(2) = v
        xi(3) = minval(ukv3_patch)
        case (6)
        xi(1) = u
        xi(2) = v
        xi(3) = maxval(ukv3_patch)
    end select

end subroutine point_on_solid_face

!! Compute 3D parameters of a point on a given mapping face
subroutine point_on_solid_face_map(u, v, iface, xi)
    use nurbspatch
    use embeddedmapping

    implicit none

    double precision, intent(in) :: u, v
    integer, intent(in) :: iface
    double precision, intent(out) :: xi
    dimension xi(3)

    select case (iface)
        case (1)
        xi(1) = minval(ukv1_map)
        xi(2) = u
        xi(3) = v
        case (2)
        xi(1) = maxval(ukv1_map)
        xi(2) = u
        xi(3) = v
        case (3)
        xi(1) = u
        xi(2) = minval(ukv2_map)
        xi(3) = v
        case (4)
        xi(1) = u
        xi(2) = maxval(ukv2_map)
        xi(3) = v
        case (5)
        xi(1) = u
        xi(2) = v
        xi(3) = minval(ukv3_map)
        case (6)
        xi(1) = u
        xi(2) = v
        xi(3) = maxval(ukv3_map)
    end select

end subroutine point_on_solid_face_map

!! Compute derivative surface of a solid
!! Parameters
!! u : 1st parametric coordinate at which derivatives are computed
!! v : 2nd parametric coordinate at which derivatives are computed
!! iface : index of solid face (yeti convention) on which projection is performed
!! coords3D : coordinates of control points (1st index : direction, 2nd index : pt index)
!! nb_cp : number of control points
!! S : return value : surface
!! Su : return value : derivative of surface w.r.t u
!! Sv : return value : derivative of surface w.r.t v
!! Suu : return value : 2nd derivative of surface w.r.t u
!! Suv : return value : 2nd derivative of surface w.r.t u and v
!! Svv : return value : 2nd derivative of surface w.r.t v
subroutine derivative_surface(u, v, iface, coords3D, nb_cp, S, Su, Sv, Suu, Suv, Svv)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: u, v
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp


    !! Output variables
    double precision, intent(out) :: S, Su, Sv, Suu, Suv, Svv
    dimension S(3), Su(3), Sv(3), Suu(3), Suv(3), Svv(3)

    !! Local variables
    integer :: icp, idim
    double precision :: xi
    dimension xi(3)
    integer :: Ni
    dimension Ni(dim_patch)
    double precision :: coords
    dimension coords(3, nnode_patch)
    double precision :: R, dRdxi, ddRddxi
    dimension R(nnode_patch), dRdxi(nnode_patch, 3), ddRddxi(nnode_patch, 6)


    call point_on_solid_face(u, v, iface, xi)

    !! Get knot span and CP coordinates
    call updateElementNumber(xi)
    do icp = 1, nnode_patch
        coords(:, icp) = COORDS3D(:3, IEN_patch(icp, current_elem))
    end do

    call evalnurbs_w2ndDerv(xi, R, dRdxi, ddRddxi)

    S(:) = zero
    Su(:) = zero
    Sv(:) = zero
    Suu(:) = zero
    Suv(:) = zero
    Svv(:) = zero

    do icp = 1, nnode_patch
        S(:) = S(:) + R(icp) * coords(:, icp)
        select case (iface)
            !! TODO, on peut fair une table de correspondance avec les 2ème indices
            !! de dRdxi et ddRddxi
            case (1, 2)   !! xi(1) const
            Su(:) = Su(:) + dRdxi(icp, 2) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 3) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 2) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 3) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 6) * coords(:, icp)
            case (3, 4)   !! xi(2) const
            Su(:) = Su(:) + dRdxi(icp, 1) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 3) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 1) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 3) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 5) * coords(:, icp)
            case (5, 6)   !! xi(3) const
            Su(:) = Su(:) + dRdxi(icp, 1) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 2) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 1) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 2) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 4) * coords(:, icp)
        end select
    end do

end subroutine derivative_surface



!! Compute derivative surface of a solid mapping
subroutine derivative_surface_mapping(u, v, iface, coords3D, nb_cp, S, Su, Sv, Suu, Suv, Svv)
    use parameters
    use nurbspatch
    use embeddedmapping

    implicit none

    !! Input variables
    double precision, intent(in) :: u, v
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp


    !! Output variables
    double precision, intent(out) :: S, Su, Sv, Suu, Suv, Svv
    dimension S(3), Su(3), Sv(3), Suu(3), Suv(3), Svv(3)

    !! Local variables
    integer :: icp, idim
    double precision :: xi
    dimension xi(3)
    double precision :: coords
    dimension coords(3, nnode_map)
    double precision :: R, dRdxi, ddRddxi
    dimension R(nnode_map), dRdxi(nnode_map, 3), ddRddxi(nnode_map, 6)


    call point_on_solid_face_map(u, v, iface, xi)

    !! Get knot span and CP coordinates
    call updateMapElementNumber(xi)
    do icp = 1, nnode_map
        coords(:, icp) = COORDS3D(:3, IEN_map(icp, current_map_elem))
    end do

    call evalnurbs_mapping_w2ndDerv(xi, R, dRdxi, ddRddxi)

    S(:) = zero
    Su(:) = zero
    Sv(:) = zero
    Suu(:) = zero
    Suv(:) = zero
    Svv(:) = zero

    do icp = 1, nnode_map
        S(:) = S(:) + R(icp) * coords(:, icp)
        select case (iface)
            !! TODO, on peut faire une table de correspondance avec les 2ème indices
            !! de dRdxi et ddRddxi
            case (1, 2)   !! xi(1) const
            Su(:) = Su(:) + dRdxi(icp, 2) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 3) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 2) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 3) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 6) * coords(:, icp)
            case (3, 4)   !! xi(2) const
            Su(:) = Su(:) + dRdxi(icp, 1) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 3) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 1) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 3) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 5) * coords(:, icp)
            case (5, 6)   !! xi(3) const
            Su(:) = Su(:) + dRdxi(icp, 1) * coords(:, icp)
            Sv(:) = Sv(:) + dRdxi(icp, 2) * coords(:, icp)
            Suu(:) = Suu(:) + ddRddxi(icp, 1) * coords(:, icp)
            Svv(:) = Svv(:) + ddRddxi(icp, 2) * coords(:, icp)
            Suv(:) = Suv(:) + ddRddxi(icp, 4) * coords(:, icp)
        end select
    end do

end subroutine derivative_surface_mapping


!! Project a point at the limit of a 2D solid
!! Parameters
!! point : coordinates of point to be projected
!! iface : index of solid face (yeti convention) on which projection is performed
!! coords3D : coordinates of control points (1st index : direction, 2nd index : pt index)
!! nb_cp : number of control points
!! is_hull : logical indicating if solid is a hulle of not
!! maxstep : max step allowed for parametric coordinates at each Newton iteration
!! maxiter : max number of iterations for Newton methode
!! u0 : initial value for 1st parametric coordinate
!! xi : return value : parametric coordinates of the projected point
!! info : return code :    0 = standard exit (point coincidence or zero cosine)
!! 1 = maximum number of iterations reached
!! 2 = exit because parameters do not change significantly
!! 3 = exit for any other reason
subroutine projection_curve(point, iface, coords3D, nb_cp, is_hull, maxstep, maxiter, u0, &
        &   xi, info)

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
    double precision, intent(in) :: maxstep
    integer, intent(in) :: maxiter
    double precision, intent(in) :: u0

    !! Output variable
    double precision, intent(out) :: xi
    dimension xi(2)
    integer, intent(out) :: info

    !! Local variables
    ! Tolerance for point coincidence condition
    double precision, parameter :: tol1 = 1.D-8
    ! Tolerance for zero cosine condition
    double precision, parameter :: tol2 = 1.D-12
    double precision :: u, uprev
    double precision, dimension(3) :: C
    double precision, dimension(3) :: Cu
    double precision, dimension(3) :: Cuu
    double precision, dimension(3) :: r
    double precision :: err
    integer :: i

    info = 3
    u = u0

    do i = 1, maxiter
        !! Evaluation at u
        call derivative_curve(u, iface, coords3D, nb_cp, C, Cu, Cuu)

        !! Check condition
        ! 1. Point coincidence
        r(:) = C(:) - point(:)
        if (sqrt(dot_product(r, r)) <= tol1) then
            info = 0
            exit
        end if
        ! 2. Zero cosine
        err = abs(dot_product(Cu, r)) / sqrt(dot_product(Cu, Cu) * dot_product(r, r))
        if ((err <= tol2) .and. (i /= 1)) then
            info = 4
            exit
        end if

        uprev = u
        u = u - dot_product(Cu, r) / (dot_product(Cuu, r) + dot_product(Cu, Cu))

        select case (iface)
            case (1, 2)
            if ( u < minval(Ukv2_patch)) u = minval(Ukv2_patch)
            if ( u > maxval(Ukv2_patch)) u = maxval(Ukv2_patch)
            case (3, 4)
            if ( u < minval(Ukv1_patch)) u = minval(Ukv1_patch)
            if ( u > maxval(Ukv1_patch)) u = maxval(Ukv1_patch)
        end select

        if (norm2((u - uprev) * Cu) <= tol1) then
            info = 2
            exit
        end if
    end do

    if (i == maxiter) then
        write(*, *) "Warning : max number of iterations reached during point to surface projection"
        info = 1
    end if

    call point_on_2D_solid_curve(u, iface, xi)

end subroutine projection_curve

!! Compute derivative curve of a 2D solid
!! Parameters
!! u : 1st parametric coordinate at which derivatives are computed
!! iface : index of solid face (yeti convention) on which projection is performed
!! coords3D : coordinates of control points (1st index : direction, 2nd index : pt index)
!! nb_cp : number of control points
!! C : return value : curve
!! Cu : return value : derivative of curve w.r.t u
!! Cuu : return value : 2nd derivative of curve w.r.t u
subroutine derivative_curve(u, iface, coords3D, nb_cp, C, Cu, Cuu)
    use parameters
    use nurbspatch

    implicit none

    !! Input variables
    double precision, intent(in) :: u
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    integer, intent(in) :: nb_cp

    !! Output variables
    double precision, intent(out) :: C, Cu, Cuu
    dimension C(3), Cu(3), Cuu(3)

    !! Local variables
    integer :: icp
    double precision, dimension(3) :: xi
    double precision, dimension(3, nnode_patch) :: coords
    double precision, dimension(nnode_patch) :: R
    double precision, dimension(nnode_patch, 3) :: dRdxi
    double precision, dimension(nnode_patch, 6) :: ddRddxi

    call point_on_2D_solid_curve(u, iface, xi)

    !! Get knot span and CP coordinates
    call updateElementNumber(xi)
    do icp = 1, nnode_patch
        coords(:, icp) = COORDS3D(:3, IEN_patch(icp, current_elem))
    end do

    call evalnurbs_w2ndDerv(xi, R, dRdxi, ddRddxi)

    C(:) = zero
    Cu(:) = zero
    Cuu(:) = zero

    do icp = 1, nnode_patch
        C(:) = C(:) + R(icp) * coords(:, icp)
        select case (iface)
            !! TODO, on peut fair une table de correspondance avec les 2ème indices
            !! de dRdxi et ddRddxi
            case (1, 2)   !! xi(1) const
            Cu(:) = Cu(:) + dRdxi(icp, 2) * coords(:, icp)
            Cuu(:) = Cuu(:) + ddRddxi(icp, 2) * coords(:, icp)
            case (3, 4)   !! xi(2) const
            Cu(:) = Cu(:) + dRdxi(icp, 1) * coords(:, icp)
            Cuu(:) = Cuu(:) + ddRddxi(icp, 1) * coords(:, icp)
        end select
    end do

end subroutine derivative_curve

!! Compute 2D parameters of a point on a given edge of a 2D solid
subroutine point_on_2D_solid_curve(u, iface, xi)
    use nurbspatch

    implicit none

    double precision, intent(in) :: u
    integer, intent(in) :: iface
    double precision, intent(out) :: xi
    dimension xi(2)

    select case (iface)
        case (1)
        xi(1) = minval(ukv1_patch)
        xi(2) = u
        case (2)
        xi(1) = maxval(ukv1_patch)
        xi(2) = u
        case (3)
        xi(1) = u
        xi(2) = minval(ukv2_patch)
        case (4)
        xi(1) = u
        xi(2) = maxval(ukv2_patch)
    end select

end subroutine point_on_2D_solid_curve

