!! Perform trivarate point inversion
!! Referenc: https://doi.org/10.1007/s12008-020-00719-z
!! exit status :
!!  - 0 Point coincidence
!!  - 1 Zero cosine
!!  - 2 Non evolution
!!  - 3 max number of iteration reached
!! conergence parameters :
!!  - eps1 : point coincidence tolerance. Depends on geometry size
!!  - eps2 : zero cosine tolerance
!!  - eps3 : non evolution tolerance

subroutine trivariate_point_inversion(param, status, distance, point, i_patch, &
    & eps1, eps2, eps3, maxiter, nnode, ien,      &
    & nb_elem_patch, nkv, jpqr, nijk, ukv, weight, coords, nb_patch, &
    & nb_elem, nb_cp, mcrd, len_ukv, len_ien, len_weight)

    use, intrinsic :: iso_c_binding

    use nurbspatch
    use linalg3x3

    implicit none

    real(c_double), dimension(3), intent(out) :: param
    integer(c_int), intent(out) :: status
    real(c_double), intent(out) :: distance

    real(c_double), dimension(mcrd) , intent(in):: point
    integer(c_int), intent(in) :: i_patch
    integer(c_int), intent(in) :: maxiter
    real(c_double), intent(in) :: eps1
    real(c_double), intent(in) :: eps2
    real(c_double), intent(in) :: eps3

    integer(c_int), dimension(nb_patch) :: nnode
    integer(c_int), dimension(len_ien), intent(in) :: ien
    integer(c_int), dimension(nb_patch), intent(in) :: nb_elem_patch
    integer(c_int), dimension(3, nb_patch), intent(in) :: nkv
    integer(c_int), dimension(3, nb_patch), intent(in) :: jpqr
    integer(c_int), dimension(3, nb_elem), intent(in) :: nijk
    real(c_double), dimension(len_ukv), intent(in) :: ukv
    real(c_double), dimension(len_weight), intent(in) :: weight
    real(c_double), dimension(3, nb_cp), intent(in) :: coords
    integer(c_int), intent(in) :: nb_patch
    integer(c_int), intent(in) :: nb_elem
    integer(c_int), intent(in) :: nb_cp
    integer(c_int), intent(in) :: mcrd
    integer(c_int), intent(in) :: len_ukv
    integer(c_int), intent(in) :: len_ien
    integer(c_int), intent(in) :: len_weight

    !! NURBS functions
    real(c_double), allocatable, dimension(:) :: R
    !! derivatives of NURBS functions
    real(c_double), allocatable, dimension(:, :) :: dRdxi
    !! 2nd derivatives of NURBS function
    !! convention: uu, vv, ww, uv, uw, vw
    real(c_double), allocatable, dimension(:, :) :: ddRddxi
    real(c_double), allocatable, dimension(:, :) :: coords_elem
    real(c_double), dimension(3) :: V, residual, delta, fgh
    real(c_double), dimension(3, 3) :: dVdu
    real(c_double), dimension(3, 6) :: ddVddu
    real(c_double), dimension(3, 3) :: jac
    logical(c_bool) :: converged
    integer(c_int) :: iter, i, j
    logical :: zero_cosines
    real(c_double) :: cosval

    integer, parameter :: kmap(3,3) = reshape([ &
    1, 4, 5, &
    4, 2, 6, &
    5, 6, 3 ], [3,3])

    !! fake input data ==> USE real data from IGAParametrization ???
    integer(c_int), dimension(1) :: jprops
    real(c_double), dimension(1) :: props
    jprops(:) = 0
    props(:) = 0.

    !! Extract data for patch of interest
    call extractNurbsPatchGeoInfos(i_patch, nkv, jpqr, nijk, ukv, &
        &   weight, nb_elem_patch)
    call extractNurbsPatchMechInfos(i_patch, ien, PROPS, JPROPS,  &
        &        nnode, nb_elem_patch, 'NONE', 'NONE')
    allocate(R(nnode_patch), dRdxi(nnode_patch, 3), ddRddxi(nnode_patch, 6))
    allocate(coords_elem(mcrd, nnode_patch))


    !! Initialization
    param(:) = 0.5
    iter = 1
    converged = .false.


    do while ((.not.converged) .and. (iter <= maxiter))
        call updateElementNumber(param)
        call evalnurbs_w2ndDerv(param, R, dRdxi, ddRddxi)

        do i = 1, nnode_patch
            coords_elem(:, i) = coords(:mcrd, IEN_patch(i, current_elem))
        enddo

        !! Compute V and its derivative and second derivative
        V(:) = 0.0
        dVdu(:, :) = 0.0
        ddVddu(:, :) = 0.0
        do i = 1, nnode_patch
            V(:) = V(:) + R(i)*coords_elem(:, i)
        enddo

        residual = V - point

        !! Test point coincidence
        distance = norm2(residual)
        ! write(*,*) param, '-->', distance

        if (distance < eps1) then
            status = 0
            exit
        endif

        do i = 1, nnode_patch
            do j = 1, 3
                dVdu(:, j) = dVdu(:, j) + dRdxi(i, j)*coords_elem(:, i)
            enddo
            do j = 1, 6
                ddVddu(:, j) = ddVddu(:, j) + ddRddxi(i, j)*coords_elem(:, i)
            enddo
        enddo

        ! Test zero cosine
        if (distance > 0.0_c_double) then
            zero_cosines = .true.
            do i = 1, 3
                if (norm2(dVdu(:,i)) <= 0.0_c_double) then
                    zero_cosines = .false.
                    exit
                end if
                cosval = abs(dot_product(dVdu(:,i), residual)) / (norm2(dVdu(:,i)) * distance)
                if (cosval >= eps2) then
                    zero_cosines = .false.
                    exit
                end if
            end do
            if (zero_cosines) then
                status = 1
                converged = .true.
            end if
        end if

        !! compute correction
        do i = 1, 3
            do j = 1, 3
                jac(i,j) = dot_product(dVdu(:,i), dVdu(:,j)) + dot_product(residual, ddVddu(:, kmap(i, j)))
            end do
        end do

        do i=1, 3
            fgh(i) = dot_product(residual, dVdu(:, i))
        enddo

        delta = - solve3x3(jac, fgh)
        !! Test non-evolution
        if (norm2(matmul(dVdu, delta)) < eps3) then
            status = 2
            converged = .true.
        endif

        param = param + delta

        param = max(0.0, min(1.0, param))


        iter = iter + 1
    enddo

    if (iter == maxiter+1) then
        status = 3
    endif

    !! Recompute distance
    call updateElementNumber(param)
    call evalnurbs_noder(param, R)

    do i = 1, nnode_patch
        coords_elem(:, i) = coords(:mcrd, IEN_patch(i, current_elem))
    enddo

    !! Compute V
    V(:) = 0.0
    do i = 1, nnode_patch
        V(:) = V(:) + R(i)*coords_elem(:, i)
    enddo

    residual = V - point
    distance = norm2(residual)

    deallocate(R, dRdxi, ddRddxi)
    call finalizeNurbsPatch()
end subroutine trivariate_point_inversion

