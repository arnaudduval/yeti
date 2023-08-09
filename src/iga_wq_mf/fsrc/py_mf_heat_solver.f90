! ====================================================
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! It is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

subroutine mf_get_cu_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, isLumped, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v    
    integer, intent(in) :: indi_u, indi_v
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
    integer, intent(in) ::  indj_u, indj_v
    dimension :: indj_u(nnz_u), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    logical, intent(in) :: isLumped
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(thermomat) :: mat
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    mat%dimen = dimen; mat%isLumped = isLumped
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, prop)
    call mf_capacity_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, &
                    indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_in, array_out)

end subroutine mf_get_cu_2d

subroutine mf_get_cu_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isLumped, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w    
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    logical, intent(in) :: isLumped
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(thermomat) :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    mat%dimen = dimen; mat%isLumped = isLumped
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, prop)
    call mf_capacity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_get_cu_3d

subroutine mf_get_ku_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    
    mat%dimen = dimen
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_conductivity_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, &
                    indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_in, array_out)
    
end subroutine mf_get_ku_2d

subroutine mf_get_ku_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    mat%dimen = dimen
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_conductivity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_in, array_out)
    
end subroutine mf_get_ku_3d

subroutine get_heatvol_2d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, detJ, prop, array_out)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v
    integer, intent(in) :: nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: coefs(nc_total)

    coefs = prop*detJ
    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs, array_out)

end subroutine get_heatvol_2d

subroutine get_heatvol_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                        detJ, prop, array_out)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    double precision :: coefs(nc_total)

    coefs = prop*detJ
    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                        coefs, array_out)

end subroutine get_heatvol_3d

subroutine get_heatsurf_2d(nc_total, nr_u, nc_u, nnz_u, indi_u, indj_u, data_W_u, JJ, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_W_u
    dimension :: data_W_u(nnz_u, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, W00
    dimension :: coefs(nc_total), v1(dimen), W00(nr_u, nc_total)
    integer :: i

    if (nc_total.ne.nc_u) stop 'Size problem'

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        dsurf = sqrt(dot_product(v1, v1))
        coefs(i) = prop(i) * dsurf
    end do

    W00 = 0.d0
    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_total, W00)
    array_out = matmul(W00, prop)
    
end subroutine get_heatsurf_2d

subroutine get_heatsurf_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, JJ, prop, array_out)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! --------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v
    integer, intent(in) :: nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(nc_total), v1(dimen), v2(dimen), v3(dimen)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total

        v1 = JJ(:, 1, i); v2 = JJ(:, 2, i)
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))
        coefs(i) = prop(i) * dsurf

    end do

    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs(:), array_out(:))

end subroutine get_heatsurf_3d

subroutine solver_steady_heat_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, ndod, dod, table, &
                            invJ, detJ, prop, Fext, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat2
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    integer, intent(in) :: ndod 
    integer, intent(in) :: dod
    dimension :: dod(ndod)
    logical, intent(in) :: table
    dimension :: table(dimen, 2)

    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs

    ! Csr format
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    if (any(dod.le.0)) stop 'Indices must be greater than 0'

    mat%dimen = dimen
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v/)

    if (methodPCG.eq.'WP') then ! CG algorithm
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, &
                indj_u, indi_v, indj_v, data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
        
    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C').or.(methodPCG.eq.'TDC')) then
        if (methodPCG.eq.'JMC') then 
            call compute_mean(mat, nc_list)
        end if

        if (methodPCG.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(solv%temp_struct, size(univMcoefs, dim=1), size(univMcoefs, dim=2), &
                                        univMcoefs, univKcoefs)
        end if

        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, table, mat%Kmean(:dimen))

        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, &
                    indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
                
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_steady_heat_2d

subroutine solver_steady_heat_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            ndod, dod, table, invJ, detJ, prop, Fext, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat3
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: ndod
    integer, intent(in) :: dod
    dimension :: dod(ndod)
    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    mat%dimen = dimen
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v, nc_w/)

    if (methodPCG.eq.'WP') then ! CG algorithm
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)

    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C').or.(methodPCG.eq.'TDC')) then

        if (methodPCG.eq.'JMC') then 
            call compute_mean(mat, nc_list)
        end if

        if (methodPCG.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(solv%temp_struct, size(univMcoefs, dim=1), size(univMcoefs, dim=2), &
                                        univMcoefs, univKcoefs)
        end if

        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                                data_W_u, data_W_v, data_W_w, table, mat%Kmean(:dimen))
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)

    else 
        stop 'Unknown method' 
    end if

end subroutine solver_steady_heat_3d

subroutine solver_lineartransient_heat_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, isLumped, ndod, dod, table, &
                                invJ, detJ, Cprop, Kprop, thetadt, Fext, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat2

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    logical, intent(in) :: isLumped
    integer, intent(in) :: ndod
    integer, intent(in) :: dod
    dimension :: dod(ndod)
    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ, Cprop, Kprop, thetadt
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Cprop(nc_total), Kprop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs

    ! Csr format
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)
    
    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    mat%dimen = dimen; mat%isLumped = isLumped
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, Cprop)
    call setup_conductivityprop(mat, nc_total, Kprop)
    mat%scalars = (/1.d0, thetadt/); nc_list = (/nc_u, nc_v/)
    solv%matrixfreetype = 3

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
        
    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C').or.(methodPCG.eq.'TDC')) then

        if (methodPCG.eq.'JMC') then 
            call compute_mean(mat, nc_list)
        end if

        if (methodPCG.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(solv%temp_struct, size(univMcoefs, dim=1), size(univMcoefs, dim=2), &
                                        univMcoefs, univKcoefs)
        end if
    
        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                                data_W_u, data_W_v, table, mat%Kmean(:dimen))
        solv%temp_struct%Deigen = mat%Cmean + thetadt*solv%temp_struct%Deigen

        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_lineartransient_heat_2d

subroutine solver_lineartransient_heat_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isLumped, ndod, dod, table, &
                                invJ, detJ, Cprop, Kprop, thetadt, Fext, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat3

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    logical, intent(in) :: isLumped
    integer, intent(in) :: ndod
    integer, intent(in) :: dod
    dimension :: dod(ndod)
    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ, Cprop, Kprop, thetadt
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Cprop(nc_total), Kprop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    mat%dimen = dimen; mat%isLumped = isLumped
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, Cprop)
    call setup_conductivityprop(mat, nc_total, Kprop)
    mat%scalars = (/1.d0, thetadt/); nc_list = (/nc_u, nc_v, nc_w/)
    solv%matrixfreetype = 3

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
        
    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C').or.(methodPCG.eq.'TDC')) then

        if (methodPCG.eq.'JMC') then 
            call compute_mean(mat, nc_list)
        end if

        if (methodPCG.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(solv%temp_struct, size(univMcoefs, dim=1), size(univMcoefs, dim=2), &
                                        univMcoefs, univKcoefs)
        end if
    
        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                                data_W_u, data_W_v, data_W_w, table, mat%Kmean(:mat%dimen))
        solv%temp_struct%Deigen = mat%Cmean + thetadt*solv%temp_struct%Deigen
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, Fext, x, resPCG)
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_lineartransient_heat_3d
