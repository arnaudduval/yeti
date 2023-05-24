! ====================================================
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! It is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

subroutine fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to steady heat problems
    !! by G. Sanaglli and M. Tani
    
    use solverheat
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigen_diag(nr_total), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call setup_preconditionerdiag(solv, nr_total, eigen_diag)
    call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)

end subroutine fd_steady_heat_3d

subroutine fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applieg in control points interpolation problems
    !! by G. Sanaglli and M. Tani
    
    use solverheat
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, array_in
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)

end subroutine fd_interpolation_3d

subroutine wq_find_capacity_diagonal_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, diag)
    !! Computes the diagonal of capacity matrix using sum-factorization algorithm
    !! IN CSR FORMAT

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    call csr_get_diag_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u(:, 1), data_B_v(:, 1), data_B_w(:, 1), &
                        data_W_u(:, 1), data_W_v(:, 1), data_W_w(:, 1), diag)

end subroutine wq_find_capacity_diagonal_3d

subroutine wq_find_conductivity_diagonal_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, diag)
    !! Computes the diagonal of conductivity matrix using sum-factorization algorithm
    !! IN CSR FORMAT

    implicit none
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)
    double precision :: diag_temp
    dimension :: diag_temp(nr_u*nr_v*nr_w)

    diag = 0.d0 
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            zeta = beta + (alpha - 1)*2
            call csr_get_diag_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                                data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), data_W_w(:, zeta(3)), diag_temp)
            diag = diag + diag_temp        
        end do
    end do

end subroutine wq_find_conductivity_diagonal_3d

subroutine mf_wq_get_cu_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(thermomat), pointer :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    allocate(mat)
    call setup_capacitycoefs(mat, nc_total, coefs)

    call mf_wq_capacity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_cu_3d

subroutine mf_wq_get_ku_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    allocate(mat)
    mat%dimen = 3
    call setup_conductivitycoefs(mat, nc_total, coefs)

    call mf_wq_conductivity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_in, array_out)
    
end subroutine mf_wq_get_ku_3d

subroutine mf_wq_get_spacetimeheat_3d(Ccoefs, Kcoefs, detG, nr_total, nc_sp, nc_time, nr_u, nc_u, nr_v, nc_v, &
                            nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, &
                            indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                            data_W_u, data_W_v, data_W_w, data_W_t, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_sp, nc_time, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    double precision, intent(in) :: Ccoefs, Kcoefs, detG
    dimension :: Ccoefs(nc_sp), Kcoefs(3, 3, nc_sp), detG(nc_time)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indi_T_t, indj_T_u, indj_T_v, indj_T_w, indj_T_t
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
    double precision :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)
    integer :: nc_total

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    call csr2csc(2, nr_t, nc_t, nnz_t, data_B_t, indj_t, indi_t, data_BT_t, indj_T_t, indi_T_t)
    
    allocate(mat)
    mat%dimen = 3
    call setup_capacitycoefs(mat, nc_sp, Ccoefs)
    call setup_conductivitycoefs(mat, nc_sp, Kcoefs)
    call setup_timediscret(mat, nc_time, detG)

    nc_total = nc_sp*nc_time
    call mf_wq_spacetimeheat_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t,indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, &
                            indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                            array_in, array_out)
    
end subroutine mf_wq_get_spacetimeheat_3d

subroutine wq_get_bodyheat_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(nr_u*nr_v*nr_w)

    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                        coefs, result)

end subroutine wq_get_bodyheat_3d

subroutine wq_get_heatflux_3d(coefs, JJ, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_out)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v
    double precision, intent(in) :: coefs, JJ
    dimension :: coefs(nc_total), JJ(d, d-1, nc_total)
    integer, intent(in) :: nnz_u, nnz_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: new_coefs, dsurf, v1, v2, v3
    dimension :: new_coefs(nc_total), v1(d), v2(d), v3(d)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total

        v1 = JJ(:, 1, i); v2 = JJ(:, 2, i)
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))
        new_coefs(i) = coefs(i) * dsurf

    end do

    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        new_coefs(:), array_out(:))

end subroutine wq_get_heatflux_3d

subroutine mf_wq_interpolate_cp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use matrixfreeheat
    use solverheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, Deigen
    double precision :: U_u, U_v, U_w
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w)

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
    
    ! Eigen decomposition
    allocate(Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w), &
            Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    call eigendecomp_heat_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, (/1.d0, 1.d0, 1.d0/), .false., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
    deallocate( Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, &
                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Set material and solver
    allocate(mat, solv)
    call setup_capacitycoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 1

    call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

end subroutine mf_wq_interpolate_cp_3d

subroutine mf_wq_steady_heat_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat
    use separatevariables

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in), target :: coefs
    dimension :: coefs(3, 3, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    type(operator), allocatable :: oper
    double precision :: kmean(3), kappa

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_total)
    double precision, dimension(:), allocatable :: Dphysical, Dparametric

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

    allocate(mat, solv)
    mat%dimen = 3
    call setup_conductivitycoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 2

    if (methodPCG.eq.'WP') then ! CG algorithm
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)
        
    else  

        ! Customize method
        allocate(Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        kmean = 1.d0

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'TDC')) then 
            call initialize_operator(oper, 3, (/nc_u, nc_v, nc_w/), (/.true., .true., .true./))
            call separatevariables_3d(oper, coefs)
            Mcoef_u = oper%MM(1, 1:nc_u); Mcoef_v = oper%MM(2, 1:nc_v); Mcoef_w = oper%MM(3, 1:nc_w)
            Kcoef_u = oper%KK(1, 1:nc_u); Kcoef_v = oper%KK(2, 1:nc_v); Kcoef_w = oper%KK(3, 1:nc_w)

        else if ((methodPCG.eq.'JMS').or.(methodPCG.eq.'JMC')) then 
            call compute_mean_3d(mat, nc_u, nc_v, nc_w)
            kmean = mat%mean(:3)
        
        end if

        ! Condition number 
        if ((methodPCG.eq.'C').or.(methodPCG.eq.'JMC')) then
            call compute_steadyheat_cond(nc_total, coefs, kmean, kappa)
            print*, 'Method: ', methodPCG, ', condition number: ', kappa

        end if

        ! Eigen decomposition
        allocate(Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
        call eigendecomp_heat_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, kmean, .true., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        call setup_preconditionerdiag(solv, nr_total, Deigen)

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then

            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, kmean, Dparametric)
            
            ! Find diagonal of real matrix
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)

            call setup_scaling(solv, nr_total, Dparametric, Dphysical)
        end if
        deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

        ! Set solver
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)
    end if

end subroutine mf_wq_steady_heat_3d

subroutine mf_wq_lineartransient_heat_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                b, thetadt, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreeheat
    use solverheat

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(3, 3, nc_total), Ccoefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG    
    double precision, intent(in) :: thetadt, threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision :: kmean(3), cmean, kappa

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_total)
    double precision, dimension(:), allocatable :: Dphysical, Dparametric, Dtemp

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

    allocate(mat, solv)
    mat%dimen = 3
    call setup_capacitycoefs(mat, nc_total, Ccoefs)
    call setup_conductivitycoefs(mat, nc_total, Kcoefs)
    mat%scalars = (/1.d0, thetadt/)
    solv%matrixfreetype = 3

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)
        
    else  

        ! Customize method  
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        kmean = 1.d0; cmean = 1.d0

        if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'JMS')) then 
            call compute_mean_3d(mat, nc_u, nc_v, nc_w)
            kmean = mat%mean(:3)
            cmean = mat%mean(4)
        
        end if

        ! Eigen decomposition
        allocate(Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
        call eigendecomp_heat_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, kmean, .true., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        Deigen = cmean + thetadt*Deigen
        call setup_preconditionerdiag(solv, nr_total, Deigen)

        if (methodPCG.eq.'JMS') then
            allocate(Dparametric(nr_total), Dphysical(nr_total), Dtemp(nr_total))

            ! Find diagonal of preconditioner
            Dtemp = 0.d0
            call kronvec3d(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, cmean)
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, kmean, Dparametric)
            Dparametric = Dtemp + thetadt*Dparametric

            ! Find diagonal of real matrix 
            Dtemp = 0.d0
            call wq_find_capacity_diagonal_3D(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dtemp)
            call wq_find_conductivity_diagonal_3D(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
            Dphysical = Dtemp + thetadt*Dphysical

            call setup_scaling(solv, nr_total, Dparametric, Dphysical)
        end if
        deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

        ! Condition number P^-1 A
        call compute_transientheat_cond(nc_total, Kcoefs, Ccoefs, kmean, cmean, kappa)
        if (methodPCG.eq.'JMS') then
            Dtemp = Dparametric/Dphysical
            kappa = kappa*maxval(Dtemp)/minval(Dtemp)
        end if
        print*, 'Method: ', methodPCG, ' condition number: ', kappa

        ! Set solver
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

    end if

end subroutine mf_wq_lineartransient_heat_3d
