! ====================================================
! module :: Matrix free methods in WQ approach 
! author :: Joaquin Cornejo
! 
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! It is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

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

! ---------------------------------
! Matrix free multiplication in 3D
! ---------------------------------
subroutine mf_wq_get_cu_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
    !! Computes C.u where C is the capacity matrix in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    
    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    double precision :: array_temp_0, array_temp_1
    dimension :: array_temp_0(nc_total), array_temp_1(nc_total)

    ! Compute B' * array_in
    call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                        array_in, array_temp_0)

    ! Compute diag(coefs) * array_temp_0
    array_temp_1 = array_temp_0 * coefs

    ! Compute W * array_temp_1
    call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                        array_temp_1, array_out)
    
end subroutine mf_wq_get_cu_3d

subroutine mf_wq_get_cu_3d_csr(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! ---------- 
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_wq_get_cu_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_cu_3d_csr

subroutine mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! -----------
    double precision :: array_temp_0, array_temp_1, array_temp_2
    dimension :: array_temp_0(nc_total), array_temp_1(nc_total), array_temp_2(nc_total)
    integer :: i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)

    array_out = 0.d0
    do j = 1, d
        beta = 1; beta(j) = 2
        call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                            array_in, array_temp_0)
        do i = 1, d
            alpha = 1; alpha(i) = 2
            zeta = beta + (alpha - 1)*2
            array_temp_1 = array_temp_0*coefs(i, j, :)
            call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                array_temp_1, array_temp_2)
            array_out = array_out + array_temp_2
        end do
    end do
    
end subroutine mf_wq_get_ku_3d

subroutine mf_wq_get_ku_3d_csr(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w,&
                        data_W_u, data_W_v, data_W_w, array_in, array_out)
    
end subroutine mf_wq_get_ku_3d_csr

subroutine mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, alpha, beta, array_in, array_out)
    !! Computes (alpha*C + beta*K).u where C and K are capacity and conductivity matrices respectively in 3D case
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: Kcoefs, Ccoefs, alpha, beta
    dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data
    ! ---------------
    double precision :: array_temp_0, array_temp_1
    dimension :: array_temp_0(nr_u*nr_v*nr_w), array_temp_1(nr_u*nr_v*nr_w)

    call mf_wq_get_cu_3d(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_temp_0)

    call mf_wq_get_ku_3d(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_temp_1)

    array_out = alpha*array_temp_0 + beta*array_temp_1
    
end subroutine mf_wq_get_kcu_3d

subroutine mf_wq_get_kcu_3d_csr(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, alpha, beta, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: Kcoefs, Ccoefs, alpha, beta
    dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)

    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, alpha, beta, array_in, array_out)
    
end subroutine mf_wq_get_kcu_3d_csr
