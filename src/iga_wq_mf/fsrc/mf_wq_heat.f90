! ====================================================
! module :: Weighted quadrature - Matrix free methods 
! author :: Joaquin Cornejo
! ====================================================

subroutine wq_find_conductivity_diagonal_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Kdiag)
    
    use tensor_methods
    implicit none
    ! Input / output 
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

    double precision, intent(out) :: Kdiag
    dimension :: Kdiag(nr_u*nr_v*nr_w)

    ! Local data
    ! ---------------
    integer :: i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)
    double precision :: Kdiag_temp
    dimension :: Kdiag_temp(nr_u*nr_v*nr_w)

    ! Initialize
    Kdiag = 0.d0 
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            zeta = beta + (alpha - 1)*2
            call find_physical_diag_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                                    data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), data_W_w(:, zeta(3)), Kdiag_temp)
            Kdiag = Kdiag + Kdiag_temp        
        end do
    end do

end subroutine wq_find_conductivity_diagonal_3d

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------
subroutine mf_wq_get_cu_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
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

    double precision, intent(in) :: array_input
    dimension :: array_input(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------------- 
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t

    ! Initialize
    allocate(array_temp_1(nc_total))

    ! Eval B.transpose * array_in
    call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                                array_input, array_temp_1)

    ! Evaluate diag(coefs) * array_temp1
    allocate(array_temp_1t(nc_total))
    array_temp_1t = array_temp_1 * coefs
    deallocate(array_temp_1)

    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                array_temp_1t, array_output)
    deallocate(array_temp_1t)
    
end subroutine mf_wq_get_cu_3d

subroutine mf_wq_get_cu_3d_csr(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
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

    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total)

    ! Local data 
    ! ----------------- 
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_wq_get_cu_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_input, array_output)

end subroutine mf_wq_get_cu_3d_csr

subroutine mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_input, array_output)
    !! Computes K.u in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
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

    double precision, intent(in) :: array_input
    dimension :: array_input(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_u*nr_v*nr_w)

    ! Local data 
    ! ------------------
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t, array_temp_1tt
    integer :: i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)

    ! Initialize
    array_output = 0.d0
    allocate(array_temp_1(nc_total), array_temp_1t(nc_total), array_temp_1tt(nc_total))
    do j = 1, d
        beta = 1; beta(j) = 2
        call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                    nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                    nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                                    array_input, array_temp_1)
        do i = 1, d
            alpha = 1; alpha(i) = 2
            zeta = beta + (alpha - 1)*2
            array_temp_1t = array_temp_1 * coefs(i, j, :)
            call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                        array_temp_1t, array_temp_1tt)
            array_output = array_output + array_temp_1tt
        end do
    end do
    deallocate(array_temp_1, array_temp_1t, array_temp_1tt)
    
end subroutine mf_wq_get_ku_3d

subroutine mf_wq_get_ku_3d_csr(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_input, array_output)
    !! Computes K.u in 3D case
    !! Indices must be in CSR format

    implicit none 
    ! Input / output data
    ! ---------------------
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
    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total)

    ! Local data
    ! ------------------
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w,&
                        data_W_u, data_W_v, data_W_w, array_input, array_output)
    
end subroutine mf_wq_get_ku_3d_csr

! ----------------------------------------
! Conjugate gradient for steady problems
! ----------------------------------------

subroutine mf_wq_steady_heat_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIter, epsilon, Method, nnz_cond, cond, JJ, directsol, x, RelRes, RelError)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
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

    character(len=10), intent(in) :: Method
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, b, JJ, directsol
    dimension :: b(nr_total), JJ(3, 3, nc_total), directsol(nr_total)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_total), RelRes(nbIter+1), RelError(nbIter+1)

    ! Local data
    ! ------------------
    ! Pre / Conjugate gradient algoritm
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, Ap, s, As, dummy, ptilde, Aptilde, stilde, Astilde
    dimension ::    r(nr_total), rhat(nr_total), p(nr_total), & 
                    Ap(nr_total), As(nr_total), s(nr_total), dummy(nr_total), &
                    ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
    integer :: iter, t_robin(2)

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Initiate variables
    x = 0.d0; RelRes = 0.d0; RelError = 0.d0

    if (Method.eq.'WP') then 
        if (nbIter.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            r = b; rhat = r; p = r
            rsold = dot_product(r, rhat)

            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do iter = 1, nbIter
                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, Ap)
                alpha = rsold/dot_product(Ap, rhat)
                s = r - alpha*Ap

                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, s, As)
                omega = dot_product(As, s)/dot_product(As, As)
                x = x + alpha*p + omega*s
                r = s - omega*As

                RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                
                if (RelRes(iter+1).le.epsilon) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew
            end do
        end if

    else  
        ! Initialize
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0

        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            ! --------------------------------------------
            ! DIAGONAL DECOMPOSITION
            ! --------------------------------------------            
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((Method.eq.'JMS').or.(Method.eq.'JMC')) then 
            ! --------------------------------------------
            ! MY METHOD
            ! --------------------------------------------
            ! Find dimensions and conductivity
            call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, Lu, Lv, Lw)
            call conductivity_mean_3d(nc_u, nc_v, nc_w, nnz_cond, cond, lamb_u, lamb_v, lamb_w)
            
            Kcoef_u = lamb_u*Lv*Lw/Lu
            Kcoef_v = lamb_v*Lw*Lu/Lv
            Kcoef_w = lamb_w*Lu*Lv/Lw

        end if

        ! --------------------------------------------
        ! EIGEN DECOMPOSITION
        ! -------------------------------------------- 
        t_robin = (/0, 0/)
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))

        allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), t_robin, D_u, U_u, Kdiag_u, Mdiag_u)

        allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), t_robin, D_v, U_v, Kdiag_v, Mdiag_v)    

        allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), t_robin, D_w, U_w, Kdiag_w, Mdiag_w)   
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        ! Find diagonal of eigen values
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
        allocate(Deigen(nr_total))
        I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        end if

        if (nbIter.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            r = b; rhat = r; p = r
            rsold = dot_product(r, rhat)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do iter = 1, nbIter

                dummy = p
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if 
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, ptilde)  
                end if

                call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)

                alpha = rsold/dot_product(Aptilde, rhat)
                s = r - alpha*Aptilde

                dummy = s
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy)  
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, stilde) 
                end if

                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)

                omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
                x = x + alpha*ptilde + omega*stilde
                r = s - omega*Astilde    
                
                RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                
                if (RelRes(iter+1).le.epsilon) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                rsold = rsnew
            end do
        end if
    end if

end subroutine mf_wq_steady_heat_3d

subroutine mf_wq_interpolate_cp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterations, epsilon, x, RelRes)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
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

    integer, intent(in) :: nbIterations
    double precision, intent(in) :: epsilon, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, RelRes
    dimension :: x(nr_total), RelRes(nbIterations+1)

    ! Local data
    ! ------------------
    ! Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, s, dummy, ptilde, Aptilde, stilde, Astilde
    dimension ::    r(nr_total), rhat(nr_total), p(nr_total), & 
                    s(nr_total), dummy(nr_total), ptilde(nr_total), &
                    Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
    integer :: iter, dorobin(2)

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Initiate variables
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    
    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    dorobin = (/0, 0/)
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))

    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), dorobin, D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), dorobin, D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), dorobin, D_w, U_w, Kdiag_w, Mdiag_w)  

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    x = 0.d0; RelRes = 0.d0
    r = b; rhat = r; p = r
    rsold = dot_product(r, rhat)
    RelRes(1) = 1.d0

    do iter = 1, nbIterations
        dummy = p
        call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, dummy, ptilde)

        call mf_wq_get_cu_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, ptilde, Aptilde)

        alpha = rsold/dot_product(Aptilde, rhat)
        s = r - alpha*Aptilde

        dummy = s
        call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, dummy, stilde)

        call mf_wq_get_cu_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, stilde, Astilde)

        omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
        x = x + alpha*ptilde + omega*stilde
        r = s - omega*Astilde    
        
        RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))
        
        if (RelRes(iter+1).le.epsilon) exit

        rsnew = dot_product(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        rsold = rsnew
    end do

end subroutine mf_wq_interpolate_cp_3d
