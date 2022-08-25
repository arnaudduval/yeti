! ====================================================
! module :: IGA - Matrix free methods 
! author :: Joaquin Cornejo
! ====================================================

subroutine iga_find_conductivity_diagonal_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, Kdiag)
    !! Find the diagonal of conductivity matrix
    !! Indices in CSR format
    
    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total,nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    
    double precision, intent(out) :: Kdiag
    dimension :: Kdiag(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------------------
    ! Find diagonal
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)
    double precision, allocatable, dimension(:, :) :: data_W_u, data_W_v, data_W_w
    double precision :: Kdiag_temp
    dimension :: Kdiag_temp(nr_u*nr_v*nr_w)

    allocate(data_W_u(nnz_u, 2))
    do i = 1, nnz_u
        data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
    end do
    
    allocate(data_W_v(nnz_v, 2))
    do i = 1, nnz_v
        data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
    end do

    allocate(data_W_w(nnz_w, 2))
    do i = 1, nnz_w
        data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
    end do

    ! Initialize
    Kdiag = 0.d0
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            call find_physical_diag_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                                data_W_u(:, alpha(1)), data_W_v(:, alpha(2)), data_W_w(:, alpha(3)), Kdiag_temp)
            Kdiag = Kdiag + Kdiag_temp
        end do
    end do

end subroutine iga_find_conductivity_diagonal_3d

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine mf_iga_get_cu_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indices must be in CSR format
    
    use tensor_methods
    use omp_lib
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)

    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total)

    ! Local data 
    ! ------------------
    integer :: ju, jv, jw, genPos, nb_tasks
    double precision :: coefs_new
    dimension :: coefs_new(nc_total)
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t

    ! Get new coefficients 
    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_total/nb_tasks) 
    ! Initialize coefficients
    do jw = 1, nc_w
        do jv = 1, nc_v
            do ju = 1, nc_u
                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                coefs_new(genPos) = coefs(genPos)*W_u(ju)*W_v(jv)*W_w(jw)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! Eval B.transpose * array_in
    allocate(array_temp_1(nc_total))
    call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), & 
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            array_input, array_temp_1)

    ! Evaluate diag(coefs) * array_temp1
    allocate(array_temp_1t(nc_total))
    array_temp_1t = array_temp_1*coefs_new
    deallocate(array_temp_1)

    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_B_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_B_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_B_w(:, 1), &
                            array_temp_1t, array_output)
    deallocate(array_temp_1t)

end subroutine mf_iga_get_cu_3d

subroutine mf_iga_get_ku_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, array_input, array_output)
    !! Computes K.u in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    use omp_lib
    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    ! Csr format
    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    
    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total)

    ! Local data 
    ! ------------------
    integer :: ju, jv, jw, genPos, nb_tasks
    double precision :: coefs_new
    dimension :: coefs_new(d, d, nc_total)
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t, array_temp_1tt
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)

    ! Get new coefficients 
    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_total/nb_tasks) 
    ! Initialize coefficients
    do jw = 1, nc_w
        do jv = 1, nc_v
            do ju = 1, nc_u
                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                coefs_new(:, :, genPos) = coefs(:, :, genPos)*W_u(ju)*W_v(jv)*W_w(jw)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

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
            array_temp_1t = array_temp_1*coefs_new(i, j, :)
            call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_B_u(:, alpha(1)), &
                                    nnz_v, indi_v, indj_v, data_B_v(:, alpha(2)), &
                                    nnz_w, indi_w, indj_w, data_B_w(:, alpha(3)), & 
                                    array_temp_1t, array_temp_1tt)
            array_output = array_output + array_temp_1tt
        end do
    end do
    deallocate(array_temp_1, array_temp_1t, array_temp_1tt)
    
end subroutine mf_iga_get_ku_3d

subroutine mf_iga_get_ku_3d_csr(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, array_input, array_output)
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
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)
    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total)

    ! Local data
    ! ------------------
    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_iga_get_ku_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, array_input, array_output)
    
end subroutine mf_iga_get_ku_3d_csr

! ----------------------------------------
! Conjugate gradient
! ----------------------------------------

subroutine mf_iga_steady_heat_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIter, epsilon, & 
                        Method, nnz_cond, cond, JJ, directsol, x, RelRes, RelError)
    !! Conjugate gradient with or without preconditioner 
    !! CSR FORMAT
                        
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
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2) 
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    character(len = 10) :: Method
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, b, JJ, directsol
    dimension :: b(nr_total), JJ(3, 3, nc_u*nc_v*nc_w), directsol(nr_total)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_total), RelRes(nbIter+1), RelError(nbIter+1)

    ! Local data
    ! ------------------
    ! Pre / Conjugate gradient algoritm
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: rsold, rsnew, alpha
    double precision :: r, p, Ap, dummy, z
    dimension :: r(nr_total), p(nr_total), Ap(nr_total), dummy(nr_total), z(nr_total)
    integer :: iter, i, t_robin(2)

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
    double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w

    ! Initialize B transpose in CSR format
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
            r = b; p = r
            rsold = dot_product(r, r)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do iter = 1, nbIter
                ! Calculate Ann xn 
                call mf_iga_get_Ku_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                    data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, p, Ap)
                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                ! Set relative value of residual 
                RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))

                if (RelRes(iter+1).le.epsilon) exit
            
                rsnew = dot_product(r, r)
                p = r + rsnew/rsold * p
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
        allocate(data_W_u(nnz_u, 2), Kdiag_u(nr_u), Mdiag_u(nr_u))
        do i = 1, nnz_u
            data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
        end do
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 2), t_robin, D_u, U_u, Kdiag_u, Mdiag_u)
        deallocate(data_W_u)
        
        allocate(data_W_v(nnz_v, 2), Kdiag_v(nr_v), Mdiag_v(nr_v))
        do i = 1, nnz_v
            data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
        end do
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 2), t_robin, D_v, U_v, Kdiag_v, Mdiag_v)    
        deallocate(data_W_v)

        allocate(data_W_w(nnz_w, 2), Kdiag_w(nr_w), Mdiag_w(nr_w))
        do i = 1, nnz_w
            data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
        end do
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 2), t_robin, D_w, U_w, Kdiag_w, Mdiag_w)  
        deallocate(data_W_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        ! Find diagonal of eigen values
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(nr_total))
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
            call iga_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, &
                                nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, Dphysical)
        end if
        
        if (nbIter.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            r = b; dummy = r 
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
            end if
            call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                call fd_sqr_scaling(nr_total, Dparametric, Dphysical, z) 
            end if
            p = z
            rsold = dot_product(r, z)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do iter = 1, nbIter
                ! Calculate Ann xn 
                call mf_iga_get_Ku_3D(coefs, nr_total, nc_total, nr_u, nc_u, &
                            nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, p, Ap)

                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                ! Set relative value of residual 
                RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))

                if (RelRes(iter+1).le.epsilon) exit       

                dummy = r
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, z) 
                end if
                rsnew = dot_product(r, z)
                                
                p = z + rsnew/rsold * p
                rsold = rsnew
            end do
        end if
    end if

end subroutine mf_iga_steady_heat_3d

subroutine mf_iga_interpolate_cp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIter, epsilon, x, RelRes)
    !! Conjugate gradient with or without preconditioner 
    !! CSR FORMAT
                        
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    integer, intent(in) :: nbIter
    double precision, intent(in) :: epsilon, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, RelRes
    dimension :: x(nr_total), RelRes(nbIter+1)

    ! Local data
    ! ------------------
    ! Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha
    double precision :: r, p, Ap, dummy, z
    dimension :: r(nr_total), p(nr_total), Ap(nr_total), dummy(nr_total), z(nr_total)
    integer :: iter, i, dorobin(2)

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w
    double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
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
    
    allocate(data_W_u(nnz_u, 2), Kdiag_u(nr_u), Mdiag_u(nr_u))
    do i = 1, nnz_u
        data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
    end do
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 2), dorobin, D_u, U_u, Kdiag_u, Mdiag_u)
    deallocate(data_W_u)
    
    allocate(data_W_v(nnz_v, 2), Kdiag_v(nr_v), Mdiag_v(nr_v))
    do i = 1, nnz_v
        data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
    end do
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 2), dorobin, D_v, U_v, Kdiag_v, Mdiag_v)    
    deallocate(data_W_v)

    allocate(data_W_w(nnz_w, 2), Kdiag_w(nr_w), Mdiag_w(nr_w))
    do i = 1, nnz_w
        data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
    end do
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 2), dorobin, D_w, U_w, Kdiag_w, Mdiag_w)  
    deallocate(data_W_w)

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    x = 0.d0; RelRes = 0.d0
    r = b; dummy = r 
    call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, dummy, z)
    p = z
    rsold = dot_product(r, z)
    RelRes(1) = 1.d0

    do iter = 1, nbIter
        ! Calculate Ann xn 
        call mf_iga_get_Cu_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_B_u, data_B_v, data_B_w, p, Ap)

        alpha = rsold/dot_product(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap

        ! Set relative value of residual 
        RelRes(iter+1) = maxval(abs(r))/maxval(abs(b))

        if (RelRes(iter+1).le.epsilon) exit
        
        dummy = r
        call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, dummy, z)
        rsnew = dot_product(r, z)
                        
        p = z + rsnew/rsold * p
        rsold = rsnew
    end do

end subroutine mf_iga_interpolate_cp_3d
