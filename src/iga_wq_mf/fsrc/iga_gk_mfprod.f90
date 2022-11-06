! ====================================================
! module :: Matrix free methods in IGA approach 
! author :: Joaquin Cornejo
!
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

subroutine iga_find_conductivity_diagonal_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, diag)
    !! Computes the diagonal of conductivity matrix using sum-factorization algorithm
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
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
    
    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 2), data_W_v(nnz_v, 2), data_W_w(nnz_w, 2)
    double precision :: diag_temp
    dimension :: diag_temp(nr_u*nr_v*nr_w)

    ! Compute equivalent weights 
    do i = 1, nnz_u
        data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
    end do
    
    do i = 1, nnz_v
        data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
    end do

    do i = 1, nnz_w
        data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
    end do

    ! Compute diagonal
    diag = 0.d0
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            call csr_get_diag_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                                data_W_u(:, alpha(1)), data_W_v(:, alpha(2)), data_W_w(:, alpha(3)), diag_temp)
            diag = diag + diag_temp
        end do
    end do

end subroutine iga_find_conductivity_diagonal_3d

! ----------------------------------------
! Matrix free multiplication in 3D
! ----------------------------------------

subroutine mf_iga_get_cu_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, array_in, array_out)
    !! Computes C.u where C is the capacity matrix in 3D 
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! -----------
    integer :: ju, jv, jw, genPos, nb_tasks
    double precision :: coefs_temp, array_temp_0, array_temp_1
    dimension :: coefs_temp(nc_total), array_temp_0(nc_total), array_temp_1(nc_total)

    ! Get new coefficients 
    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_total/nb_tasks) 
    do jw = 1, nc_w
        do jv = 1, nc_v
            do ju = 1, nc_u
                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                coefs_temp(genPos) = coefs(genPos)*W_u(ju)*W_v(jv)*W_w(jw)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! Eval B' * array_in
    call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), & 
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            array_in, array_temp_0)

    ! Evaluate diag(coefs) * array_temp_0
    array_temp_1 = array_temp_0*coefs_temp

    ! Eval W * array_temp_1
    call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_B_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_B_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_B_w(:, 1), &
                            array_temp_1, array_out)

end subroutine mf_iga_get_cu_3d

subroutine mf_iga_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! IN CSR FORMAT

    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    integer :: ju, jv, jw, genPos, nb_tasks
    double precision :: coefs_temp, array_temp_0, array_temp_1, array_temp_2
    dimension :: coefs_temp(d,d, nc_total), array_temp_0(nc_total), array_temp_1(nc_total), array_temp_2(nc_total)
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)

    ! Get new coefficients 
    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_total/nb_tasks) 
    do jw = 1, nc_w
        do jv = 1, nc_v
            do ju = 1, nc_u
                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                coefs_temp(:, :, genPos) = coefs(:, :, genPos)*W_u(ju)*W_v(jv)*W_w(jw)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

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
            array_temp_1 = array_temp_0*coefs_temp(i, j, :)
            call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_B_u(:, alpha(1)), &
                                    nnz_v, indi_v, indj_v, data_B_v(:, alpha(2)), &
                                    nnz_w, indi_w, indj_w, data_B_w(:, alpha(3)), & 
                                    array_temp_1, array_temp_2)
            array_out = array_out + array_temp_2
        end do
    end do
    
end subroutine mf_iga_get_ku_3d

subroutine mf_iga_get_ku_3d_csr(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    call mf_iga_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, array_in, array_out)
    
end subroutine mf_iga_get_ku_3d_csr
