! ==========================
! module :: Assembly for IGA-WQ 
! author :: Joaquin Cornejo
! ==========================

subroutine eval_thermal_coefficient(dime, nnz, JJ, nnz_K, KK, nnz_C, CC, Kcoef, Ccoef, info)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------    
    integer, intent(in) :: dime, nnz, nnz_K, nnz_C
    double precision, intent(in) :: JJ, KK, CC
    dimension :: JJ(dime, dime, nnz), KK(dime, dime, nnz_K), CC(nnz_C)

    integer, intent(out) :: info
    double precision, intent(out) :: Kcoef, Ccoef
    dimension :: Kcoef(dime, dime, nnz), Ccoef(nnz)

    ! Local data
    ! -----------
    integer :: i, nb_tasks
    double precision :: Jt, detJt, invJt
    dimension :: Jt(dime, dime), invJt(dime, dime)
    double precision :: Kt
    dimension :: Kt(dime, dime)  

    info = 1

    ! Verifiy nnz_KK value
    if (nnz_K.eq.1) then 
        !$OMP PARALLEL PRIVATE(Jt, invJt, detJt, Kt)
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
        do i = 1, nnz
            ! Get individual jacobien
            Jt = JJ(:, :, i)

            ! Evaluate inverse
            call MatrixInv(invJt, Jt, detJt, dime)

            ! For K = invJ * prop * detJ * transpose(invJ) 
            Kt = detJt * matmul(invJt, KK(:, :, 1)) 
            Kcoef(:, :, i) = matmul(Kt, transpose(invJt))
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 
    else if (nnz_K.eq.nnz) then 
        !$OMP PARALLEL PRIVATE(Jt, invJt, detJt, Kt)
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
        do i = 1, nnz
            ! Get individual jacobien
            Jt = JJ(:, :, i)

            ! Evaluate inverse
            call MatrixInv(invJt, Jt, detJt, dime)

            ! For K = invJ * prop * detJ * transpose(invJ) 
            Kt = detJt * matmul(invJt, KK(:, :, i)) 
            Kcoef(:, :, i) = matmul(Kt, transpose(invJt))

        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    else 
        info = 0
        print*, "Error computing thermal coefficient (Conductivity)"
    end if

    if (nnz_C.eq.1) then 
        !$OMP PARALLEL PRIVATE(Jt, detJt)
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
        do i = 1, nnz
            ! Get individual jacobien
            Jt = JJ(:, :, i)

            ! For C = detJ  * prop
            call MatrixDet(Jt, detJt, dime)
            Ccoef(i) = detJt * CC(1)
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    else if (nnz_C.eq.nnz) then
        !$OMP PARALLEL PRIVATE(Jt, detJt)
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
        do i = 1, nnz
            ! Get individual jacobien
            Jt = JJ(:, :, i)

            ! For C = detJ  * prop
            call MatrixDet(Jt, detJt, dime)
            Ccoef(i) = detJt * CC(i)
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    else
        info = 0
        print*, "Error computing thermal coefficient (Capacity) "
    end if

end subroutine eval_thermal_coefficient

subroutine jacobien_physicalposition_3d(nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        nb_rows_w, nb_cols_w, &
                                        size_data_u, size_data_v, size_data_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, &
                                        ctrlpts, jacob, physical_pos, detJ)
    !! Computes jacobien in 3D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) ::  nb_rows_u, nb_rows_v, nb_rows_w, &
                            nb_cols_u, nb_cols_v, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(size_data_u, 2), data_B_v(size_data_v, 2), data_B_w(size_data_w, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(3, nb_rows_u*nb_rows_v*nb_rows_w)

    double precision, intent(out) :: jacob
    dimension :: jacob(3, 3, nb_cols_u*nb_cols_v*nb_cols_w)
    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(3, nb_cols_u*nb_cols_v*nb_cols_w)
    double precision, intent(out) :: detJ
    dimension :: detJ(nb_cols_u*nb_cols_v*nb_cols_w)

    ! Local data
    !-----------------
    integer :: i, j, k, nb_tasks, alpha
    dimension :: alpha(3)
    double precision :: result_temp, detJt
    dimension ::  result_temp(nb_cols_u*nb_cols_v*nb_cols_w)

    ! Csr format (Transpose)
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nb_cols_u+1), indi_T_v(nb_cols_v+1), indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(size_data_u), indj_T_v(size_data_v), indj_T_w(size_data_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(size_data_u, 2), data_BT_v(size_data_v, 2), data_BT_w(size_data_w, 2)

    ! Initialize
    call csr2csc(2, nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nb_rows_w, nb_cols_w, size_data_w, data_B_w, indj_w, indi_w, &
                data_BT_w, indj_T_w, indi_T_w)

    ! Compute physical position
    do i = 1, 3
        call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                            nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                            size_data_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            size_data_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            size_data_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            ctrlpts(i, :), result_temp)
        physical_pos(i, :) = result_temp
    end do

    ! Compute jacobien matrix
    do j = 1, 3
        do i = 1, 3
            alpha = 1; alpha(j) = 2
            call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                            nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                            size_data_u, indi_T_u, indj_T_u, data_BT_u(:, alpha(1)), &
                            size_data_v, indi_T_v, indj_T_v, data_BT_v(:, alpha(2)), &
                            size_data_w, indi_T_w, indj_T_w, data_BT_w(:, alpha(3)), &
                            ctrlpts(i, :), result_temp)
            jacob(i, j, :) = result_temp
        end do
    end do

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, k), detJt, 3)

        ! Assign values
        detJ(k) = detJt
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_3d

subroutine interpolation_3d(nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, u_ctrlpts, u_interp)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in) ::  nb_rows_u, nb_rows_v, nb_rows_w, &
                            nb_cols_u, nb_cols_v, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(size_data_u), data_B_v(size_data_v), data_B_w(size_data_w)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nb_rows_u*nb_rows_v*nb_rows_w)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nb_cols_u*nb_cols_v*nb_cols_w)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nb_cols_u+1), indi_T_v(nb_cols_v+1), indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(size_data_u), indj_T_v(size_data_v), indj_T_w(size_data_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(size_data_u), data_BT_v(size_data_v), data_BT_w(size_data_w)

    ! Initialize
    call csr2csc(1, nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(1, nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(1, nb_rows_w, nb_cols_w, size_data_w, data_B_w, indj_w, indi_w, &
                data_BT_w, indj_T_w, indi_T_w)

    ! Interpolation
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_BT_u, &
                                size_data_v, indi_T_v, indj_T_v, data_BT_v, &
                                size_data_w, indi_T_w, indj_T_w, data_BT_w, &
                                u_ctrlpts, u_interp)

end subroutine interpolation_3d

subroutine jacobien_physicalposition_2d(nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        size_data_u, size_data_v, &
                                        indi_u, indj_u, indi_v, indj_v, &
                                        data_B_u, data_B_v, &
                                        ctrlpts, jacob, physical_pos, detJ)
    !! Computes jacobien in 2D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) :: nb_rows_u, nb_rows_v, nb_cols_u, nb_cols_v, size_data_u, size_data_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(size_data_u, 2), data_B_v(size_data_v, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(2, nb_rows_u*nb_rows_v)

    double precision, intent(out) :: jacob
    dimension :: jacob(2, 2, nb_cols_u*nb_cols_v)
    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(2, nb_cols_u*nb_cols_v)
    double precision, intent(out) :: detJ
    dimension :: detJ(nb_cols_u*nb_cols_v)

    ! Local data
    !-----------------
    integer :: i, j, k, nb_tasks, alpha
    dimension :: alpha(2)
    double precision :: result_temp
    dimension ::  result_temp(nb_cols_u*nb_cols_v)
    double precision :: detJt

    ! Csr format (Transpose)
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nb_cols_u+1), indi_T_v(nb_cols_v+1)
    integer ::  indj_T_u, indj_T_v
    dimension :: indj_T_u(size_data_u), indj_T_v(size_data_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(size_data_u, 2), data_BT_v(size_data_v, 2)

    ! Initialize
    call csr2csc(2, nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)
    
    ! Compute physical position
    do i = 1, 3
        call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, nb_cols_v, nb_rows_v, &
                            size_data_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            size_data_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            ctrlpts(i, :), result_temp)
        physical_pos(i, :) = result_temp
    end do

    ! Compute jacobien matrix
    do j = 1, 2
        do i = 1, 2
            alpha = 1; alpha(j) = 2
            call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, nb_cols_v, nb_rows_v, &
                            size_data_u, indi_T_u, indj_T_u, data_BT_u(:, alpha(1)), &
                            size_data_v, indi_T_v, indj_T_v, data_BT_v(:, alpha(2)), &
                            ctrlpts(i, :), result_temp)
            jacob(i, j, :) = result_temp
        end do
    end do

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, k), detJt, 2)

        ! Assign values
        detJ(k) = detJt
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_2d

subroutine interpolation_2d(nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in) ::  nb_rows_u, nb_rows_v, &
                            nb_cols_u, nb_cols_v, &
                            size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(size_data_u), data_B_v(size_data_v)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nb_rows_u*nb_rows_v)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nb_cols_u*nb_cols_v)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nb_cols_u+1), indi_T_v(nb_cols_v+1)
    integer :: indj_T_u, indj_T_v
    dimension :: indj_T_u(size_data_u), indj_T_v(size_data_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(size_data_u), data_BT_v(size_data_v)

    ! Initialize
    call csr2csc(1, nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(1, nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)

    ! Interpolation
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_BT_u, &
                                size_data_v, indi_T_v, indj_T_v, data_BT_v, &
                                u_ctrlpts, u_interp)

end subroutine interpolation_2d