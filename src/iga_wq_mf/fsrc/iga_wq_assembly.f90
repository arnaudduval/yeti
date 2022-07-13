! ==========================
! module :: assembly for IGA-WQ 
! author :: Joaquin Cornejo
! modules :: operateurs.f90 (MatrixInv and MatrixDet)
! ==========================

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine wq_get_capacity_3d(nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, &
                            nnz_I_u, nnz_I_v, nnz_I_w, & 
                            data_result, indi_result, indj_result)
    !! Computes a matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(size_data_u, 2), data_W_u(size_data_u, 4), &
                    data_B_v(size_data_v, 2), data_W_v(size_data_v, 4), &
                    data_B_w(size_data_w, 2), data_W_w(size_data_w, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ---------------
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: nnz_result, dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nb_rows_w, nb_rows_w, nnz_I_w, indi_I_w, indj_I_w, &
                                nb_rows_v, nb_rows_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data 
    nnz_result = nnz_I_u*nnz_I_v*nnz_I_w
    call csr_get_matrix_3d(capacity_coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u(:, 1), data_B_v(:, 1), data_B_w(:, 1), &
                            data_W_u(:, 1), data_W_v(:, 1), data_W_w(:, 1), &
                            nnz_I_u, nnz_I_v, nnz_I_w, & 
                            indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, nnz_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_3d

subroutine wq_get_conductivity_3d(nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, & 
                                nnz_I_u, nnz_I_v, nnz_I_w, & 
                                data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(size_data_u, 2), data_W_u(size_data_u, 4), &
                    data_B_v(size_data_v, 2), data_W_v(size_data_v, 4), &
                    data_B_w(size_data_w, 2), data_W_w(size_data_w, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w 

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(nnz_I_u*nnz_I_v*nnz_I_w)
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: nnz_result, dummy1, dummy2, dummy3, i, j, alpha, beta, zeta
    dimension :: alpha(3), beta(3), zeta(3)

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nb_rows_w, nb_rows_w, nnz_I_w, indi_I_w, indj_I_w, &
                                nb_rows_v, nb_rows_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v*nnz_I_w
    data_result = 0.d0
    do j = 1, 3
        do i = 1, 3
            beta = 1; beta(j) = 2
            alpha = 1; alpha(i) = 2
            zeta = alpha + (beta - 1)*2
            call csr_get_matrix_3d(cond_coefs(i, j, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                        data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), data_W_w(:, zeta(3)), &
                        nnz_I_u, nnz_I_v, nnz_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, nnz_result, data_result_temp)
            data_result = data_result + data_result_temp
        end do
    end do
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_conductivity_3d

subroutine wq_get_source_3d(nb_cols_total, source_coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, &
                            source_vector)
    !! Computes source vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(size_data_u, 4), data_W_v(size_data_v, 4), data_W_w(size_data_w, 4)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Find vector
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, indi_u, indj_u, data_W_u(:, 1), &
                        size_data_v, indi_v, indj_v, data_W_v(:, 1), &
                        size_data_w, indi_w, indj_w, data_W_w(:, 1), &
                        source_coefs, source_vector)

end subroutine wq_get_source_3d

! ----------------------------------------
! Assembly in 2D
! ----------------------------------------

subroutine wq_get_capacity_2d(nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_W_u, data_B_v, data_W_v, &
                            nnz_I_u, nnz_I_v, & 
                            data_result, indi_result, indj_result)
    !! Computes a matrix in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(size_data_u, 2), data_W_u(size_data_u, 4), &
                    data_B_v(size_data_v, 2), data_W_v(size_data_v, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v+1), &
                    indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ---------------
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: nnz_result, dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nb_rows_v, nb_rows_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v
    call csr_get_matrix_2d(capacity_coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u(:, 1), data_B_v(:, 1), data_W_u(:, 1), data_W_v(:, 1), &
                            nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                            indi_result, nnz_result, data_result)
    deallocate(indj_I_u, indj_I_v)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_2d

subroutine wq_get_conductivity_2d(nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                                size_data_u, size_data_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_W_u, data_B_v, data_W_v, &
                                nnz_I_u, nnz_I_v, & 
                                data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(2, 2, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(size_data_u, 2), data_W_u(size_data_u, 4), &
                    data_B_v(size_data_v, 2), data_W_v(size_data_v, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v 

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v+1), &
                    indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(nnz_I_u*nnz_I_v)
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: nnz_result, dummy1, dummy2, dummy3, i, j, alpha, beta, zeta
    dimension :: alpha(2), beta(2), zeta(2)

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nb_rows_v, nb_rows_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v
    data_result = 0.d0
    do j = 1, 3
        do i = 1, 3
            beta = 1; beta(j) = 2
            alpha = 1; alpha(i) = 2
            zeta = alpha + (beta - 1)*2
            call csr_get_matrix_2d(cond_coefs(i, j, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v,&
                        data_B_u(:, beta(1)), data_B_v(:, beta(2)), &
                        data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), &
                        nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, nnz_result, data_result_temp)
            data_result = data_result + data_result_temp
        end do
    end do
    
    deallocate(indj_I_u, indj_I_v)   
    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_conductivity_2d

subroutine wq_get_source_2d(nb_cols_total, source_coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, source_vector)
    !! Computes source vector in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(size_data_u, 4), data_W_v(size_data_v, 4)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_rows_u*nb_rows_v)

    ! Find vector
    call tensor2d_dot_vector_sp(nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                        size_data_u, indi_u, indj_u, data_W_u(:, 1), &
                        size_data_v, indi_v, indj_v, data_W_v(:, 1), &
                        source_coefs, source_vector)

end subroutine wq_get_source_2d