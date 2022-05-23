! ==========================
! module :: assembly for IGA (using sum-factorization)
! author :: Joaquin Cornejo
! ==========================

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine iga_get_capacity_3d(nb_cols_total, capacity_coefs, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, W_u, data_B_v, W_v, data_B_w, W_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, & 
                                data_result, indi_result, indj_result)

    !! Computes a matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, W_u, data_B_v, W_v, data_B_w, W_w
    dimension ::    data_B_u(size_data_u), W_u(nb_cols_u), &
                    data_B_v(size_data_v), W_v(nb_cols_v), &
                    data_B_w(size_data_w), W_w(nb_cols_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer :: i, size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    double precision, dimension(:), allocatable :: data_W00_u, data_W00_v, data_W00_w

    integer :: nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, dummy1, dummy2, dummy3
    integer, allocatable, dimension(:) :: indi_I_temp, indj_I_temp

    ! Get indexes of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, data_B_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    allocate(indi_I_temp(nb_rows_u*nb_rows_v+1), indj_I_temp(size_data_I_u*size_data_I_v))
    call get_indexes_kron_product(nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
            nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
            nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp)

    call get_indexes_kron_product(nb_rows_w, nb_rows_w, size_data_I_w, indi_I_w, indj_I_w, &
            nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp, &
            dummy1, dummy2, dummy3, indi_result, indj_result)
    deallocate(indi_I_temp, indj_I_temp)

    ! Initialize 
    data_result = 0.d0
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w

    allocate(data_W00_u(size_data_u))
    do i = 1, size_data_u
        data_W00_u(i) = data_B_u(i) * W_u(indj_u(i))
    end do

    allocate(data_W00_v(size_data_v))
    do i = 1, size_data_v
        data_W00_v(i) = data_B_v(i) * W_v(indj_v(i))
    end do

    allocate(data_W00_w(size_data_w))
    do i = 1, size_data_w
        data_W00_w(i) = data_B_w(i) * W_w(indj_w(i))
    end do

    ! Get values
    call csr_get_matrix_3d(capacity_coefs, nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_W00_u, data_B_v, data_W00_v, data_B_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)
    deallocate(data_W00_u, data_W00_v, data_W00_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine iga_get_capacity_3d

subroutine iga_get_conductivity_3d(nb_cols_total, cond_coefs, &
                                    nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, &
                                    nb_rows_w, nb_cols_w, &
                                    size_data_u, size_data_v, size_data_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B0_u, data_B1_u, W_u, &
                                    data_B0_v, data_B1_v, W_v, &
                                    data_B0_w, data_B1_w, W_w, &
                                    size_data_I_u, size_data_I_v, size_data_I_w, & 
                                    data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 3D case
    !! IN CSR FORMAT
                
    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    double precision, intent(out) :: data_result
    dimension :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer :: i, size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    double precision, dimension(:), allocatable :: data_W00_u, data_W00_v, data_W00_w
    double precision, dimension(:), allocatable :: data_W11_u, data_W11_v, data_W11_w

    integer :: nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, dummy1, dummy2, dummy3
    integer, allocatable, dimension(:) :: indi_I_temp, indj_I_temp

    ! Get indexes of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    allocate(indi_I_temp(nb_rows_u*nb_rows_v+1), indj_I_temp(size_data_I_u*size_data_I_v))
    call get_indexes_kron_product(nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
            nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
            nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp)

    call get_indexes_kron_product(nb_rows_w, nb_rows_w, size_data_I_w, indi_I_w, indj_I_w, &
            nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp, &
            dummy1, dummy2, dummy3, indi_result, indj_result)
    deallocate(indi_I_temp, indj_I_temp)

    ! Initialize 
    data_result = 0.d0
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w

    allocate(data_W00_u(size_data_u), data_W11_u(size_data_u))
    do i = 1, size_data_u
        data_W00_u(i) = data_B0_u(i) * W_u(indj_u(i))
        data_W11_u(i) = data_B1_u(i) * W_u(indj_u(i))
    end do

    allocate(data_W00_v(size_data_v), data_W11_v(size_data_v))
    do i = 1, size_data_v
        data_W00_v(i) = data_B0_v(i) * W_v(indj_v(i))
        data_W11_v(i) = data_B1_v(i) * W_v(indj_v(i))
    end do

    allocate(data_W00_w(size_data_w), data_W11_w(size_data_w))
    do i = 1, size_data_w
        data_W00_w(i) = data_B0_w(i) * W_w(indj_w(i))
        data_W11_w(i) = data_B1_w(i) * W_w(indj_w(i))
    end do

    ! Get values
    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    ! Get W = W00_w x W00_v x W11_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B1_u, data_W11_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W10_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B1_u, data_W00_u, data_B0_v, data_W11_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_w x W00_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B1_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W11_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! Get W = W00_w x W01_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W11_u, data_B1_v, data_W00_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W00_u, data_B1_v, data_W11_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_w x W01_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W00_u, data_B1_v, data_W00_v, data_B0_w, data_W11_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! Get W = W01_w x W00_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W11_u, data_B0_v, data_W00_v, data_B1_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W01_w x W10_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W00_u, data_B0_v, data_W11_v, data_B1_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B1_w, data_W11_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)
    deallocate(data_W00_u, data_W00_v, data_W00_w)
    deallocate(data_W11_u, data_W11_v, data_W11_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine iga_get_conductivity_3d

subroutine iga_get_source_3D(nb_cols_total, source_coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, W_u, data_B_v, W_v, data_B_w, W_w, &
                            source_vector)
    !! Computes source vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in):: nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(size_data_u), data_B_v(size_data_v), data_B_w(size_data_w)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Local data
    ! ------------------
    double precision, dimension(:), allocatable :: data_W00_u, data_W00_v, data_W00_w
    integer :: i

    ! Calculate equivalent weight
    allocate(data_W00_u(size_data_u))
    do i = 1, size_data_u
        data_W00_u(i) = data_B_u(i) * W_u(indj_u(i))
    end do

    allocate(data_W00_v(size_data_v))
    do i = 1, size_data_v
        data_W00_v(i) = data_B_v(i) * W_v(indj_v(i))
    end do

    allocate(data_W00_w(size_data_w))
    do i = 1, size_data_w
        data_W00_w(i) = data_B_w(i) * W_w(indj_w(i))
    end do

    ! Find vector 
    source_vector = 0.d0
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, indi_u, indj_u, data_W00_u, &
                                size_data_v, indi_v, indj_v, data_W00_v, &
                                size_data_w, indi_w, indj_w, data_W00_w, &
                                source_coefs, source_vector)

end subroutine iga_get_source_3D