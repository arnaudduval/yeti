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
                            data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
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
    double precision, intent(in) :: data_B0_u, data_W00_u, &
                                    data_B0_v, data_W00_v, &
                                    data_B0_w, data_W00_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nb_rows_w, nb_rows_w, size_data_I_w, indi_I_w, indj_I_w, &
                                nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Initialize 
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    call csr_get_matrix_3d(capacity_coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_B0_v, data_B0_w, data_W00_u, data_W00_v, data_W00_w, &
                            size_data_I_u, size_data_I_v, size_data_I_w, & 
                            indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, size_data_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_3d

subroutine wq_get_conductivity_3d(nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_B1_u, &
                                data_W00_u, data_W01_u, &
                                data_W10_u, data_W11_u, &
                                data_B0_v, data_B1_v, &
                                data_W00_v, data_W01_v, &
                                data_W10_v, data_W11_v, &
                                data_B0_w, data_B1_w, &
                                data_W00_w, data_W01_w, &
                                data_W10_w, data_W11_w, &
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
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
                                    data_B0_v, data_B1_v, &
                                    data_W00_v, data_W01_v, data_W10_v, data_W11_v, &
                                    data_B0_w, data_B1_w, &
                                    data_W00_w, data_W01_w, data_W10_w, data_W11_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_W00_u(size_data_u), data_W01_u(size_data_u), &
                    data_W10_u(size_data_u), data_W11_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_W00_v(size_data_v), data_W01_v(size_data_v), &
                    data_W10_v(size_data_v), data_W11_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w), &
                    data_W00_w(size_data_w), data_W01_w(size_data_w), &
                    data_W10_w(size_data_w), data_W11_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w 

    double precision, intent(out) :: data_result
    dimension :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(size_data_I_u*size_data_I_v*size_data_I_w)
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nb_rows_w, nb_rows_w, size_data_I_w, indi_I_w, indj_I_w, &
                                nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Initialize 
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    
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
                        data_B1_u, data_B0_v, data_B0_w, data_W11_u, data_W00_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W10_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B1_u, data_B0_v, data_B0_w, data_W01_u, data_W10_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W10_w x W00_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B1_u, data_B0_v, data_B0_w, data_W01_u, data_W00_v, data_W10_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! Get W = W00_w x W01_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B1_v, data_B0_w, data_W10_u, data_W01_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W11_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W10_w x W01_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W01_v, data_W10_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! Get W = W01_w x W00_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B1_w, data_W10_u, data_W00_v, data_W01_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W01_w x W10_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W10_v, data_W01_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W00_v, data_W11_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

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
                            data_W00_u, data_W00_v, data_W00_w, &
                            source_vector)
    !! Computes source vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W00_v, data_W00_w
    dimension :: data_W00_u(size_data_u), data_W00_v(size_data_v), data_W00_w(size_data_w)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Find vector
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, indi_u, indj_u, data_W00_u, &
                        size_data_v, indi_v, indj_v, data_W00_v, &
                        size_data_w, indi_w, indj_w, data_W00_w, &
                        source_coefs, source_vector)

end subroutine wq_get_source_3d

subroutine wq_get_force_3d(nb_cols_total, force_coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w, &
                            data_W00_u, data_W00_v, data_W00_w, &
                            force_vector)
    !! Computes force vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: force_coefs
    dimension :: force_coefs(nb_cols_total, 3)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W00_v, data_W00_w
    dimension :: data_W00_u(size_data_u), data_W00_v(size_data_v), data_W00_w(size_data_w)

    integer, parameter :: dimen = 3
    double precision, intent(out) :: force_vector
    dimension :: force_vector(nb_rows_u*nb_rows_v*nb_rows_w*dimen)

    ! Local data
    ! -------------
    integer :: nb_rows_total
    double precision, allocatable, dimension(:) :: force_vector_temp
    integer :: i, init, fin

    ! Initialize
    nb_rows_total = nb_rows_u*nb_rows_v*nb_rows_w*dimen
    allocate(force_vector_temp(nb_rows_total))

    do i = 1, dimen
            
        call wq_get_source_3D(nb_cols_total, force_coefs(:, i), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, data_W00_u, data_W00_v, data_W00_w, &
                        force_vector_temp)

        init = (i-1)*nb_rows_total + 1
        fin = init + nb_rows_total - 1
        force_vector(init : fin) = force_vector_temp   
    end do
    deallocate(force_vector_temp)

end subroutine wq_get_force_3d

! ----------------------------------------
! Assembly in 2D
! ----------------------------------------

subroutine wq_get_capacity_2d(nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B0_u, data_W00_u, data_B0_v, data_W00_v, &
                            size_data_I_u, size_data_I_v, & 
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
    double precision, intent(in) :: data_B0_u, data_W00_u, &
                                    data_B0_v, data_W00_v
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v)
    integer, intent(in) :: size_data_I_u, size_data_I_v

    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v+1), &
                    indj_result(size_data_I_u*size_data_I_v)

    ! Local data
    ! ---------------
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Initialize 
    size_data_result = size_data_I_u*size_data_I_v
    call csr_get_matrix_2d(capacity_coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B0_u, data_B0_v, data_W00_u, data_W00_v, &
                            size_data_I_u, size_data_I_v, & 
                            indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                            indi_result, size_data_result, data_result)
    deallocate(indj_I_u, indj_I_v)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_2d

subroutine wq_get_conductivity_2d(nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                                size_data_u, size_data_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B0_u, data_B1_u, &
                                data_W00_u, data_W01_u, &
                                data_W10_u, data_W11_u, &
                                data_B0_v, data_B1_v, &
                                data_W00_v, data_W01_v, &
                                data_W10_v, data_W11_v, &
                                size_data_I_u, size_data_I_v, & 
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
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
                                    data_B0_v, data_B1_v, &
                                    data_W00_v, data_W01_v, data_W10_v, data_W11_v
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_W00_u(size_data_u), data_W01_u(size_data_u), &
                    data_W10_u(size_data_u), data_W11_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_W00_v(size_data_v), data_W01_v(size_data_v), &
                    data_W10_v(size_data_v), data_W11_v(size_data_v)
    integer, intent(in) :: size_data_I_u, size_data_I_v 

    double precision, intent(out) :: data_result
    dimension :: data_result(size_data_I_u*size_data_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nb_rows_u*nb_rows_v+1), &
                    indj_result(size_data_I_u*size_data_I_v)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(size_data_I_u*size_data_I_v)
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v))
    call get_I_csr(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nb_rows_v, nb_rows_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_rows_u, nb_rows_u, size_data_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Initialize 
    size_data_result = size_data_I_u*size_data_I_v
    
    ! Get values
    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_v x B1_u (Kronecker product)
    ! Get W = W00_v x W11_u (Kronecker produt)

    call csr_get_matrix_2d(cond_coefs(1, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B1_u, data_B0_v, data_W11_u, data_W00_v, &
                        size_data_I_u, size_data_I_v, &
                        indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_v x W01_u (Kronecker produt)
    call csr_get_matrix_2d(cond_coefs(2, 1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B1_u, data_B0_v, data_W01_u, data_W10_v, &
                        size_data_I_u, size_data_I_v, &
                        indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B1_v x B0_u (Kronecker product)
    ! Get W = W01_v x W10_u (Kronecker produt)
    call csr_get_matrix_2d(cond_coefs(1, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B0_u, data_B1_v, data_W10_u, data_W01_v, &
                        size_data_I_u, size_data_I_v, &
                        indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL

    ! Get W = W11_v x W00_u (Kronecker produt)
    call csr_get_matrix_2d(cond_coefs(2, 2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B0_u, data_B1_v,  data_W00_u, data_W11_v, &
                        size_data_I_u, size_data_I_v, &
                        indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL
    
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
                            data_W00_u, data_W00_v, &
                            source_vector)
    !! Computes source vector in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_W00_u, data_W00_v
    dimension :: data_W00_u(size_data_u), data_W00_v(size_data_v)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_rows_u*nb_rows_v)

    ! Find vector
    call tensor2d_dot_vector_sp(nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, indi_u, indj_u, data_W00_u, &
                        size_data_v, indi_v, indj_v, data_W00_v, &
                        source_coefs, source_vector)

end subroutine wq_get_source_2d