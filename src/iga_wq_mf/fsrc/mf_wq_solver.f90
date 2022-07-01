! ====================================================
! module :: Weighted quadrature - Matrix free methods 
! author :: Joaquin Cornejo
! ====================================================

subroutine wq_find_conductivity_diagonal_3d(nb_cols_total, cond_coefs, &
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
                            Kdiag)
    
    use omp_lib
    use tensor_methods
    implicit none
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w

    ! Csr format
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension ::    indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), indi_w(nb_rows_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension ::    indj_u(size_data_u), indj_v(size_data_v), indj_w(size_data_w)
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

    double precision, intent(out) :: Kdiag
    dimension :: Kdiag(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Initialize
    Kdiag = 0.d0

    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    ! ---------------------
    ! Get W = W = W00_w x W00_v x W11_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 1, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W11_u, data_W00_v, data_W00_w, Kdiag)

    ! Get W = W00_w x W10_v x W01_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 1, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W01_u, data_W10_v, data_W00_w, Kdiag)

    ! Get W = W10_w x W00_v x W01_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 1, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W01_u, data_W00_v, data_W10_w, Kdiag)

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! ---------------------
    ! Get W = W00_w x W01_v x W10_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W10_u, data_W01_v, data_W00_w, Kdiag)

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W11_v, data_W00_w, Kdiag)

    ! Get W = W10_w x W01_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W01_v, data_W10_w, Kdiag)

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! ---------------------
    ! Get W = W01_w x W00_v x W10_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W10_u, data_W00_v, data_W01_w, Kdiag)

    ! Get W = W01_w x W10_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W10_v, data_W01_w, Kdiag)

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W00_v, data_W11_w, Kdiag)

end subroutine wq_find_conductivity_diagonal_3d

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------
subroutine wq_diagonal_dot_vector(nb_cols_total, coefs, array_in, array_out)
    
    use omp_lib
    implicit  none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    double precision, intent(in) :: coefs, array_in
    dimension :: coefs(nb_cols_total), array_in(nb_cols_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nb_cols_total)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_out = coefs * array_in
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    
end subroutine wq_diagonal_dot_vector

subroutine mf_wq_get_cu_3d( nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B0T_v, data_B0T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W00_u, data_W00_v, data_W00_w, &
                            array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    
    ! Csr format
    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision, intent(in) :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension ::    indi_u(nb_rows_u+1), &
                    indi_v(nb_rows_v+1), &
                    indi_w(nb_rows_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension ::    indj_u(size_data_u), &
                    indj_v(size_data_v), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W00_v, data_W00_w
    dimension ::    data_W00_u(size_data_u), &
                    data_W00_v(size_data_v), &
                    data_W00_w(size_data_w)

    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_u*nb_rows_v*nb_rows_w)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Local data 
    ! ----------------- 
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t
    double precision :: start1, finish1, start2, finish2, start3, finish3

    ! Initialize
    allocate(array_temp_1(nb_cols_total))

    ! Eval B.transpose * array_in
    call cpu_time(start1)
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_T_u, indj_T_u, data_B0T_u, & 
    size_data_v, indi_T_v, indj_T_v, data_B0T_v, size_data_w, indi_T_w, indj_T_w,  &
    data_B0T_w, array_input, array_temp_1)
    call cpu_time(finish1)

    ! Evaluate diag(coefs) * array_temp1
    call cpu_time(start2)
    allocate(array_temp_1t(nb_cols_total))
    call wq_diagonal_dot_vector(size(array_temp_1), capacity_coefs, array_temp_1, array_temp_1t)
    deallocate(array_temp_1)
    call cpu_time(finish2)

    ! Eval W * array_temp1
    call cpu_time(start3)
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W00_u, &
    size_data_v, indi_v, indj_v, data_W00_v, size_data_w, indi_w, indj_w, & 
    data_W00_w, array_temp_1t, array_output)
    deallocate(array_temp_1t)
    call cpu_time(finish3)
    
    ! print*, finish1 - start1, finish2 -start2, finish3 - start3

end subroutine mf_wq_get_cu_3d

subroutine mf_wq_get_cu_3d_csr( nb_rows_total, nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                            array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    
    ! Csr format
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension ::    indi_u(nb_rows_u+1), &
                    indi_v(nb_rows_v+1), &
                    indi_w(nb_rows_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension ::    indj_u(size_data_u), &
                    indj_v(size_data_v), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, &
                                    data_B0_v, data_W00_v, &
                                    data_B0_w, data_W00_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w)

    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_total)

    ! Local data 
    ! ----------------- 
    ! Local data
    ! ------------------
    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)

    ! ====================================================
    ! Initialize B transpose in CSR format
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w, &
                    indj_T_w, indi_T_w)
    ! ====================================================

    call mf_wq_get_cu_3d(nb_cols_total, capacity_coefs, &
                    nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                    size_data_u, size_data_v, size_data_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_B0T_u, data_B0T_v, data_B0T_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W00_u, data_W00_v, data_W00_w, &
                    array_input, array_output)

end subroutine mf_wq_get_cu_3d_csr

subroutine mf_wq_get_ku_3d( nb_cols_total, cond_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w,&
                            data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
                            data_W00_v, data_W01_v, data_W10_v, data_W11_v, &
                            data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                            array_input, array_output)
    !! Computes K.u in 3D case
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision, intent(in) :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)
    double precision, intent(in):: data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B1T_u(size_data_u), &
                    data_B1T_v(size_data_v), &
                    data_B1T_w(size_data_w)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension ::    indi_u(nb_rows_u+1), &
                    indi_v(nb_rows_v+1), &
                    indi_w(nb_rows_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension ::    indj_u(size_data_u), &
                    indj_v(size_data_v), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W01_u, data_W00_v, data_W01_v, data_W00_w, data_W01_w
    dimension ::    data_W00_u(size_data_u), data_W01_u(size_data_u), &
                    data_W00_v(size_data_v), data_W01_v(size_data_v), &
                    data_W00_w(size_data_w), data_W01_w(size_data_w)
    double precision, intent(in) :: data_W10_u, data_W11_u, data_W10_v, data_W11_v, data_W10_w, data_W11_w
    dimension ::    data_W10_u(size_data_u), data_W11_u(size_data_u), &
                    data_W10_v(size_data_v), data_W11_v(size_data_v), &
                    data_W10_w(size_data_w), data_W11_w(size_data_w)

    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_u*nb_rows_v*nb_rows_w)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Local data 
    ! ------------------
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1t, array_temp_1tt

    ! Initialize
    allocate(array_temp_1(nb_cols_total))
    allocate(array_temp_1t(nb_cols_total))
    allocate(array_temp_1tt(nb_cols_total))

    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    ! Eval B.transpose * array_in
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_T_u, indj_T_u, data_B1T_u, & 
    size_data_v, indi_T_v, indj_T_v, data_B0T_v, size_data_w, indi_T_w, indj_T_w,  &
    data_B0T_w, array_input, array_temp_1)

    ! ---------------------
    ! Get W = W00_w x W00_v x W11_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(1, 1, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W11_u, &
    size_data_v, indi_v, indj_v, data_W00_v, size_data_w, indi_w, indj_w, & 
    data_W00_w, array_temp_1t, array_output)
        
    ! ---------------------
    ! Get W = W00_w x W10_v x W01_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(2, 1, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W01_u, &
    size_data_v, indi_v, indj_v, data_W10_v, size_data_w, indi_w, indj_w, & 
    data_W00_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ---------------------
    ! Get W = W10_w x W00_v x W01_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(3, 1, :), array_temp_1, array_temp_1t)

    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W01_u, &
    size_data_v, indi_v, indj_v, data_W00_v, size_data_w, indi_w, indj_w, & 
    data_W10_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! Eval B.transpose * array_in
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_T_u, indj_T_u, data_B0T_u, & 
    size_data_v, indi_T_v, indj_T_v, data_B1T_v, size_data_w, indi_T_w, indj_T_w,  &
    data_B0T_w, array_input, array_temp_1)
    
    ! ---------------------
    ! Get W = W00_w x W01_v x W10_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(1, 2, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W10_u, &
    size_data_v, indi_v, indj_v, data_W01_v, size_data_w, indi_w, indj_w, & 
    data_W00_w, array_temp_1t, array_temp_1tt)
    
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ---------------------
    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(2, 2, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W00_u, &
    size_data_v, indi_v, indj_v, data_W11_v, size_data_w, indi_w, indj_w, & 
    data_W00_w, array_temp_1t, array_temp_1tt)
    
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ---------------------
    ! Get W = W10_w x W01_v x W00_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(3, 2, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W00_u, &
    size_data_v, indi_v, indj_v, data_W01_v, size_data_w, indi_w, indj_w, & 
    data_W10_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! Eval B.transpose * array_in
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_T_u, indj_T_u, data_B0T_u, & 
    size_data_v, indi_T_v, indj_T_v, data_B0T_v, size_data_w, indi_T_w, indj_T_w,  &
    data_B1T_w, array_input, array_temp_1)
    
    ! ---------------------
    ! Get W = W01_w x W00_v x W10_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(1, 3, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W10_u, &
    size_data_v, indi_v, indj_v, data_W00_v, size_data_w, indi_w, indj_w, & 
    data_W01_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 
    
    ! ---------------------
    ! Get W = W01_w x W10_v x W00_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(2, 3, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W00_u, &
    size_data_v, indi_v, indj_v, data_W10_v, size_data_w, indi_w, indj_w, & 
    data_W01_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! ---------------------
    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    ! Evaluate diag(coefs) * array_temp1
    call wq_diagonal_dot_vector(size(array_temp_1), cond_coefs(3, 3, :), array_temp_1, array_temp_1t)
    
    ! Eval W * array_temp1
    call tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_u, indj_u, data_W00_u, &
    size_data_v, indi_v, indj_v, data_W00_v, size_data_w, indi_w, indj_w, & 
    data_W11_w, array_temp_1t, array_temp_1tt)

    !$OMP PARALLEL
    !$OMP WORKSHARE 
    array_output = array_output + array_temp_1tt
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    deallocate(array_temp_1tt)
    deallocate(array_temp_1t)
    deallocate(array_temp_1)
    
end subroutine mf_wq_get_ku_3d

subroutine mf_wq_get_ku_3d_csr( nb_rows_total, nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_B1_u, data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
                                data_B0_v, data_B1_v, data_W00_v, data_W01_v, data_W10_v, data_W11_v, &
                                data_B0_w, data_B1_w, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                                array_input, array_output)
    
    !! Computes K.u in 3D case
    !! Indices must be in CSR format
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
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
    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_total)

    ! Local data
    ! ------------------
    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)
    double precision :: data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B1T_u(size_data_u), &
                    data_B1T_v(size_data_v), &
                    data_B1T_w(size_data_w)

    ! ====================================================
    ! Initialize B transpose in CSR format
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w, &
                    indj_T_w, indi_T_w)
    
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B1_w, indj_w, indi_w, data_B1T_w, &
                    indj_T_w, indi_T_w)
    ! ====================================================
    call mf_wq_get_ku_3d(nb_cols_total, cond_coefs, &
                        nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w,&
                        data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
                        data_W00_v, data_W01_v, data_W10_v, data_W11_v, &
                        data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                        array_input, array_output)
    
end subroutine mf_wq_get_ku_3d_csr

! ----------------------------------------
! Conjugate gradient
! ----------------------------------------
subroutine test_precondfd(nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                    size_data_u, size_data_v, size_data_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_B0_u, data_B1_u, data_W00_u, data_W11_u, &
                    data_B0_v, data_B1_v, data_W00_v, data_W11_v, &
                    data_B0_w, data_B1_w, data_W00_w, data_W11_w)
    !! Solve Sylvester equation P s = r
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, data_W00_u, data_W11_u, &
                                    data_B0_v, data_B1_v, data_W00_v, data_W11_v, &
                                    data_B0_w, data_B1_w, data_W00_w, data_W11_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_W00_u(size_data_u), data_W11_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_W00_v(size_data_v), data_W11_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w), &
                    data_W00_w(size_data_w), data_W11_w(size_data_w)

    character(len=10) :: Method = 'C'

    ! Local data
    ! ------------------
    double precision, dimension(:), allocatable :: capacity_coefs
    double precision, dimension(:, :, :), allocatable :: cond_coefs
    double precision :: s, r1
    dimension :: s(nb_rows_u*nb_rows_v*nb_rows_w), r1(nb_rows_u*nb_rows_v*nb_rows_w)
    
    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Kdiag_dummy, Mdiag_dummy
    double precision, dimension(:), allocatable :: Mcoef_dummy, Kcoef_dummy
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w

    double precision, dimension(:), allocatable :: I_u, I_v, I_w, Deigen
    double precision :: start, finish
    integer :: nb_rows_total, nb_cols_total

    ! Initialize
    s = 1.d0
    nb_rows_total = nb_rows_u*nb_rows_v*nb_rows_w
    nb_cols_total = nb_cols_u*nb_cols_v*nb_cols_w

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    allocate(U_u(nb_rows_u, nb_rows_u), D_u(nb_rows_u))
    allocate(U_v(nb_rows_v, nb_rows_v), D_v(nb_rows_v))
    allocate(U_w(nb_rows_w, nb_rows_w), D_w(nb_rows_w))
    
    allocate(Kdiag_dummy(nb_rows_u), Mdiag_dummy(nb_rows_u))
    call eigen_decomposition(nb_rows_u, nb_cols_u, Mcoef_dummy, Kcoef_dummy, size_data_u, &
                            indi_u, indj_u, data_B0_u, data_W00_u, data_B1_u, &
                            data_W11_u, Method, D_u, U_u, Kdiag_dummy, Mdiag_dummy)
    deallocate(Kdiag_dummy, Mdiag_dummy)

    allocate(Kdiag_dummy(nb_rows_v), Mdiag_dummy(nb_rows_v))
    call eigen_decomposition(nb_rows_v, nb_cols_v, Mcoef_dummy, Kcoef_dummy, size_data_v, &
                            indi_v, indj_v, data_B0_v, data_W00_v, data_B1_v, &
                            data_W11_v, Method, D_v, U_v, Kdiag_dummy, Mdiag_dummy)    
    deallocate(Kdiag_dummy, Mdiag_dummy)

    allocate(Kdiag_dummy(nb_rows_w), Mdiag_dummy(nb_rows_w))
    call eigen_decomposition(nb_rows_w, nb_cols_w, Mcoef_dummy, Kcoef_dummy, size_data_w, &
                            indi_w, indj_w, data_B0_w, data_W00_w, data_B1_w, &
                            data_W11_w, Method, D_w, U_w, Kdiag_dummy, Mdiag_dummy)  
    deallocate(Kdiag_dummy, Mdiag_dummy)

    ! Find diagonal of eigen values
    allocate(I_u(nb_rows_u), I_v(nb_rows_v), I_w(nb_rows_w))
    allocate(Deigen(nb_rows_total))
    I_u = 1.d0
    I_v = 1.d0
    I_w = 1.d0
    call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, 1.d0, 1.d0, 1.d0, &
                            I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
    deallocate(I_u, I_v, I_w)

    ! ! =============================
    ! ! It is already optimized
    ! ! Do fast diagonalization direct method
    ! call cpu_time(start)
    ! call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
    !             U_u, U_v, U_w, Deigen, s, r1)
    ! call cpu_time(finish)
    ! print *, finish-start

    ! ! =============================
    ! ! It is already optimized
    ! call cpu_time(start)
    ! call tensor3d_dot_vector(nb_rows_u, nb_rows_u, nb_rows_v, nb_rows_v, nb_rows_w, nb_rows_w, &
    !             U_u, U_v, U_w, s, r1)
    ! call cpu_time(finish)
    ! print*, finish-start

    ! ! =============================
    ! allocate(capacity_coefs(nb_cols_total))
    ! capacity_coefs = 1.d0

    ! call cpu_time(start)
    ! call mf_wq_get_cu_3d_csr( nb_rows_total, nb_cols_total, capacity_coefs, &
    ! nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    ! size_data_u, size_data_v, size_data_w, &
    ! indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    ! data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
    ! s, r1)
    ! call cpu_time(finish)
    ! print*, finish-start

    ! =============================
    allocate(cond_coefs(3, 3, nb_cols_total))
    cond_coefs = 1.d0

    call cpu_time(start)
    call mf_wq_get_ku_3d_csr( nb_rows_total, nb_cols_total, cond_coefs, &
    nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_u, data_W00_u, data_W00_u, data_W11_u, data_W11_u, &
    data_B0_v, data_B1_v, data_W00_v, data_W00_v, data_W11_v, data_W11_v, &
    data_B0_w, data_B1_w, data_W00_w, data_W00_w, data_W11_w, data_W11_w, &
    s, r1)
    call cpu_time(finish)
    print*, finish-start

end subroutine test_precondfd

subroutine wq_mf_cg_3d(nb_rows_total, nb_cols_total, coefs, &
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
                        b, nbIterations, epsilon, & 
                        Method, size_cond, conductivity, &
                        Jacob, directsol, x, RelRes, RelError)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nb_cols_total)
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

    character(len=10) :: Method
    integer, intent(in) :: nbIterations
    integer, intent(in) :: size_cond
    double precision, intent(in) :: conductivity
    dimension :: conductivity(3, 3, size_cond)
    double precision, intent(in) :: epsilon
    double precision, intent(in) :: b, Jacob, directsol
    dimension :: b(nb_rows_total), &
                 Jacob(3, 3, nb_cols_total), &
                 directsol(nb_rows_total)

    double precision, intent(out) :: x, RelRes, RelError
    dimension ::    x(nb_rows_total), &
                    RelRes(nbIterations+1), &
                    RelError(nbIterations+1)

    ! Local data
    ! ------------------
    ! Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha
    double precision :: r, p, Ap, dummy
    dimension :: r(nb_rows_total), p(nb_rows_total), Ap(nb_rows_total), dummy(nb_rows_total)
    integer :: k

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: preconddiag, matrixdiag
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Deigen
    double precision :: Lu, Lv, Lw
    integer :: iter

    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w
    double precision, dimension(:), allocatable :: I_u, I_v, I_w

    ! Preconditioned conjugate gradient
    double precision :: z
    dimension :: z(nb_rows_total)

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)
    double precision :: data_B1T_u, data_B1T_v, data_B1T_w
    dimension :: data_B1T_u(size_data_u), &
                data_B1T_v(size_data_v), &
                data_B1T_w(size_data_w)

    double precision :: c_u, c_v, c_w
    double precision :: lambda1, lambda2, lambda3

    ! ====================================================
    ! Initialize B transpose in CSR format
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w, &
                    indj_T_w, indi_T_w)
    
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B1_w, indj_w, indi_w, data_B1T_w, &
                    indj_T_w, indi_T_w)
    ! ====================================================
    
    ! Initiate variables
    x = 0.d0
    RelRes = 0.d0
    RelError = 0.d0

    if (Method.eq.'WP') then 
        if (nbIterations.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            r = b
            p = r
            rsold = dot_product(r, r)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do k = 1, nbIterations
                ! Calculate Ann xn 
                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                        data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                        p, Ap)
                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                ! Set relative value of residual 
                RelRes(k+1) = maxval(abs(r))/maxval(abs(b))
                RelError(k+1) = maxval(abs(directsol - x))/maxval(abs(directsol))

                if (RelRes(k+1).le.epsilon) then 
                    exit
                end if
                                
                rsnew = dot_product(r, r)
                p = r + rsnew/rsold * p
                rsold = rsnew
            end do
        end if
    else  
        ! Dimensions
        c_u = 1.d0
        c_v = 1.d0
        c_w = 1.d0

        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            ! --------------------------------------------
            ! DIAGONAL DECOMPOSITION
            ! --------------------------------------------
            ! Initialize coefficients
            allocate(Mcoef_u(nb_cols_u), Kcoef_u(nb_cols_u), &
                    Mcoef_v(nb_cols_v), Kcoef_v(nb_cols_v), &
                    Mcoef_w(nb_cols_w), Kcoef_w(nb_cols_w))
            Mcoef_u = 1.d0; Kcoef_u = 1.d0
            Mcoef_v = 1.d0; Kcoef_v = 1.d0
            Mcoef_w = 1.d0; Kcoef_w = 1.d0

            do iter = 1, 2
                call tensor_decomposition_3d(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((Method.eq.'JMS').or.(Method.eq.'JMC')) then 
            ! --------------------------------------------
            ! NEW METHOD
            ! --------------------------------------------
            ! Find dimensions and conductivity
            call jacobien_mean_3d(nb_cols_u, nb_cols_v, nb_cols_w, nb_cols_total, Jacob, &
                                size_cond, conductivity, Lu, Lv, Lw, lambda1, lambda2, lambda3)
                    
            c_u = lambda1*Lv*Lw/Lu
            c_v = lambda2*Lw*Lu/Lv
            c_w = lambda3*Lu*Lv/Lw
        
        end if

        ! --------------------------------------------
        ! EIGEN DECOMPOSITION
        ! -------------------------------------------- 
        allocate(U_u(nb_rows_u, nb_rows_u), D_u(nb_rows_u))
        allocate(U_v(nb_rows_v, nb_rows_v), D_v(nb_rows_v))
        allocate(U_w(nb_rows_w, nb_rows_w), D_w(nb_rows_w))
        
        allocate(Kdiag_u(nb_rows_u), Mdiag_u(nb_rows_u))
        call eigen_decomposition(nb_rows_u, nb_cols_u, Mcoef_u, Kcoef_u, size_data_u, &
                                indi_u, indj_u, data_B0_u, data_W00_u, data_B1_u, &
                                data_W11_u, Method, D_u, U_u, Kdiag_u, Mdiag_u)

        allocate(Kdiag_v(nb_rows_v), Mdiag_v(nb_rows_v))
        call eigen_decomposition(nb_rows_v, nb_cols_v, Mcoef_v, Kcoef_v, size_data_v, &
                                indi_v, indj_v, data_B0_v, data_W00_v, data_B1_v, &
                                data_W11_v, Method, D_v, U_v, Kdiag_v, Mdiag_v)    

        allocate(Kdiag_w(nb_rows_w), Mdiag_w(nb_rows_w))
        call eigen_decomposition(nb_rows_w, nb_cols_w, Mcoef_w, Kcoef_w, size_data_w, &
                                indi_w, indj_w, data_B0_w, data_W00_w, data_B1_w, &
                                data_W11_w, Method, D_w, U_w, Kdiag_w, Mdiag_w)  

        ! Find diagonal of eigen values
        allocate(I_u(nb_rows_u), I_v(nb_rows_v), I_w(nb_rows_w))
        allocate(Deigen(nb_rows_total))
        I_u = 1.d0
        I_v = 1.d0
        I_w = 1.d0
        call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, c_u, c_v, c_w, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        end if

        if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then 
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
            ! Find diagonal of preconditioner
            allocate(preconddiag(nb_rows_total))
            call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, &
                                    c_u, c_v, c_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, preconddiag)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(matrixdiag(nb_rows_total))
            call wq_find_conductivity_diagonal_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, nb_rows_v, &
            nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, size_data_v, size_data_w, &
            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B0_u, data_B1_u, data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_B0_v, data_B1_v, &
            data_W00_v, data_W01_v, data_W10_v, data_W11_v, data_B0_w, data_B1_w, data_W00_w, data_W01_w, &
            data_W10_w, data_W11_w, matrixdiag)

        end if

        if (nbIterations.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            r = b
            dummy = r
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then 
                call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, dummy) 
            end if
            call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                        U_u, U_v, U_w, Deigen, dummy, z)
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, z) 
            end if
            p = z
            rsold = dot_product(r, z)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do k = 1, nbIterations
                ! Calculate Ann xn 
                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                            data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                            p, Ap)

                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                ! Set relative value of residual 
                RelRes(k+1) = maxval(abs(r))/maxval(abs(b))
                RelError(k+1) = maxval(abs(directsol - x))/maxval(abs(directsol))

                if (RelRes(k+1).le.epsilon) then 
                    exit
                end if
                
                ! Calculate z
                dummy = r
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, dummy)   
                end if
                call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                            U_u, U_v, U_w, Deigen, dummy, z)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, z) 
                end if
                rsnew = dot_product(r, z)
                                
                p = z + rsnew/rsold * p
                rsold = rsnew
            end do

        end if
    end if

end subroutine wq_mf_cg_3d

subroutine wq_mf_bicgstab_3d(nb_rows_total, nb_cols_total, coefs, &
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
                            b, nbIterations, epsilon, & 
                            Method, size_cond, conductivity, & 
                            Jacob, directsol, x, RelRes, RelError)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nb_cols_total)
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

    character(len=10), intent(in) :: Method
    integer, intent(in) :: nbIterations
    integer, intent(in) :: size_cond
    double precision, intent(in) :: conductivity
    dimension :: conductivity(3, 3, size_cond)
    double precision, intent(in) :: epsilon
    double precision, intent(in) :: b, Jacob, directsol
    dimension ::    b(nb_rows_total), &
                    Jacob(3, 3, nb_cols_total), &
                    directsol(nb_rows_total)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension ::    x(nb_rows_total), &
                    RelRes(nbIterations+1), &
                    RelError(nbIterations+1)

    ! Local data
    ! ------------------
    ! Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, Ap, s, As, dummy
    dimension ::    r(nb_rows_total), rhat(nb_rows_total), p(nb_rows_total), & 
                    Ap(nb_rows_total), As(nb_rows_total), s(nb_rows_total), dummy(nb_rows_total)
    integer :: k

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: preconddiag, matrixdiag
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Deigen
    integer :: iter
    double precision :: Lu, Lv, Lw

    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w
    double precision, dimension(:), allocatable :: I_u, I_v, I_w

    ! Preconditioned conjugate gradient
    double precision :: ptilde, Aptilde, stilde, Astilde
    dimension :: ptilde(nb_rows_total), Aptilde(nb_rows_total), Astilde(nb_rows_total), stilde(nb_rows_total)

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)
    double precision :: data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B1T_u(size_data_u), &
                    data_B1T_v(size_data_v), &
                    data_B1T_w(size_data_w)
                    
    double precision :: c_u, c_v, c_w
    double precision :: lambda1, lambda2, lambda3
    
    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w, &
                    indj_T_w, indi_T_w)
    
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B1_w, indj_w, indi_w, data_B1T_w, &
                    indj_T_w, indi_T_w)
    ! =======================================================

    ! Initiate variables
    x = 0.d0
    RelRes = 0.d0
    RelError = 0.d0

    if (Method.eq.'WP') then 
        if (nbIterations.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            r = b
            rhat = r
            p = r
            rsold = dot_product(r, rhat)

            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do k = 1, nbIterations
                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                        data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                        p, Ap)
                alpha = rsold/dot_product(Ap, rhat)
                s = r - alpha*Ap

                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                        data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                        s, As)
                omega = dot_product(As, s)/dot_product(As, As)
                x = x + alpha*p + omega*s
                r = s - omega*As

                RelRes(k+1) = maxval(abs(r))/maxval(abs(b))
                RelError(k+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                
                if (RelRes(k+1).le.epsilon) then 
                    exit
                end if

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew
            end do
        end if

    else  
        ! Dimensions
        c_u = 1.d0
        c_v = 1.d0
        c_w = 1.d0

        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            ! --------------------------------------------
            ! DIAGONAL DECOMPOSITION
            ! --------------------------------------------
            ! Initialize coefficients 
            allocate(Mcoef_u(nb_cols_u), Kcoef_u(nb_cols_u), &
                    Mcoef_v(nb_cols_v), Kcoef_v(nb_cols_v), &
                    Mcoef_w(nb_cols_w), Kcoef_w(nb_cols_w))            
            Mcoef_u = 1.d0; Kcoef_u = 1.d0
            Mcoef_v = 1.d0; Kcoef_v = 1.d0
            Mcoef_w = 1.d0; Kcoef_w = 1.d0

            do iter = 1, 2
                call tensor_decomposition_3d(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((Method.eq.'JMS').or.(Method.eq.'JMC')) then 
            ! --------------------------------------------
            ! NEW METHOD
            ! --------------------------------------------
            ! Find dimensions and conductivity
            call jacobien_mean_3d(nb_cols_u, nb_cols_v, nb_cols_w, nb_cols_total, Jacob, &
                                size_cond, conductivity, Lu, Lv, Lw, lambda1, lambda2, lambda3)
                                
            c_u = lambda1*Lv*Lw/Lu
            c_v = lambda2*Lw*Lu/Lv
            c_w = lambda3*Lu*Lv/Lw

        end if

        ! --------------------------------------------
        ! EIGEN DECOMPOSITION
        ! -------------------------------------------- 
        allocate(U_u(nb_rows_u, nb_rows_u), D_u(nb_rows_u))
        allocate(U_v(nb_rows_v, nb_rows_v), D_v(nb_rows_v))
        allocate(U_w(nb_rows_w, nb_rows_w), D_w(nb_rows_w))
        
        allocate(Kdiag_u(nb_rows_u), Mdiag_u(nb_rows_u))
        call eigen_decomposition(nb_rows_u, nb_cols_u, Mcoef_u, Kcoef_u, size_data_u, &
                                indi_u, indj_u, data_B0_u, data_W00_u, data_B1_u, &
                                data_W11_u, Method, D_u, U_u, Kdiag_u, Mdiag_u)

        allocate(Kdiag_v(nb_rows_v), Mdiag_v(nb_rows_v))
        call eigen_decomposition(nb_rows_v, nb_cols_v, Mcoef_v, Kcoef_v, size_data_v, &
                                indi_v, indj_v, data_B0_v, data_W00_v, data_B1_v, &
                                data_W11_v, Method, D_v, U_v, Kdiag_v, Mdiag_v)    

        allocate(Kdiag_w(nb_rows_w), Mdiag_w(nb_rows_w))
        call eigen_decomposition(nb_rows_w, nb_cols_w, Mcoef_w, Kcoef_w, size_data_w, &
                                indi_w, indj_w, data_B0_w, data_W00_w, data_B1_w, &
                                data_W11_w, Method, D_w, U_w, Kdiag_w, Mdiag_w)  

        ! Find diagonal of eigen values
        allocate(I_u(nb_rows_u), I_v(nb_rows_v), I_w(nb_rows_w))
        allocate(Deigen(nb_rows_total))
        I_u = 1.d0
        I_v = 1.d0
        I_w = 1.d0
        call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, c_u, c_v, c_w, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        end if

        if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
            ! Find diagonal of preconditioner
            allocate(preconddiag(nb_rows_total))
            call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, &
                                    c_u, c_v, c_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, preconddiag)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(matrixdiag(nb_rows_total))
            call wq_find_conductivity_diagonal_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, nb_rows_v, &
            nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, size_data_v, size_data_w, &
            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B0_u, data_B1_u, data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_B0_v, data_B1_v, &
            data_W00_v, data_W01_v, data_W10_v, data_W11_v, data_B0_w, data_B1_w, data_W00_w, data_W01_w, &
            data_W10_w, data_W11_w, matrixdiag)

        end if

        if (nbIterations.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            r = b
            rhat = r
            p = r
            dummy = p
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, dummy) 
            end if 
            call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                        U_u, U_v, U_w, Deigen, dummy, ptilde)
            if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, ptilde)  
            end if
            rsold = dot_product(r, rhat)
            RelRes(1) = 1.d0
            RelError(1) = 1.d0

            do k = 1, nbIterations
                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                            data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                            ptilde, Aptilde)

                alpha = rsold/dot_product(Aptilde, rhat)
                s = r - alpha*Aptilde

                dummy = s
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, dummy)  
                end if
                call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                            U_u, U_v, U_w, Deigen, dummy, stilde)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, stilde) 
                end if

                call mf_wq_get_Ku_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W00_u, data_W01_u, data_W10_u, data_W11_u, data_W00_v, data_W01_v, &
                            data_W10_v, data_W11_v, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
                            stilde, Astilde)

                omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
                x = x + alpha*ptilde + omega*stilde
                r = s - omega*Astilde    
                
                RelRes(k+1) = maxval(abs(r))/maxval(abs(b))
                RelError(k+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                
                if (RelRes(k+1).le.epsilon) then 
                    exit
                end if

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                
                dummy = p
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, dummy) 
                end if
                call fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                            U_u, U_v, U_w, Deigen, dummy, ptilde)
                if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
                    call scaling_FastDiag(nb_rows_total, preconddiag, matrixdiag, ptilde) 
                end if

                rsold = rsnew
            end do
        end if
    end if

end subroutine wq_mf_bicgstab_3d

subroutine wq_mf_interp_3d(nb_rows_total, nb_cols_total, coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_B1_u, &
                            data_W00_u, data_W11_u, &
                            data_B0_v, data_B1_v, &
                            data_W00_v, data_W11_v, &
                            data_B0_w, data_B1_w, &
                            data_W00_w, data_W11_w, &
                            b, nbIterations, epsilon, & 
                            x, RelRes)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_W00_u, data_W11_u, &
                                    data_B0_v, data_B1_v, &
                                    data_W00_v, data_W11_v, &
                                    data_B0_w, data_B1_w, &
                                    data_W00_w, data_W11_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_W00_u(size_data_u), data_W11_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_W00_v(size_data_v), data_W11_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w), &
                    data_W00_w(size_data_w), data_W11_w(size_data_w)

    character(len=10) :: Method = 'C'
    integer, intent(in) :: nbIterations
    double precision, intent(in) :: epsilon
    double precision, intent(in) :: b
    dimension :: b(nb_rows_total)
    
    double precision, intent(out) :: x, RelRes
    dimension :: x(nb_rows_total), RelRes(nbIterations+1)

    ! Local data
    ! ------------------
    ! Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, s, dummy
    dimension ::    r(nb_rows_total), rhat(nb_rows_total), p(nb_rows_total), & 
                    s(nb_rows_total), dummy(nb_rows_total)
    integer :: k

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w

    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w

    ! Preconditioned conjugate gradient
    double precision :: ptilde, Aptilde, stilde, Astilde
    dimension :: ptilde(nb_rows_total), Aptilde(nb_rows_total), Astilde(nb_rows_total), stilde(nb_rows_total)

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension ::    data_B0T_u(size_data_u), &
                    data_B0T_v(size_data_v), &
                    data_B0T_w(size_data_w)
    double precision :: data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B1T_u(size_data_u), &
                    data_B1T_v(size_data_v), &
                    data_B1T_w(size_data_w)

    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w, &
                    indj_T_w, indi_T_w)
    
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u, &
                    indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v, &
                    indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B1_w, indj_w, indi_w, data_B1T_w, &
                    indj_T_w, indi_T_w)
    ! =======================================================

    ! Initiate variables
    x = 0.d0
    RelRes = 0.d0

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    allocate(U_u(nb_rows_u, nb_rows_u), D_u(nb_rows_u))
    allocate(U_v(nb_rows_v, nb_rows_v), D_v(nb_rows_v))
    allocate(U_w(nb_rows_w, nb_rows_w), D_w(nb_rows_w))
    
    allocate(Kdiag_u(nb_rows_u), Mdiag_u(nb_rows_u))
    call eigen_decomposition(nb_rows_u, nb_cols_u, Mcoef_u, Kcoef_u, size_data_u, &
                            indi_u, indj_u, data_B0_u, data_W00_u, data_B1_u, &
                            data_W11_u, Method, D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nb_rows_v), Mdiag_v(nb_rows_v))
    call eigen_decomposition(nb_rows_v, nb_cols_v, Mcoef_v, Kcoef_v, size_data_v, &
                            indi_v, indj_v, data_B0_v, data_W00_v, data_B1_v, &
                            data_W11_v, Method, D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nb_rows_w), Mdiag_w(nb_rows_w))
    call eigen_decomposition(nb_rows_w, nb_cols_w, Mcoef_w, Kcoef_w, size_data_w, &
                            indi_w, indj_w, data_B0_w, data_W00_w, data_B1_w, &
                            data_W11_w, Method, D_w, U_w, Kdiag_w, Mdiag_w)  

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    r = b
    rhat = r
    p = r
    dummy = p
    call precond_interp_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                U_u, U_v, U_w, dummy, ptilde)
    rsold = dot_product(r, rhat)
    RelRes(1) = 1.d0

    do k = 1, nbIterations
        call mf_wq_get_Cu_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                    size_data_u, size_data_v, size_data_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_B0T_u, data_B0T_v, data_B0T_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W00_u, data_W00_v, data_W00_w, &
                    ptilde, Aptilde)

        alpha = rsold/dot_product(Aptilde, rhat)
        s = r - alpha*Aptilde

        dummy = s
        call precond_interp_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                    U_u, U_v, U_w, dummy, stilde)

        call mf_wq_get_Cu_3D(nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                    size_data_u, size_data_v, size_data_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_B0T_u, data_B0T_v, data_B0T_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W00_u, data_W00_v, data_W00_w, &
                    stilde, Astilde)

        omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
        x = x + alpha*ptilde + omega*stilde
        r = s - omega*Astilde    
        
        RelRes(k+1) = maxval(abs(r))/maxval(abs(b))
        
        if (RelRes(k+1).le.epsilon) then 
            exit
        end if

        rsnew = dot_product(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        
        dummy = p
        call precond_interp_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                    U_u, U_v, U_w, dummy, ptilde)
        rsold = rsnew
    end do

end subroutine wq_mf_interp_3d