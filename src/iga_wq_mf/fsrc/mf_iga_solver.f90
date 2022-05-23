! ====================================================
! module :: IGA - Matrix free methods 
! author :: Joaquin Cornejo
! ====================================================
subroutine iga_find_conductivity_diagonal_3d(nb_rows_total, nb_cols_total, cond_coefs, &
                                        nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                        size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
                                        Kdiag)
    !! Find the diagonal of conductivity matrix
    !! Indexes in CSR format
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in):: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)
    ! Csr format
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), indi_w(nb_rows_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension :: indj_u(size_data_u), indj_v(size_data_v), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B0_v, data_B0_w, &
                                    data_B1_u, data_B1_v, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w)
    
    double precision, intent(out) :: Kdiag
    dimension :: Kdiag(nb_rows_total)

    ! Local data
    ! ----------------------
    ! Find diagonal
    integer :: i
    double precision :: Kdiagtemp
    dimension :: Kdiagtemp(nb_rows_total)
    double precision, allocatable, dimension(:) :: data_W00_u, data_W11_u, data_W00_v, data_W11_v, data_W00_w, data_W11_w

    ! Initialize
    Kdiag = 0.d0

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
    
    !! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    ! ---------------------
    ! Get W = W = W00_w x W00_v x W11_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 1,:), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W11_u, data_W00_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 1,:), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W00_u, data_W11_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 1,:), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B1_u, data_B0_v, data_B0_w, data_W00_u, data_W00_v, data_W11_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! ---------------------
    ! Get W = W00_w x W00_v x W11_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W11_u, data_W00_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W11_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 2, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B1_v, data_B0_w, data_W00_u, data_W00_v, data_W11_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! ---------------------
    ! Get W = W00_w x W00_v x W11_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(1, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W11_u, data_W00_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(2, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W11_v, data_W00_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call find_physical_diag_3d(cond_coefs(3, 3, :), nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
    size_data_u, size_data_v, size_data_w, &
    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
    data_B0_u, data_B0_v, data_B1_w, data_W00_u, data_W00_v, data_W11_w, Kdiagtemp)

    Kdiag = Kdiag + Kdiagtemp

end subroutine iga_find_conductivity_diagonal_3d

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------
subroutine iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, coefs, array_in, array_out)

    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) ::  nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)
    double precision, intent(in) :: coefs
    dimension :: coefs(nb_cols_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nb_cols_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nb_cols_total)

    ! Local data
    ! ---------------
    integer :: nb_tasks, j1, j2, j3, genPos

    ! Initialize
    array_out = 0.d0 

    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(array_in)/nb_tasks) 
    ! Initialize coefficients
    do j3 = 1, nb_cols_w
        do j2 = 1, nb_cols_v
            do j1 = 1, nb_cols_u
                genPos = j1 + (j2-1)*nb_cols_u + (j3-1)*nb_cols_u*nb_cols_v
                array_out(genPos) = coefs(genPos)*array_in(genPos)*W_u(j1)*W_v(j2)*W_w(j3)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

end subroutine iga_diagonal_dot_vector

subroutine mf_iga_get_cu_3d(nb_rows_total, nb_cols_total, capacity_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                            indi_BT_u, indj_BT_u, indi_BT_v, indj_BT_v, indi_BT_w, indj_BT_w, &
                            data_B0T_u, data_B0T_v, data_B0T_w, &
                            indi_B_u, indj_B_u, indi_B_v, indj_B_v, indi_B_w, indj_B_w, &
                            data_B0_u, data_B0_v, data_B0_w, &
                            array_input, array_output)
    !! Computes capacity matrix in 3D case
    !! Indexes must be in CSR format
    
    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_cols_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)

    ! Csr format
    integer, intent(in) :: indi_BT_u, indi_BT_v, indi_BT_w
    dimension ::    indi_BT_u(nb_cols_u+1), &
                    indi_BT_v(nb_cols_v+1), &
                    indi_BT_w(nb_cols_w+1)
    integer, intent(in) :: indj_BT_u, indj_BT_v, indj_BT_w
    dimension ::    indj_BT_u(size_data_u), &
                    indj_BT_v(size_data_v), &
                    indj_BT_w(size_data_w)
    double precision, intent(in) :: data_B0T_u, data_B0T_v, data_B0T_w
    dimension :: data_B0T_u(size_data_u), data_B0T_v(size_data_v), data_B0T_w(size_data_w)
    integer, intent(in) :: indi_B_u, indi_B_v, indi_B_w
    dimension ::    indi_B_u(nb_rows_u+1), &
                    indi_B_v(nb_rows_v+1), &
                    indi_B_w(nb_rows_w+1)
    integer, intent(in) ::  indj_B_u, indj_B_v, indj_B_w
    dimension ::    indj_B_u(size_data_u), &
                    indj_B_v(size_data_v), &
                    indj_B_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B0_v, data_B0_w
    dimension :: data_B0_u(size_data_u), data_B0_v(size_data_v), data_B0_w(size_data_w)

    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_total)

    ! Local data 
    ! ------------------
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1tt

    ! Initialize
    allocate(array_temp_1(nb_cols_total))
    array_temp_1 = 0.d0

    ! Eval B.transpose * array_in
    call tensor3d_sparsedot_vector(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_BT_u, indj_BT_u, data_B0T_u, & 
    size_data_v, indi_BT_v, indj_BT_v, data_B0T_v, size_data_w, indi_BT_w, indj_BT_w, &
    data_B0T_w, array_input, array_temp_1)

    ! Evaluate diag(coefs) * array_temp1
    allocate(array_temp_1tt(nb_cols_total))
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            capacity_coefs, array_temp_1, array_temp_1tt)
    deallocate(array_temp_1)

    ! Eval W * array_temp1
    array_output = 0.d0
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)
    deallocate(array_temp_1tt)

end subroutine mf_iga_get_cu_3d

subroutine mf_iga_get_ku_3d(nb_rows_total, nb_cols_total, cond_coefs, &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                            indi_BT_u, indj_BT_u, indi_BT_v, indj_BT_v, indi_BT_w, indj_BT_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_B_u, indj_B_u, indi_B_v, indj_B_v, indi_B_w, indj_B_w, &
                            data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
                            array_input, array_output)
    !! Computes K.u in 3D case
    !! Indexes must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_cols_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)

    ! Csr format
    integer, intent(in) :: indi_BT_u, indi_BT_v, indi_BT_w
    dimension ::    indi_BT_u(nb_cols_u+1), &
                    indi_BT_v(nb_cols_v+1), &
                    indi_BT_w(nb_cols_w+1)
    integer, intent(in) :: indj_BT_u, indj_BT_v, indj_BT_w
    dimension ::    indj_BT_u(size_data_u), &
                    indj_BT_v(size_data_v), &
                    indj_BT_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w, &
                        data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B0T_u(size_data_u), data_B0T_v(size_data_v), data_B0T_w(size_data_w), &
                    data_B1T_u(size_data_u), data_B1T_v(size_data_v), data_B1T_w(size_data_w)

    integer, intent(in) :: indi_B_u, indi_B_v, indi_B_w
    dimension ::    indi_B_u(nb_rows_u+1), &
                    indi_B_v(nb_rows_v+1), &
                    indi_B_w(nb_rows_w+1)
    integer, intent(in) ::  indj_B_u, indj_B_v, indj_B_w
    dimension ::    indj_B_u(size_data_u), &
                    indj_B_v(size_data_v), &
                    indj_B_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B0_v, data_B0_w, &
                                    data_B1_u, data_B1_v, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B0_v(size_data_v), data_B0_w(size_data_w), &
                    data_B1_u(size_data_u), data_B1_v(size_data_v), data_B1_w(size_data_w)
    
    double precision, intent(in) :: array_input
    dimension :: array_input(nb_rows_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(nb_rows_total)

    ! Local data 
    ! ------------------
    double precision, allocatable, dimension(:) :: array_temp_1, array_temp_1tt

    ! Initialize
    array_output = 0.d0 

    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)

    ! Initialize
    allocate(array_temp_1(nb_cols_total))
    array_temp_1 = 0.d0

    ! Eval B.transpose * array_in
    call tensor3d_sparsedot_vector(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_BT_u, indj_BT_u, data_B1T_u, & 
    size_data_v, indi_BT_v, indj_BT_v, data_B0T_v, size_data_w, indi_BT_w, indj_BT_w,  &
    data_B0T_w, array_input, array_temp_1)

    ! ---------------------
    ! Get W = B0_w x B0_v x B1_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(1, 1, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B1_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)

    ! ---------------------
    ! Get W = B0_w x B1_v x B0_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(2, 1, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B1_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)

    ! ---------------------
    ! Get W = B1_w x B0_v x B0_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(3, 1, :), array_temp_1, array_temp_1tt)

    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B1_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)
    deallocate(array_temp_1)

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! Initialize
    allocate(array_temp_1(nb_cols_total))
    array_temp_1 = 0.d0

    ! Eval B.transpose * array_in
    call tensor3d_sparsedot_vector(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_BT_u, indj_BT_u, data_B0T_u, & 
    size_data_v, indi_BT_v, indj_BT_v, data_B1T_v, size_data_w, indi_BT_w, indj_BT_w,  &
    data_B0T_w, array_input, array_temp_1)
    
    ! ---------------------
    ! Get W = B0_w x B0_v x B1_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(1, 2, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B1_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)

    ! ---------------------
    ! Get W = B0_w x B1_v x B0_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(2, 2, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B1_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)
    
    ! ---------------------
    ! Get W = B1_w x B0_v x B0_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(3, 2, :), array_temp_1, array_temp_1tt) 

    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B1_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)
    deallocate(array_temp_1)

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! Initialize
    allocate(array_temp_1(nb_cols_total))
    array_temp_1 = 0.d0

    ! Eval B.transpose * array_in
    call tensor3d_sparsedot_vector(nb_cols_u, nb_rows_u, &
    nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, size_data_u, indi_BT_u, indj_BT_u, data_B0T_u, & 
    size_data_v, indi_BT_v, indj_BT_v, data_B0T_v, size_data_w, indi_BT_w, indj_BT_w,  &
    data_B1T_w, array_input, array_temp_1)
    
    ! ---------------------
    ! Get W = B0_w x B0_v x B1_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(1, 3, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B1_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)

    ! ---------------------
    ! Get W = B0_w x B1_v x B0_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(2, 3, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B1_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B0_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)

    ! ---------------------
    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    allocate(array_temp_1tt(nb_cols_total))
    array_temp_1tt = 0.d0

    ! Evaluate diag(coefs) * array_temp1
    call iga_diagonal_dot_vector(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, W_u, W_v, W_w, &
                            cond_coefs(3, 3, :), array_temp_1, array_temp_1tt)
    
    ! Eval W * array_temp1
    call tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, size_data_u, indi_B_u, indj_B_u, data_B0_u, &
    size_data_v, indi_B_v, indj_B_v, data_B0_v, size_data_w, indi_B_w, indj_B_w, & 
    data_B1_w, array_temp_1tt, array_output)

    deallocate(array_temp_1tt)
    deallocate(array_temp_1)
    
end subroutine mf_iga_get_ku_3d

subroutine mf_iga_get_ku_3d_csr(nb_rows_total, nb_cols_total, cond_coefs, &
                                nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_B1_u, W_u, &
                                data_B0_v, data_B1_v, W_v, &
                                data_B0_w, data_B1_w, W_w, &
                                array_input, array_output)
    
    !! Computes K.u in 3D case
    !! Indexes must be in CSR format
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
    double precision, intent(in) :: data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nb_cols_u), W_v(nb_cols_v), W_w(nb_cols_w)
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

    call mf_iga_get_ku_3d(nb_rows_total, nb_cols_total, cond_coefs, &
                        nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
                        array_input, array_output)
    
end subroutine mf_iga_get_ku_3d_csr

! ----------------------------------------
! Conjugate gradient
! ----------------------------------------
subroutine iga_mf_cg_3d(nb_rows_total, nb_cols_total, coefs, &
                        nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B1_u, W_u, &
                        data_B0_v, data_B1_v, W_v, &
                        data_B0_w, data_B1_w, W_w, &
                        b, nbIterations, epsilon, & 
                        Method, Jacob, directsol, x, RelRes, RelError)
    !! Conjugate gradient with ot without preconditioner 
    !! CSR FORMAT
                        
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nb_cols_total)
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

    character(len = 10) :: Method
    integer, intent(in) :: nbIterations
    double precision, intent(in) :: epsilon
    double precision, intent(in) :: b, Jacob, directsol
    dimension :: b(nb_rows_total), &
                Jacob(3, 3, nb_cols_u*nb_cols_v*nb_cols_w), &
                directsol(nb_rows_total)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nb_rows_total), &
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
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: data_W00_u, data_W11_u, data_W00_v, data_W11_v, data_W00_w, data_W11_w
    double precision, dimension(:), allocatable :: preconddiag, matrixdiag
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Deigen
    double precision :: Lu, Lv, Lw
    integer :: iter, i

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
                call mf_iga_get_Ku_3D(nb_rows_total, nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
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
        Lu = 1.d0
        Lv = 1.d0
        Lw = 1.d0

        if ((Method.eq.'TDS').or.(Method.eq.'TD')) then 
            ! --------------------------------------------
            ! DIAGONAL DECOMPOSITION
            ! -------------------------------------------- 
            allocate(Mcoef_u(nb_cols_u), Kcoef_u(nb_cols_u), &
                    Mcoef_v(nb_cols_v), Kcoef_v(nb_cols_v), &
                    Mcoef_w(nb_cols_w), Kcoef_w(nb_cols_w))

            ! Initialize coefficients 
            Mcoef_u = 1.d0
            Kcoef_u = 1.d0
            Mcoef_v = 1.d0
            Kcoef_v = 1.d0
            Mcoef_w = 1.d0
            Kcoef_w = 1.d0

            do iter = 1, 2
                call tensor_decomposition_3d(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((Method.eq.'JMS').or.(Method.eq.'JM')) then 
            ! --------------------------------------------
            ! NEW METHOD
            ! -------------------------------------------- 
            call jacobien_mean_3d(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, Jacob, Lu, Lv, Lw)
        end if

        ! --------------------------------------------
        ! EIGEN DECOMPOSITION
        ! -------------------------------------------- 
        allocate(U_u(nb_rows_u, nb_rows_u), D_u(nb_rows_u))
        allocate(U_v(nb_rows_v, nb_rows_v), D_v(nb_rows_v))
        allocate(U_w(nb_rows_w, nb_rows_w), D_w(nb_rows_w))
        
        allocate(data_W00_u(size_data_u), data_W11_u(size_data_u))
        allocate(Kdiag_u(nb_rows_u), Mdiag_u(nb_rows_u))
        do i = 1, size_data_u
            data_W00_u(i) = data_B0_u(i) * W_u(indj_u(i))
            data_W11_u(i) = data_B1_u(i) * W_u(indj_u(i))
        end do
        call eigen_decomposition(nb_rows_u, nb_cols_u, Mcoef_u, Kcoef_u, size_data_u, &
                                indi_u, indj_u, data_B0_u, data_W00_u, data_B1_u, &
                                data_W11_u, Method, D_u, U_u, Kdiag_u, Mdiag_u)
        deallocate(data_W00_u, data_W11_u)
        
        allocate(data_W00_v(size_data_v), data_W11_v(size_data_v))
        allocate(Kdiag_v(nb_rows_v), Mdiag_v(nb_rows_v))
        do i = 1, size_data_v
            data_W00_v(i) = data_B0_v(i) * W_v(indj_v(i))
            data_W11_v(i) = data_B1_v(i) * W_v(indj_v(i))
        end do
        call eigen_decomposition(nb_rows_v, nb_cols_v, Mcoef_v, Kcoef_v, size_data_v, &
                                indi_v, indj_v, data_B0_v, data_W00_v, data_B1_v, &
                                data_W11_v, Method, D_v, U_v, Kdiag_v, Mdiag_v)    
        deallocate(data_W00_v, data_W11_v)

        allocate(data_W00_w(size_data_w), data_W11_w(size_data_w))
        allocate(Kdiag_w(nb_rows_w), Mdiag_w(nb_rows_w))
        do i = 1, size_data_w
            data_W00_w(i) = data_B0_w(i) * W_w(indj_w(i))
            data_W11_w(i) = data_B1_w(i) * W_w(indj_w(i))
        end do
        call eigen_decomposition(nb_rows_w, nb_cols_w, Mcoef_w, Kcoef_w, size_data_w, &
                                indi_w, indj_w, data_B0_w, data_W00_w, data_B1_w, &
                                data_W11_w, Method, D_w, U_w, Kdiag_w, Mdiag_w)  
        deallocate(data_W00_w, data_W11_w)

        ! Find diagonal of eigen values
        allocate(I_u(nb_rows_u), I_v(nb_rows_v), I_w(nb_rows_w))
        allocate(Deigen(nb_rows_total))
        I_u = 1.d0
        I_v = 1.d0
        I_w = 1.d0
        call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, Lu, Lv, Lw, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((Method.eq.'TDS').or.(Method.eq.'TD')) then 
            deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        end if

        if ((Method.eq.'TDS').or.(Method.eq.'JMS')) then
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
            ! Find diagonal of preconditioner
            allocate(preconddiag(nb_rows_total))
            call find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, &
                                    Lu, Lv, Lw, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, preconddiag)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(matrixdiag(nb_rows_total))
            call iga_find_conductivity_diagonal_3D(nb_rows_total, nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
                                matrixdiag)
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
                call mf_iga_get_Ku_3D(nb_rows_total, nb_cols_total, coefs, nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_B0T_u, data_B1T_u, data_B0T_v, data_B1T_v, data_B0T_w, data_B1T_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_B1_u, data_B0_v, data_B1_v, data_B0_w, data_B1_w, &
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

end subroutine iga_mf_cg_3d
