! ==========================
! module :: assembly for IGA-WQ 
! author :: Joaquin Cornejo
! modules :: operateurs.f90 (MatrixInv and MatrixDet)
! ==========================

subroutine eval_thermal_coefficient(dime, nnz, JJ, KK, CC, Kcoef, Ccoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------    
    integer, intent(in) :: dime, nnz
    double precision, intent(in) :: JJ, KK
    dimension :: JJ(dime, dime, nnz), KK(dime, dime)
    double precision, intent(in) :: CC

    double precision, intent(out) :: Kcoef, Ccoef
    dimension :: Kcoef(dime, dime, nnz), Ccoef(nnz)

    ! Local data
    ! -----------
    integer :: i, nb_tasks
    double precision :: Jt, detJt, invJt
    dimension :: Jt(dime, dime), invJt(dime, dime)
    double precision :: Kt
    dimension :: Kt(dime, dime)  

    !$OMP PARALLEL PRIVATE(Jt, invJt, detJt, Kt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
    do i = 1, nnz
        ! Get individual jacobien
        Jt = JJ(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJt, Jt, detJt, dime)

        ! For K = invJ * prop * detJ * transpose(invJ) 
        Kt = detJt * matmul(invJt, KK) 
        Kcoef(:, :, i) = matmul(Kt, transpose(invJt))

        ! For C = detJ  * prop
        Ccoef(i) = detJt * CC
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_thermal_coefficient

subroutine eval_mech_coefficient(dime, nnz, JJ, DD, Scoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in):: dime, nnz
    double precision, intent(in) :: JJ, DD
    dimension :: JJ(dime, dime, nnz), DD(3*(dime-1), 3*(dime-1))
    double precision, intent(out) :: Scoef
    dimension :: Scoef(dime*dime, dime*dime, nnz)

    ! Local data
    ! -----------
    integer :: i, k, nb_tasks
    double precision :: Jt, detJt, invJt
    dimension :: Jt(dime, dime), invJt(dime, dime)
    double precision :: St1, St2
    dimension :: St1(dime*dime, dime*dime), St2(dime*dime, dime*dime)  

    double precision, allocatable, dimension(:, :) :: M
    double precision, allocatable, dimension(:, :) :: MDMt
    double precision :: MDM(dime*dime, dime*dime)
    double precision :: invJextend(dime*dime, dime*dime)

    ! Define MM transformation matrix 
    if (dime.eq.2) then 
        allocate(M(3, 4))
        allocate(MDMt(4, 3))
        M = 0.d0
        M(1, 1) = 1.d0
        M(2, 4) = 1.d0
        M(3, 2) = 1.d0
        M(3, 3) = 1.d0

    else if (dime.eq.3) then 
        allocate(M(6, 9))
        allocate(MDMt(9, 6))
        M = 0.d0
        M(1, 1) = 1.d0
        M(2, 5) = 1.d0
        M(3, 9) = 1.d0
        M(4, 2) = 1.d0
        M(4, 4) = 1.d0
        M(5, 6) = 1.d0
        M(5, 8) = 1.d0
        M(6, 3) = 1.d0
        M(6, 7) = 1.d0
    end if

    ! Evaluate MM.transpose * DD * MM
    MDMt = matmul(transpose(M), DD)
    MDM = matmul(MDMt, M)
    deallocate(M, MDMt)

    !$OMP PARALLEL PRIVATE(Jt, invJt, detJt, invJextend, St1, St2, k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
    do i = 1, nnz

        ! Initialize
        invJextend = 0.d0

        ! Get individual jacobien
        Jt = JJ(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJt, Jt, detJt, dime)

        do k = 1, dime
            invJextend((k-1)*dime+1:k*dime, (k-1)*dime+1:k*dime) = invJt
        end do
        
        St1 = matmul(transpose(invJextend), MDM)
        St2 = matmul(St1, invJextend) * detJt

        Scoef(:, :, i) = St2
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_mech_coefficient

subroutine eval_thermomech_coefficient(dime, nnz, JJ, DD, alpha, Tcoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dime, nnz
    double precision, intent(in) :: JJ, DD
    dimension :: JJ(dime, dime, nnz), DD(3*(dime-1), 3*(dime-1))
    double precision, intent(in) :: alpha
            
    double precision, intent(out) :: Tcoef
    dimension :: Tcoef(dime*dime, nnz)

    ! Local data
    ! -----------
    integer :: i, k, nb_tasks
    double precision :: Jt, detJt, invJt
    dimension :: Jt(dime, dime), invJt(dime, dime)
    double precision :: Rt
    dimension :: Rt(dime*dime)

    double precision, allocatable, dimension(:, :) :: Mt
    double precision :: MDV(dime*dime), VV(3*(dime-1)), DV(3*(dime-1))
    double precision :: invJextend(dime*dime, dime*dime)

    ! Define MM transformation matrix 
    if (dime.eq.2) then 
        allocate(Mt(3, 4))
        Mt = 0.d0
        Mt(1, 1) = 1.d0
        Mt(2, 4) = 1.d0
        Mt(3, 2) = 1.d0
        Mt(3, 3) = 1.d0

        VV = 0.d0
        VV(1) = 1.d0
        VV(2) = 1.d0

    else if (dime.eq.3) then 
        allocate(Mt(6, 9))
        Mt = 0.d0
        Mt(1, 1) = 1.d0
        Mt(2, 5) = 1.d0
        Mt(3, 9) = 1.d0
        Mt(4, 2) = 1.d0
        Mt(4, 4) = 1.d0
        Mt(5, 6) = 1.d0
        Mt(5, 8) = 1.d0
        Mt(6, 3) = 1.d0
        Mt(6, 7) = 1.d0

        VV = 0.d0
        VV(1) = 1.d0
        VV(2) = 1.d0
        VV(3) = 1.d0
    end if

    ! Evaluate MM.transpose * DD * MM
    DV = matmul(DD, VV)
    MDV = matmul(transpose(Mt), DV)
    deallocate(Mt)

    !$OMP PARALLEL PRIVATE(Jt, invJt, detJt, invJextend, Rt, k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
    do i = 1, nnz

        ! Initialize
        invJextend = 0.d0

        ! Get individual jacobien
        Jt = JJ(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJt, Jt, detJt, dime)

        do k = 1, dime
            invJextend((k-1)*dime+1:k*dime, (k-1)*dime+1:k*dime) = invJt
        end do
        
        Rt = matmul(transpose(invJextend), MDV) * detJt * alpha
        Tcoef(:, i) = Rt
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_thermomech_coefficient

subroutine jacobien_physicalposition_3d(nb_rows_total, nb_cols_total, &
                                        nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        nb_rows_w, nb_cols_w, &
                                        size_data_u, size_data_v, size_data_w, &
                                        ctrlpts_x, ctrlpts_y, ctrlpts_z, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B0_u, data_B1_u, &
                                        data_B0_v, data_B1_v, &
                                        data_B0_w, data_B1_w, &
                                        jacob, physical_pos, detJ)
    !! Computes jacobien in 3D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_rows_v, nb_rows_w, &
                            nb_cols_u, nb_cols_v, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_B0_v, data_B1_v, &
                                    data_B0_w, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w)
    double precision, intent(in) :: ctrlpts_x, ctrlpts_y, ctrlpts_z
    dimension :: ctrlpts_x(nb_rows_total), ctrlpts_y(nb_rows_total), ctrlpts_z(nb_rows_total)

    double precision, intent(out) :: jacob
    dimension :: jacob(3, 3, nb_cols_total)

    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(3, nb_cols_total)

    double precision, intent(out) :: detJ
    dimension :: detJ(nb_cols_total)

    ! Local data
    !-----------------
    double precision :: result_temp
    dimension ::  result_temp(nb_cols_total)
    integer :: nb_tasks, i
    double precision :: detJt

    ! Csr format (Transpose)
    integer ::  indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer ::  indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_B0T_u, data_B0T_v, data_B0T_w, &
                        data_B1T_u, data_B1T_v, data_B1T_w
    dimension ::    data_B0T_u(size_data_u), data_B0T_v(size_data_v), data_B0T_w(size_data_w), &
                    data_B1T_u(size_data_u), data_B1T_v(size_data_v), data_B1T_w(size_data_w)

    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, &
                data_B0T_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, &
                data_B0T_v, indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B0_w, indj_w, indi_w, &
                data_B0T_w, indj_T_w, indi_T_w)
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, &
                data_B1T_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, &
                data_B1T_v, indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B1_w, indj_w, indi_w, &
                data_B1T_w, indj_T_w, indi_T_w)
    ! ====================================================
    ! ---------------------------------------------------
    ! For J00, J10 and J20
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B1T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_x, result_temp)

    jacob(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B1T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_y, result_temp)
    jacob(2, 1, :) = result_temp

    ! Compute B.Transpose . CP_z
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B1T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_z, result_temp)
    jacob(3, 1, :) = result_temp

    ! ---------------------------------------------------
    ! For J01, J11, and J21
    ! ---------------------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)

    ! Compute B.Transpose . CP_x
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B1T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_x, result_temp)
    jacob(1, 2, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B1T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_y, result_temp)
    jacob(2, 2, :) = result_temp

    ! Compute B.Transpose . CP_z
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B1T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_z, result_temp)
    jacob(3, 2, :) = result_temp

    ! ---------------------------------------------------
    ! For J02, J12, and J22
    ! ---------------------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B1T_w, &
                                ctrlpts_x, result_temp)
    jacob(1, 3, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B1T_w, &
                                ctrlpts_y, result_temp)
    jacob(2, 3, :) = result_temp

    ! Compute B.Transpose . CP_z
    result_temp = 0.d0
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B1T_w, &
                                ctrlpts_z, result_temp)
    jacob(3, 3, :) = result_temp

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_x, result_temp)
    physical_pos(1, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_y, result_temp)
    physical_pos(2, :) = result_temp

    ! Compute B.Transpose . CP_z
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                size_data_w, indi_T_w, indj_T_w, data_B0T_w, &
                                ctrlpts_z, result_temp)
    physical_pos(3, :) = result_temp

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_cols_total/nb_tasks) 
    do i = 1, nb_cols_total
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, i), detJt, 3)

        ! Assign values
        detJ(i) = detJt
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_3d

subroutine interpolation_3d(nb_rows_total, nb_cols_total, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            ctrlpts, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, &
                            interpolation)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_rows_v, nb_rows_w, &
                            nb_cols_u, nb_cols_v, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension ::    data_B_u(size_data_u), &
                    data_B_v(size_data_v), &
                    data_B_w(size_data_w)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(nb_rows_total)

    double precision, intent(out) :: interpolation
    dimension :: interpolation(nb_cols_total)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer ::  indi_T_u, indi_T_v, indi_T_w
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1), &
                    indi_T_w(nb_cols_w+1)
    integer ::  indj_T_u, indj_T_v, indj_T_w
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v), &
                    indj_T_w(size_data_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension ::    data_BT_u(size_data_u), &
                    data_BT_v(size_data_v), &
                    data_BT_w(size_data_w)

    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(nb_rows_w, nb_cols_w, size_data_w, data_B_w, indj_w, indi_w, &
                data_BT_w, indj_T_w, indi_T_w)
    ! ====================================================

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B_w x B_v x B_u (Kronecker product)
    ! Compute B.Transpose . CP
    call tensor3d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, nb_cols_w, nb_rows_w, &
                                size_data_u, indi_T_u, indj_T_u, data_BT_u, &
                                size_data_v, indi_T_v, indj_T_v, data_BT_v, &
                                size_data_w, indi_T_w, indj_T_w, data_BT_w, &
                                ctrlpts, interpolation)

end subroutine interpolation_3d

subroutine jacobien_physicalposition_2d(nb_rows_total, nb_cols_total, &
                                        nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        size_data_u, size_data_v, &
                                        ctrlpts_x, ctrlpts_y, &
                                        indi_u, indj_u, indi_v, indj_v, &
                                        data_B0_u, data_B1_u, &
                                        data_B0_v, data_B1_v, &
                                        jacob, physical_pos, detJ)
    !! Computes jacobien in 2D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_rows_v, &
                            nb_cols_u, nb_cols_v, &
                            size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_B0_v, data_B1_v
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v)
    double precision, intent(in) :: ctrlpts_x, ctrlpts_y
    dimension ::    ctrlpts_x(nb_rows_total), &
                    ctrlpts_y(nb_rows_total)

    double precision, intent(out) :: jacob
    dimension ::  jacob(2, 2, nb_cols_total)

    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(2, nb_cols_total)

    double precision, intent(out) :: detJ
    dimension :: detJ(nb_cols_total)

    ! Local data
    !-----------------
    double precision :: result_temp
    dimension ::  result_temp(nb_cols_total)
    integer :: nb_tasks, i
    double precision :: detJt

    ! Csr format (Transpose)
    integer ::  indi_T_u, indi_T_v
    dimension ::    indi_T_u(nb_cols_u+1), &
                    indi_T_v(nb_cols_v+1)
    integer ::  indj_T_u, indj_T_v
    dimension ::    indj_T_u(size_data_u), &
                    indj_T_v(size_data_v)
    double precision :: data_B0T_u, data_B0T_v, &
                        data_B1T_u, data_B1T_v
    dimension ::    data_B0T_u(size_data_u), data_B0T_v(size_data_v), &
                    data_B1T_u(size_data_u), data_B1T_v(size_data_v)

    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B0_u, indj_u, indi_u, &
                data_B0T_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B0_v, indj_v, indi_v, &
                data_B0T_v, indj_T_v, indi_T_v)
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B1_u, indj_u, indi_u, &
                data_B1T_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B1_v, indj_v, indi_v, &
                data_B1T_v, indj_T_v, indi_T_v)
    ! ====================================================
    ! ---------------------------------------------------
    ! For J00, J10 and J20
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B1T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                ctrlpts_x, result_temp)

    jacob(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B1T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                ctrlpts_y, result_temp)
    jacob(2, 1, :) = result_temp

    ! ---------------------------------------------------
    ! For J01, J11, and J21
    ! ---------------------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)

    ! Compute B.Transpose . CP_x
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B1T_v, &
                                ctrlpts_x, result_temp)
    jacob(1, 2, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B1T_v, &
                                ctrlpts_y, result_temp)
    jacob(2, 2, :) = result_temp

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                ctrlpts_x, result_temp)
    physical_pos(1, :) = result_temp

    ! Compute B.Transpose . CP_y
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_B0T_u, &
                                size_data_v, indi_T_v, indj_T_v, data_B0T_v, &
                                ctrlpts_y, result_temp)
    physical_pos(2, :) = result_temp

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_cols_total/nb_tasks) 
    do i = 1, nb_cols_total
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, i), detJt, 2)

        ! Assign values
        detJ(i) = detJt
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_2d

subroutine interpolation_2d(nb_rows_total, nb_cols_total, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            size_data_u, size_data_v, &
                            ctrlpts, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            data_B_u, data_B_v, &
                            interpolation)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in) :: nb_rows_total, nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_rows_v, &
                            nb_cols_u, nb_cols_v, &
                            size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(size_data_u), data_B_v(size_data_v)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(nb_rows_total)

    double precision, intent(out) :: interpolation
    dimension :: interpolation(nb_cols_total)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nb_cols_u+1), indi_T_v(nb_cols_v+1)
    integer :: indj_T_u, indj_T_v
    dimension :: indj_T_u(size_data_u), indj_T_v(size_data_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(size_data_u), data_BT_v(size_data_v)

    ! ====================================================
    ! Initialize
    call csr2csc(nb_rows_u, nb_cols_u, size_data_u, data_B_u, indj_u, indi_u, &
                data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(nb_rows_v, nb_cols_v, size_data_v, data_B_v, indj_v, indi_v, &
                data_BT_v, indj_T_v, indi_T_v)
    ! ====================================================

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B_w x B_v x B_u (Kronecker product)  
    ! Compute B.Transpose . CP
    call tensor2d_dot_vector_sp(nb_cols_u, nb_rows_u, &
                                nb_cols_v, nb_rows_v, &
                                size_data_u, indi_T_u, indj_T_u, data_BT_u, &
                                size_data_v, indi_T_v, indj_T_v, data_BT_v, &
                                ctrlpts, interpolation)

end subroutine interpolation_2d

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

subroutine wq_get_advention_3d(nb_cols_total, adv_coefs, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_W00_u, data_W10_u, &
                                data_B0_v, data_W00_v, data_W10_v, &
                                data_B0_w, data_W00_w, data_W10_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, & 
                                data_result, indi_result, indj_result)
    !! Computes advention matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: adv_coefs
    dimension :: adv_coefs(3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, data_W10_u, &
                                    data_B0_v, data_W00_v, data_W10_v, &
                                    data_B0_w, data_W00_w, data_W10_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), data_W10_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), data_W10_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w), data_W10_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w 

    double precision, intent(out) :: data_result
    dimension :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1), &
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
    size_data_result =  size_data_I_u*size_data_I_v*size_data_I_w

    ! Get values
    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B0_u (Kronecker product)
    ! Get W = W00_w x W00_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(1, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B0_w, data_W10_u, data_W00_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W10_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(2, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B0_w, data_W00_u, data_W10_v, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result_temp)
    !$OMP PARALLEL
    !$OMP WORKSHARE 
    data_result = data_result + data_result_temp
    !$OMP END WORKSHARE NOWAIT
    !$OMP END PARALLEL 

    ! Get W = W10_w x W00_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(3, :), nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B0_u, data_B0_v, data_B0_w, data_W00_u, data_W00_v, data_W10_w, &
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

end subroutine wq_get_advention_3d

subroutine wq_get_stiffness_3d(nb_cols_total, stiff_coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
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
    !! Computes stiffness matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: stiff_coefs
    dimension :: stiff_coefs(9, 9, nb_cols_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
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
    integer, intent(in) ::  size_data_I_u, size_data_I_v, size_data_I_w 

    integer, parameter :: dimen = 3
    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen*dimen)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen*dimen), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen*dimen)

    ! Local data 
    !-------------
    integer :: nb_rows_total, size_data_result
    double precision, allocatable, dimension(:) :: data_result_temp
    integer, allocatable, dimension(:) :: indi_result_temp_csr, indj_result_temp, indi_result_temp_coo
    integer :: i, j, k, l, nnz, count
    integer :: init, fin

    ! Initialize
    nb_rows_total = nb_rows_u*nb_rows_v*nb_rows_w
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    allocate(data_result_temp(size_data_result))
    allocate(indi_result_temp_csr(nb_rows_total+1), &
            indj_result_temp(size_data_result), indi_result_temp_coo(size_data_result))
    
    do i = 1, dimen
        do j = 1, dimen
            
            call wq_get_conductivity_3D( nb_cols_total, &
            stiff_coefs((i-1)*dimen+1:i*dimen, (j-1)*dimen+1:j*dimen, :), &
            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
            size_data_u, size_data_v, size_data_w, &
            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B0_u, data_B1_u, data_W00_u, data_W01_u, data_W10_u, data_W11_u, &
            data_B0_v, data_B1_v, data_W00_v, data_W01_v, data_W10_v, data_W11_v, &
            data_B0_w, data_B1_w, data_W00_w, data_W01_w, data_W10_w, data_W11_w, &
            size_data_I_u, size_data_I_v, size_data_I_w, &
            data_result_temp, indi_result_temp_csr, indj_result_temp)

            if ((i.eq.1).and.(j.eq.1)) then
                count = 1
                indi_result_temp_coo = 0
                do k = 1, nb_rows_total
                    nnz = indi_result_temp_csr(k+1) - indi_result_temp_csr(k)
                    do l = 1, nnz
                        indi_result_temp_coo(count) = k-1
                        count = count + 1
                    end do
                end do
            end if

            init = (j + (i-1)*dimen - 1)*size_data_result + 1
            fin = (j + (i-1)*dimen)*size_data_result

            data_result(init : fin) = data_result_temp    
            indi_result(init : fin) = indi_result_temp_coo + (i-1)*nb_rows_total
            indj_result(init : fin) = indj_result_temp + (j-1)*nb_rows_total

        end do 
    end do

    deallocate(data_result_temp, indi_result_temp_csr, indi_result_temp_coo, indj_result_temp)

end subroutine wq_get_stiffness_3d

subroutine wq_get_thermalstiffness_3d(nb_cols_total, therstiff_coefs, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, &
                                indi_v, indj_v, &
                                indi_w, indj_w, &
                                data_B0_u, data_W00_u, data_W10_u, &
                                data_B0_v, data_W00_v, data_W10_v, &
                                data_B0_w, data_W00_w, data_W10_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, & 
                                data_result, indi_result, indj_result)
    !! Computes pseudo thermal-stiffness matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_cols_total
    integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
    double precision, intent(in) :: therstiff_coefs
    dimension :: therstiff_coefs(9, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v), &
                    indi_w(nb_rows_w+1), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, data_W10_u, &
                                    data_B0_v, data_W00_v, data_W10_v, &
                                    data_B0_w, data_W00_w, data_W10_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), data_W10_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), data_W10_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w), data_W10_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    integer, parameter :: dimen = 3
    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen)

    ! Local data 
    !-------------
    integer :: nb_rows_total, size_data_result
    double precision, allocatable, dimension(:) :: data_result_temp
    integer, allocatable, dimension(:) :: indi_result_temp_csr, indj_result_temp, indi_result_temp_coo
    integer :: i, k, l, nnz, count
    integer :: init, fin

    ! Initialize
    nb_rows_total = nb_rows_u*nb_rows_v*nb_rows_w
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    allocate(data_result_temp(size_data_result))
    allocate(indi_result_temp_csr(nb_rows_total+1), &
            indj_result_temp(size_data_result), indi_result_temp_coo(size_data_result))

    do i = 1, dimen
        call wq_get_advention_3D(nb_cols_total, therstiff_coefs((i-1)*dimen+1:i*dimen, :), &
                            nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_W00_u, data_W10_u, &
                            data_B0_v, data_W00_v, data_W10_v, &
                            data_B0_w, data_W00_w, data_W10_w, &
                            size_data_I_u, size_data_I_v, size_data_I_w, &
                            data_result_temp, indi_result_temp_csr, indj_result_temp)

        count = 1
        indi_result_temp_coo = 0
        do k = 1, nb_rows_total
            nnz = indi_result_temp_csr(k+1) - indi_result_temp_csr(k)
            do l = 1, nnz
                indi_result_temp_coo(count) = k-1
                count = count + 1
            end do
        end do

        init = (i-1)*size(data_result_temp) + 1
        fin = i*size(data_result_temp)

        data_result(init : fin) = data_result_temp    
        indi_result(init : fin) = indi_result_temp_coo + (i-1)*nb_rows_total
        indj_result(init : fin) = indj_result_temp 

    end do

    deallocate(data_result_temp, indi_result_temp_csr, indi_result_temp_coo, indj_result_temp)

end subroutine wq_get_thermalstiffness_3d

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

subroutine wq_get_flux_3d(nb_cols_total, u, v, coefs, jacob, &
                        nb_rows_u, nb_cols_u, &
                        nb_rows_v, nb_cols_v, &
                        size_data_u, size_data_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_W00_u, data_W00_v, &
                        flux_vector)
    !! Computes flux vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_cols_total, u, v
    integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
    double precision, intent(in) :: coefs, jacob
    dimension :: coefs(nb_cols_total), jacob(3, 3, nb_cols_total)
    integer, intent(in) :: size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                    indi_v(nb_rows_v+1), indj_v(size_data_v)
    double precision, intent(in) :: data_W00_u, data_W00_v
    dimension :: data_W00_u(size_data_u), data_W00_v(size_data_v)

    double precision, intent(out) :: flux_vector
    dimension :: flux_vector(nb_rows_u*nb_rows_v)

    ! Local data
    ! ----------------------
    integer :: i
    double precision :: cross_vector(3), coefs_final(nb_cols_total)

    ! Evaluate determinant of surface transformation
    do i = 1, nb_cols_total
        ! Eval cross product 
        call crossproduct(jacob(:, u, i), jacob(:, v, i), cross_vector)

        ! Eval final coefficient
        coefs_final(i) = norm2(cross_vector)*coefs(i)

    end do

    ! Find vector
    call tensor2d_dot_vector_sp(nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                        size_data_u, indi_u, indj_u, data_W00_u, &
                        size_data_v, indi_v, indj_v, data_W00_v, &
                        coefs_final, flux_vector)

end subroutine wq_get_flux_3d

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