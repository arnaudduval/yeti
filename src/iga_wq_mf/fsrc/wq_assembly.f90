! ==========================
! module :: assembly for IGA-WQ 
! author :: Joaquin Cornejo
! modules :: operateurs.f90 (MatrixInv and MatrixDet)
! ==========================

subroutine eval_thermal_coefficient(dime, nb_pts, J_pts, Kprop, Cprop, Kcoef, Ccoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------    
    integer, intent(in) :: dime, nb_pts
    double precision, intent(in) :: J_pts, Kprop
    dimension :: J_pts(dime, dime, nb_pts), Kprop(dime, dime)
    double precision, intent(in) :: Cprop

    double precision, intent(out) :: Kcoef, Ccoef
    dimension :: Kcoef(dime, dime, nb_pts), Ccoef(nb_pts)

    ! Local data
    ! -----------
    integer :: i, nb_tasks
    double precision :: J, detJ, invJ
    dimension :: J(dime, dime), invJ(dime, dime)
    double precision :: Ktemp
    dimension :: Ktemp(dime, dime)  

    !$OMP PARALLEL PRIVATE(J, invJ, detJ, Ktemp)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_pts/nb_tasks) 
    do i = 1, nb_pts
        ! Get individual jacobien
        J = J_pts(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJ, J, detJ, dime)

        ! For K = invJ * detJ * transpose(invJ) * prop 
        Ktemp = matmul(invJ, transpose(invJ)) * detJ
        Kcoef(:, :, i) = matmul(Ktemp, Kprop)

        ! For C = detJ  * prop
        Ccoef(i) = detJ * Cprop
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_thermal_coefficient

subroutine eval_mech_coefficient(dime, nb_pts, J_pts, DD, Scoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in):: dime, nb_pts
    double precision, intent(in) :: J_pts, DD
    dimension ::    J_pts(dime, dime, nb_pts), &
                    DD(3*(dime-1), 3*(dime-1))
    double precision, intent(out) :: Scoef
    dimension ::    Scoef(dime*dime, dime*dime, nb_pts)

    ! Local data
    ! -----------
    external :: MatrixInv
    integer :: i, k, nb_tasks
    double precision :: J, detJ, invJ
    dimension :: J(dime, dime), &
                invJ(dime, dime)
    double precision :: Stemp1, Stemp2
    dimension :: Stemp1(dime*dime, dime*dime), &
                    Stemp2(dime*dime, dime*dime)  

    double precision, allocatable, dimension(:, :) :: MM
    double precision, allocatable, dimension(:, :) :: MDM_temp
    double precision :: MDM(dime*dime, dime*dime)
    double precision :: invJextend(dime*dime, dime*dime)

    ! Define MM transformation matrix 
    if (dime.eq.2) then 
        allocate(MM(3, 4))
        allocate(MDM_temp(4, 3))
        MM = 0.0d0
        MM(1, 1) = 1.0d0
        MM(2, 4) = 1.0d0
        MM(3, 2) = 1.0d0
        MM(3, 3) = 1.0d0

    else if (dime.eq.3) then 
        allocate(MM(6, 9))
        allocate(MDM_temp(9, 6))
        MM = 0.0d0
        MM(1, 1) = 1.0d0
        MM(2, 5) = 1.0d0
        MM(3, 9) = 1.0d0
        MM(4, 2) = 1.0d0
        MM(4, 4) = 1.0d0
        MM(5, 6) = 1.0d0
        MM(5, 8) = 1.0d0
        MM(6, 3) = 1.0d0
        MM(6, 7) = 1.0d0
    end if

    ! Evaluate MM.transpose * DD * MM
    MDM_temp = matmul(transpose(MM), DD)
    MDM = matmul(MDM_temp, MM)
    deallocate(MM, MDM_temp)

    !$OMP PARALLEL PRIVATE(J, invJ, detJ, invJextend, Stemp1, Stemp2, k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_pts/nb_tasks) 
    do i = 1, nb_pts

        ! Initialize
        invJextend = 0.0d0

        ! Get individual jacobien
        J = J_pts(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJ, J, detJ, dime)

        do k = 1, dime
            invJextend((k-1)*dime+1:k*dime, (k-1)*dime+1:k*dime) = invJ
        end do
        
        Stemp1 = matmul(transpose(invJextend), MDM)
        Stemp2 = matmul(Stemp1, invJextend) * detJ

        Scoef(:, :, i) = Stemp2
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_mech_coefficient

subroutine eval_thermomech_coefficient(dime, nb_pts, J_pts, DD, alpha, Tcoef)
    !! Computes coefficient for K, C matrices and F vector
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dime, nb_pts
    double precision, intent(in) :: J_pts, DD
    dimension ::    J_pts(dime, dime, nb_pts), &
                    DD(3*(dime-1), 3*(dime-1))
    double precision, intent(in) :: alpha
            
    double precision, intent(out) :: Tcoef
    dimension ::    Tcoef(dime*dime, nb_pts)

    ! Local data
    ! -----------
    external :: MatrixInv
    integer :: i, k, nb_tasks
    double precision :: J, detJ, invJ
    dimension :: J(dime, dime), &
                invJ(dime, dime)
    double precision :: Rtemp1
    dimension :: Rtemp1(dime*dime)

    double precision, allocatable, dimension(:, :) :: MM
    double precision :: MDV(dime*dime), VV(3*(dime-1)), DV(3*(dime-1))
    double precision :: invJextend(dime*dime, dime*dime)

    ! Define MM transformation matrix 
    if (dime.eq.2) then 
        allocate(MM(3, 4))
        MM = 0.0d0
        MM(1, 1) = 1.0d0
        MM(2, 4) = 1.0d0
        MM(3, 2) = 1.0d0
        MM(3, 3) = 1.0d0

        VV = 0.0d0
        VV(1) = 1.0d0
        VV(2) = 1.0d0

    else if (dime.eq.3) then 
        allocate(MM(6, 9))
        MM = 0.0d0
        MM(1, 1) = 1.0d0
        MM(2, 5) = 1.0d0
        MM(3, 9) = 1.0d0
        MM(4, 2) = 1.0d0
        MM(4, 4) = 1.0d0
        MM(5, 6) = 1.0d0
        MM(5, 8) = 1.0d0
        MM(6, 3) = 1.0d0
        MM(6, 7) = 1.0d0

        VV = 0.0d0
        VV(1) = 1.0d0
        VV(2) = 1.0d0
        VV(3) = 1.0d0
    end if

    ! Evaluate MM.transpose * DD * MM
    DV = matmul(DD, VV)
    MDV = matmul(transpose(MM), DV)
    deallocate(MM)

    !$OMP PARALLEL PRIVATE(J, invJ, detJ, invJextend, Rtemp1, k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_pts/nb_tasks) 
    do i = 1, nb_pts

        ! Initialize
        invJextend = 0.0d0

        ! Get individual jacobien
        J = J_pts(:, :, i)

        ! Evaluate inverse
        call MatrixInv(invJ, J, detJ, dime)

        do k = 1, dime
            invJextend((k-1)*dime+1:k*dime, (k-1)*dime+1:k*dime) = invJ
        end do
        
        Rtemp1 = matmul(transpose(invJextend), MDV) * detJ * alpha
        Tcoef(:, i) = Rtemp1
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_thermomech_coefficient

subroutine jacobien_physicalposition_3d(nb_ctrlpts_total, &
                                        nb_ctrlpts_u, nb_qp_u, &
                                        nb_ctrlpts_v, nb_qp_v, &
                                        nb_ctrlpts_w, nb_qp_w, &
                                        size_data_u, size_data_v, size_data_w, &
                                        ctrlpts_x, ctrlpts_y, ctrlpts_z, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B0_u, data_B1_u, &
                                        data_B0_v, data_B1_v, &
                                        data_B0_w, data_B1_w, &
                                        jacob, physical_pos, detJ)
    !! Computes jacobien in 3D case
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) :: nb_ctrlpts_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w, &
                nb_qp_u, nb_qp_v, nb_qp_w, &
                size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_B0_v, data_B1_v, &
                                    data_B0_w, data_B1_w
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v), &
                    data_B0_w(size_data_w), data_B1_w(size_data_w)
    double precision, intent(in) :: ctrlpts_x, ctrlpts_y, ctrlpts_z
    dimension ::    ctrlpts_x(nb_ctrlpts_total), &
                    ctrlpts_y(nb_ctrlpts_total), &
                    ctrlpts_z(nb_ctrlpts_total)

    double precision, intent(out) :: jacob
    dimension ::  jacob(3, 3, nb_qp_u*nb_qp_v*nb_qp_w)

    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(1, 3, nb_qp_u*nb_qp_v*nb_qp_w)

    double precision, intent(out) :: detJ
    dimension :: detJ(nb_qp_u*nb_qp_v*nb_qp_w)

    ! Local data
    !-----------------
    double precision :: result_temp
    dimension ::  result_temp(nb_qp_u*nb_qp_v*nb_qp_w)
    integer :: nb_tasks, i
    double precision :: detJ_temp

    ! Csr format (Transpose)
    integer ::  indi_T_u_csr, indi_T_v_csr, indi_T_w_csr
    dimension ::    indi_T_u_csr(nb_qp_u+1), &
                    indi_T_v_csr(nb_qp_v+1), &
                    indi_T_w_csr(nb_qp_w+1)
    integer ::  indj_T_u_csr, indj_T_v_csr, indj_T_w_csr
    dimension ::    indj_T_u_csr(size_data_u), &
                    indj_T_v_csr(size_data_v), &
                    indj_T_w_csr(size_data_w)
    double precision :: data_B0T_u_csr, data_B0T_v_csr, data_B0T_w_csr, &
                        data_B1T_u_csr, data_B1T_v_csr, data_B1T_w_csr
    dimension ::    data_B0T_u_csr(size_data_u), &
                    data_B0T_v_csr(size_data_v), &
                    data_B0T_w_csr(size_data_w), &
                    data_B1T_u_csr(size_data_u), &
                    data_B1T_v_csr(size_data_v), &
                    data_B1T_w_csr(size_data_w)
    ! ====================================================
    ! Initialize
    call coo2csr(nb_qp_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    call coo2csr(nb_qp_w, size_data_w, data_B0_w, indj_w, indi_w, data_B0T_w_csr, &
                    indj_T_w_csr, indi_T_w_csr)
    call coo2csr(nb_qp_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    call coo2csr(nb_qp_w, size_data_w, data_B1_w, indj_w, indi_w, data_B1T_w_csr, &
                    indj_T_w_csr, indi_T_w_csr)
    ! ====================================================
    ! ---------------------------------------------------
    ! For J00, J10 and J20
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B1T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_x, result_temp)

    jacob(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B1T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_y, result_temp)
    jacob(2, 1, :) = result_temp

    ! Compute B.Transpose . CP_z
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B1T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_z, result_temp)
    jacob(3, 1, :) = result_temp

    ! ---------------------------------------------------
    ! For J01, J11, and J21
    ! ---------------------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)

    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B1T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_x, result_temp)
    jacob(1, 2, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B1T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_y, result_temp)
    jacob(2, 2, :) = result_temp

    ! Compute B.Transpose . CP_z
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B1T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_z, result_temp)
    jacob(3, 2, :) = result_temp

    ! ---------------------------------------------------
    ! For J02, J12, and J22
    ! ---------------------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B1T_w_csr, &
                                ctrlpts_x, result_temp)
    jacob(1, 3, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B1T_w_csr, &
                                ctrlpts_y, result_temp)
    jacob(2, 3, :) = result_temp

    ! Compute B.Transpose . CP_z
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B1T_w_csr, &
                                ctrlpts_z, result_temp)
    jacob(3, 3, :) = result_temp

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B0_w x B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_x, result_temp)
    physical_pos(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_y, result_temp)
    physical_pos(1, 2, :) = result_temp

    ! Compute B.Transpose . CP_z
    result_temp = 0.0d0
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_B0T_w_csr, &
                                ctrlpts_z, result_temp)
    physical_pos(1, 3, :) = result_temp

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJ_temp)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_qp_u*nb_qp_v*nb_qp_w/nb_tasks) 
    do i = 1, nb_qp_u*nb_qp_v*nb_qp_w
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, i), detJ_temp, 3)

        ! Assign values
        detJ(i) = detJ_temp
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_3d

subroutine interpolation_3d(nb_ctrlpts_total, &
                            nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, &
                            nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            ctrlpts, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, &
                            interpolation)

    !! Computes interpolation in 3D case (from parametric space to physical space)

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in)  :: nb_ctrlpts_total
    integer, intent(in)  ::  nb_ctrlpts_u, nb_ctrlpts_v, nb_ctrlpts_w, &
                nb_qp_u, nb_qp_v, nb_qp_w, &
                size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension ::    data_B_u(size_data_u), &
                    data_B_v(size_data_v), &
                    data_B_w(size_data_w)
    double precision, intent(in) :: ctrlpts
    dimension ::    ctrlpts(nb_ctrlpts_total)

    double precision, intent(out) :: interpolation
    dimension :: interpolation(nb_qp_u*nb_qp_v*nb_qp_w)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer ::  indi_T_u_csr, indi_T_v_csr, indi_T_w_csr
    dimension ::    indi_T_u_csr(nb_qp_u+1), &
                    indi_T_v_csr(nb_qp_v+1), &
                    indi_T_w_csr(nb_qp_w+1)
    integer ::  indj_T_u_csr, indj_T_v_csr, indj_T_w_csr
    dimension ::    indj_T_u_csr(size_data_u), &
                    indj_T_v_csr(size_data_v), &
                    indj_T_w_csr(size_data_w)
    double precision :: data_BT_u_csr, data_BT_v_csr, data_BT_w_csr
    dimension ::    data_BT_u_csr(size_data_u), &
                    data_BT_v_csr(size_data_v), &
                    data_BT_w_csr(size_data_w)

    ! ====================================================
    ! Initialize
    call coo2csr(nb_qp_u, size_data_u, data_B_u, indj_u, indi_u, data_BT_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B_v, indj_v, indi_v, data_BT_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    call coo2csr(nb_qp_w, size_data_w, data_B_w, indj_w, indi_w, data_BT_w_csr, &
                    indj_T_w_csr, indi_T_w_csr)
    ! ====================================================

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B_w x B_v x B_u (Kronecker product)
    interpolation = 0.0d0
    
    ! Compute B.Transpose . CP
    call tensor3d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, nb_qp_w, nb_ctrlpts_w, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_BT_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_BT_v_csr, &
                                size_data_w, indi_T_w_csr, indj_T_w_csr, data_BT_w_csr, &
                                ctrlpts, interpolation)

end subroutine interpolation_3d

subroutine jacobien_physicalposition_2d(nb_ctrlpts_total, &
                                        nb_ctrlpts_u, nb_qp_u, &
                                        nb_ctrlpts_v, nb_qp_v, &
                                        size_data_u, size_data_v, &
                                        ctrlpts_x, ctrlpts_y, &
                                        indi_u, indj_u, indi_v, indj_v, &
                                        data_B0_u, data_B1_u, &
                                        data_B0_v, data_B1_v, &
                                        jacob, physical_pos, detJ)
    !! Computes jacobien in 3D case
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) :: nb_ctrlpts_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_ctrlpts_v, &
                nb_qp_u, nb_qp_v, &
                size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v)
    double precision, intent(in) :: data_B0_u, data_B1_u, &
                                    data_B0_v, data_B1_v
    dimension ::    data_B0_u(size_data_u), data_B1_u(size_data_u), &
                    data_B0_v(size_data_v), data_B1_v(size_data_v)
    double precision, intent(in) :: ctrlpts_x, ctrlpts_y
    dimension ::    ctrlpts_x(nb_ctrlpts_total), &
                    ctrlpts_y(nb_ctrlpts_total)

    double precision, intent(out) :: jacob
    dimension ::  jacob(2, 2, nb_qp_u*nb_qp_v)

    double precision, intent(out) :: physical_pos
    dimension :: physical_pos(1, 2, nb_qp_u*nb_qp_v)

    double precision, intent(out) :: detJ
    dimension :: detJ(nb_qp_u*nb_qp_v)

    ! Local data
    !-----------------
    double precision :: result_temp
    dimension ::  result_temp(nb_qp_u*nb_qp_v)
    integer :: nb_tasks, i
    double precision :: detJ_temp

    ! Csr format (Transpose)
    integer ::  indi_T_u_csr, indi_T_v_csr
    dimension ::    indi_T_u_csr(nb_qp_u+1), &
                    indi_T_v_csr(nb_qp_v+1)
    integer ::  indj_T_u_csr, indj_T_v_csr
    dimension ::    indj_T_u_csr(size_data_u), &
                    indj_T_v_csr(size_data_v)
    double precision :: data_B0T_u_csr, data_B0T_v_csr,&
                        data_B1T_u_csr, data_B1T_v_csr
    dimension ::    data_B0T_u_csr(size_data_u), &
                    data_B0T_v_csr(size_data_v), &
                    data_B1T_u_csr(size_data_u), &
                    data_B1T_v_csr(size_data_v)
    ! ====================================================
    ! Initialize
    call coo2csr(nb_qp_u, size_data_u, data_B0_u, indj_u, indi_u, data_B0T_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B0_v, indj_v, indi_v, data_B0T_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    call coo2csr(nb_qp_u, size_data_u, data_B1_u, indj_u, indi_u, data_B1T_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B1_v, indj_v, indi_v, data_B1T_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    ! ====================================================
    ! ---------------------------------------------------
    ! For J00, J10 and J20
    ! ---------------------------------------------------
    ! Get B = B0_v x B1_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B1T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                ctrlpts_x, result_temp)

    jacob(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B1T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                ctrlpts_y, result_temp)
    jacob(2, 1, :) = result_temp

    ! ---------------------------------------------------
    ! For J01, J11, and J21
    ! ---------------------------------------------------
    ! Get B = B1_v x B0_u (Kronecker product)

    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B1T_v_csr, &
                                ctrlpts_x, result_temp)
    jacob(1, 2, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B1T_v_csr, &
                                ctrlpts_y, result_temp)
    jacob(2, 2, :) = result_temp

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B0_v x B0_u (Kronecker product)
    
    ! Compute B.Transpose . CP_x
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                ctrlpts_x, result_temp)
    physical_pos(1, 1, :) = result_temp

    ! Compute B.Transpose . CP_y
    result_temp = 0.0d0
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v,&
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_B0T_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_B0T_v_csr, &
                                ctrlpts_y, result_temp)
    physical_pos(1, 2, :) = result_temp

    ! ---------------------------------------------------
    ! For det J 
    ! ---------------------------------------------------
    !$OMP PARALLEL PRIVATE(detJ_temp)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nb_qp_u*nb_qp_v/nb_tasks) 
    do i = 1, nb_qp_u*nb_qp_v
        ! Evaluate determinant
        call MatrixDet(jacob(:, :, i), detJ_temp, 2)

        ! Assign values
        detJ(i) = detJ_temp
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_2d

subroutine interpolation_2d(nb_ctrlpts_total, &
                            nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, &
                            size_data_u, size_data_v, ctrlpts, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            data_B_u, data_B_v, &
                            interpolation)

    !! Computes interpolation in 3D case (from parametric space to physical space)

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------    
    integer, intent(in) :: nb_ctrlpts_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_ctrlpts_v, &
                nb_qp_u, nb_qp_v, &
                size_data_u, size_data_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension ::    data_B_u(size_data_u), &
                    data_B_v(size_data_v)
    double precision, intent(in) :: ctrlpts
    dimension ::    ctrlpts(nb_ctrlpts_total)

    double precision, intent(out) :: interpolation
    dimension :: interpolation(nb_qp_u*nb_qp_v)

    ! Local data
    !-----------------
    ! Csr format (Transpose)
    integer ::  indi_T_u_csr, indi_T_v_csr
    dimension ::    indi_T_u_csr(nb_qp_u+1), &
                    indi_T_v_csr(nb_qp_v+1)
    integer ::  indj_T_u_csr, indj_T_v_csr
    dimension ::    indj_T_u_csr(size_data_u), &
                    indj_T_v_csr(size_data_v)
    double precision :: data_BT_u_csr, data_BT_v_csr
    dimension ::    data_BT_u_csr(size_data_u), &
                    data_BT_v_csr(size_data_v)

    ! ====================================================
    ! Initialize
    call coo2csr(nb_qp_u, size_data_u, data_B_u, indj_u, indi_u, data_BT_u_csr, &
                    indj_T_u_csr, indi_T_u_csr)
    call coo2csr(nb_qp_v, size_data_v, data_B_v, indj_v, indi_v, data_BT_v_csr, &
                    indj_T_v_csr, indi_T_v_csr)
    ! ====================================================

    ! ---------------------------------------------------
    ! For position 
    ! ---------------------------------------------------
    ! Get B = B_v x B_u (Kronecker product)
    interpolation = 0.0d0
    
    ! Compute B.Transpose . CP
    call tensor2d_sparsedot_vector(nb_qp_u, nb_ctrlpts_u, &
                                nb_qp_v, nb_ctrlpts_v, &
                                size_data_u, indi_T_u_csr, indj_T_u_csr, data_BT_u_csr, &
                                size_data_v, indi_T_v_csr, indj_T_v_csr, data_BT_v_csr, &
                                ctrlpts, interpolation)

end subroutine interpolation_2d

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine wq_get_capacity_3d(nb_qp_total, capacity_coefs, &
                            nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                            size_data_I_u, size_data_I_v, size_data_I_w, & 
                            data_result, indi_result, indj_result)

    !! Computes a matrix in 3D case

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: capacity_coefs
    dimension :: capacity_coefs(nb_qp_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, &
                                    data_B0_v, data_W00_v, &
                                    data_B0_w, data_W00_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w+1), &
                indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer ::  size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension ::    indi_I_u(nb_ctrlpts_u+1), &
                    indi_I_v(nb_ctrlpts_v+1), &
                    indi_I_w(nb_ctrlpts_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w

    integer :: nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp
    integer :: dummy1, dummy2, dummy3
    integer :: indi_I_temp(nb_ctrlpts_u*nb_ctrlpts_v+1)
    integer, allocatable, dimension(:) :: indj_I_temp

    integer :: indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr
    dimension :: indi_u_csr(nb_ctrlpts_u+1), indj_u_csr(size_data_u), &
                indi_v_csr(nb_ctrlpts_v+1), indj_v_csr(size_data_v), &
                indi_w_csr(nb_ctrlpts_w+1), indj_w_csr(size_data_w)
    double precision, dimension(:), allocatable :: data_dump_csr

    ! Get indexes of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_ctrlpts_u, nb_qp_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_ctrlpts_v, nb_qp_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_ctrlpts_w, nb_qp_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    allocate(indj_I_temp(size_data_I_u*size_data_I_v))
    call get_indexes_kron_product(nb_ctrlpts_v, nb_ctrlpts_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_ctrlpts_u, nb_ctrlpts_u, size_data_I_u, indi_I_u, indj_I_u, &
                                nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp)

    call get_indexes_kron_product(nb_ctrlpts_w, nb_ctrlpts_w, size_data_I_w, indi_I_w, indj_I_w, &
                                nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    deallocate(indj_I_temp)

    ! Get Indexes of B and W in each dimension
    allocate(data_dump_csr(size_data_u))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_u, size_data_u, data_B0_u, indi_u, indj_u, data_dump_csr, indj_u_csr, indi_u_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_v))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_v, size_data_v, data_B0_v, indi_v, indj_v, data_dump_csr, indj_v_csr, indi_v_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_w))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_w, size_data_w, data_B0_w, indi_w, indj_w, data_dump_csr, indj_w_csr, indi_w_csr)
    deallocate(data_dump_csr)

    ! Initialize 
    data_result = 0.0d0
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    call csr_get_matrix_3d(capacity_coefs, nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                            data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                            size_data_I_u, size_data_I_v, size_data_I_w, & 
                            indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, size_data_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_3d

subroutine wq_get_conductivity_3d(nb_qp_total, cond_coefs, &
                                nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
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

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: cond_coefs
    dimension :: cond_coefs(3, 3, nb_qp_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
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
    dimension :: indi_result(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w+1), &
                indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension ::    indi_I_u(nb_ctrlpts_u+1), &
                    indi_I_v(nb_ctrlpts_v+1), &
                    indi_I_w(nb_ctrlpts_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w

    integer :: nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp 
    integer :: dummy1, dummy2, dummy3
    integer :: indi_I_temp(nb_ctrlpts_u*nb_ctrlpts_v+1)
    integer, allocatable, dimension(:) :: indj_I_temp

    integer :: indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr
    dimension :: indi_u_csr(nb_ctrlpts_u+1), indj_u_csr(size_data_u), &
                indi_v_csr(nb_ctrlpts_v+1), indj_v_csr(size_data_v), &
                indi_w_csr(nb_ctrlpts_w+1), indj_w_csr(size_data_w)
    double precision, dimension(:), allocatable :: data_dump_csr

    ! Get indexes of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_ctrlpts_u, nb_qp_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_ctrlpts_v, nb_qp_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_ctrlpts_w, nb_qp_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    allocate(indj_I_temp(size_data_I_u*size_data_I_v))
    call get_indexes_kron_product(nb_ctrlpts_v, nb_ctrlpts_v, size_data_I_v, indi_I_v, indj_I_v, &
                                nb_ctrlpts_u, nb_ctrlpts_u, size_data_I_u, indi_I_u, indj_I_u, &
                                nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp)

    call get_indexes_kron_product(nb_ctrlpts_w, nb_ctrlpts_w, size_data_I_w, indi_I_w, indj_I_w, &
                                nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    deallocate(indj_I_temp)

    ! Get Indexes of B and W in each dimension
    allocate(data_dump_csr(size_data_u))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_u, size_data_u, data_B0_u, indi_u, indj_u, data_dump_csr, indj_u_csr, indi_u_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_v))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_v, size_data_v, data_B0_v, indi_v, indj_v, data_dump_csr, indj_v_csr, indi_v_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_w))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_w, size_data_w, data_B0_w, indi_w, indj_w, data_dump_csr, indj_w_csr, indi_w_csr)
    deallocate(data_dump_csr)

    ! Initialize 
    data_result = 0.0d0
    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    
    ! Get values
    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B1_u (Kronecker product)
    ! Get W = W00_w x W00_v x W11_u (Kronecker produt)

    call csr_get_matrix_3d(cond_coefs(1, 1, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B1_u, data_W11_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W10_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 1, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B1_u, data_W01_u, data_B0_v, data_W10_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_w x W00_v x W01_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 1, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B1_u, data_W01_u, data_B0_v, data_W00_v, data_B0_w, data_W10_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! ----------------------------------------
    ! For c01, c11 and c21
    ! ----------------------------------------
    ! Get B = B0_w x B1_v x B0_u (Kronecker product)
    ! Get W = W00_w x W01_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 2, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W10_u, data_B1_v, data_W01_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W11_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 2, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B1_v, data_W11_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_w x W01_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 2, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B1_v, data_W01_v, data_B0_w, data_W10_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! ----------------------------------------
    ! For c02, c12 and c22
    ! ----------------------------------------
    ! Get B = B1_w x B0_v x B0_u (Kronecker product)
    ! Get W = W01_w x W00_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(1, 3, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W10_u, data_B0_v, data_W00_v, data_B1_w, data_W01_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W01_w x W10_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(2, 3, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B0_v, data_W10_v, data_B1_w, data_W01_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W11_w x W00_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(cond_coefs(3, 3, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B1_w, data_W11_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, &
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    deallocate(indj_I_u, indj_I_v, indj_I_w)
    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_conductivity_3d

subroutine wq_get_advention_3d(nb_qp_total, adv_coefs, &
                                nb_ctrlpts_u, nb_qp_u, &
                                nb_ctrlpts_v, nb_qp_v, &
                                nb_ctrlpts_w, nb_qp_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B0_u, data_W00_u, data_W10_u, &
                                data_B0_v, data_W00_v, data_W10_v, &
                                data_B0_w, data_W00_w, data_W10_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, & 
                                data_result, indi_result, indj_result)
    !! Computes advention matrix in 3D case

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in)  :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: adv_coefs
    dimension :: adv_coefs(3, nb_qp_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, data_W10_u, &
                                    data_B0_v, data_W00_v, data_W10_v, &
                                    data_B0_w, data_W00_w, data_W10_w
    dimension ::    data_B0_u(size_data_u), data_W00_u(size_data_u), &
                    data_W10_u(size_data_u), &
                    data_B0_v(size_data_v), data_W00_v(size_data_v), &
                    data_W10_v(size_data_v), &
                    data_B0_w(size_data_w), data_W00_w(size_data_w), &
                    data_W10_w(size_data_w)
    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w 

    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w+1), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w)

    ! Local data
    ! ---------------
    integer :: size_data_result
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension ::    indi_I_u(nb_ctrlpts_u+1), &
                    indi_I_v(nb_ctrlpts_v+1), &
                    indi_I_w(nb_ctrlpts_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w

    integer :: nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp 
    integer :: dummy1, dummy2, dummy3
    integer :: indi_I_temp(nb_ctrlpts_u*nb_ctrlpts_v+1)
    integer, allocatable, dimension(:) :: indj_I_temp

    integer :: indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr
    dimension :: indi_u_csr(nb_ctrlpts_u+1), indj_u_csr(size_data_u), &
                indi_v_csr(nb_ctrlpts_v+1), indj_v_csr(size_data_v), &
                indi_w_csr(nb_ctrlpts_w+1), indj_w_csr(size_data_w)
    double precision, dimension(:), allocatable :: data_dump_csr

    ! Get indexes of I in each dimension
    allocate(indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w))
    call get_I_csr(nb_ctrlpts_u, nb_qp_u, size_data_u, data_B0_u, indi_u, indj_u, size_data_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nb_ctrlpts_v, nb_qp_v, size_data_v, data_B0_v, indi_v, indj_v, size_data_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nb_ctrlpts_w, nb_qp_w, size_data_w, data_B0_w, indi_w, indj_w, size_data_I_w, indi_I_w, indj_I_w)

    allocate(indj_I_temp(size_data_I_u*size_data_I_v))
    call get_indexes_kron_product(nb_ctrlpts_v, nb_ctrlpts_v, size_data_I_v, indi_I_v, indj_I_v, &
        nb_ctrlpts_u, nb_ctrlpts_u, size_data_I_u, indi_I_u, indj_I_u, &
        nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp)

    call get_indexes_kron_product(nb_ctrlpts_w, nb_ctrlpts_w, size_data_I_w, indi_I_w, indj_I_w, &
        nb_ctrlpts_temp1, nb_ctrlpts_temp2, size_data_I_temp, indi_I_temp, indj_I_temp, &
        dummy1, dummy2, dummy3, indi_result, indj_result)

    deallocate(indj_I_temp)

    ! Get Indexes of B and W in each dimension
    allocate(data_dump_csr(size_data_u))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_u, size_data_u, data_B0_u, indi_u, indj_u, data_dump_csr, indj_u_csr, indi_u_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_v))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_v, size_data_v, data_B0_v, indi_v, indj_v, data_dump_csr, indj_v_csr, indi_v_csr)
    deallocate(data_dump_csr)

    allocate(data_dump_csr(size_data_w))
    data_dump_csr = 0.0d0
    call coo2csr(nb_ctrlpts_w, size_data_w, data_B0_w, indi_w, indj_w, data_dump_csr, indj_w_csr, indi_w_csr)
    deallocate(data_dump_csr)

    ! Initialize 
    data_result = 0.0d0
    size_data_result =  size_data_I_u*size_data_I_v*size_data_I_w

    ! Get values
    ! ----------------------------------------
    ! For c00, c10 and c20
    ! ----------------------------------------
    ! Get B = B0_w x B0_v x B0_u (Kronecker product)
    ! Get W = W00_w x W00_v x W10_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(1, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W10_u, data_B0_v, data_W00_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W00_w x W10_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(2, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B0_v, data_W10_v, data_B0_w, data_W00_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    ! Get W = W10_w x W00_v x W00_u (Kronecker produt)
    call csr_get_matrix_3d(adv_coefs(3, :), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u_csr, indj_u_csr, indi_v_csr, indj_v_csr, indi_w_csr, indj_w_csr, &
                        data_B0_u, data_W00_u, data_B0_v, data_W00_v, data_B0_w, data_W10_w, &
                        size_data_I_u, size_data_I_v, size_data_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size_data_result, data_result)

    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_advention_3d

subroutine wq_get_stiffness_3d(nb_qp_total, stiff_coefs, &
                            nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, &
                            nb_ctrlpts_w, nb_qp_w, &
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

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: stiff_coefs
    dimension :: stiff_coefs(9, 9, nb_qp_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(size_data_u), indj_u(size_data_u), &
                    indi_v(size_data_v), indj_v(size_data_v), &
                    indi_w(size_data_w), indj_w(size_data_w)
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
    integer :: size_data_result
    double precision, allocatable, dimension(:) :: data_result_temp
    integer, allocatable, dimension(:) :: indi_result_temp_csr, indj_result_temp, indi_result_temp_coo
    integer :: i, j, k, l, nnz, count
    integer :: init, fin

    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    allocate(data_result_temp(size_data_result))
    allocate(indi_result_temp_csr(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w+1), &
            indj_result_temp(size_data_result), indi_result_temp_coo(size_data_result))
    
    do i = 1, dimen
        do j = 1, dimen
            
            call wq_get_conductivity_3D( nb_qp_total, &
            stiff_coefs((i-1)*dimen+1:i*dimen, (j-1)*dimen+1:j*dimen, :), &
            nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
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
                do k = 1, nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w
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
            indi_result(init : fin) = indi_result_temp_coo + (i-1)*nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w
            indj_result(init : fin) = indj_result_temp + (j-1)*nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w

        end do 
    end do

    deallocate(data_result_temp, indi_result_temp_csr, indi_result_temp_coo, indj_result_temp)

end subroutine wq_get_stiffness_3d

subroutine wq_get_thermalstiffness_3d(nb_qp_total, therstiff_coefs, &
                                nb_ctrlpts_u, nb_qp_u, &
                                nb_ctrlpts_v, nb_qp_v, &
                                nb_ctrlpts_w, nb_qp_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, &
                                indi_v, indj_v, &
                                indi_w, indj_w, &
                                data_B0_u, data_W00_u, data_W10_u, &
                                data_B0_v, data_W00_v, data_W10_v, &
                                data_B0_w, data_W00_w, data_W10_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, & 
                                data_result, indi_result, indj_result)
    !! Computes pseudo thermal-stiffness matrix in 2D case

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: therstiff_coefs
    dimension :: therstiff_coefs(9, nb_qp_total)
    integer, intent(in) ::  size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w
    dimension ::    indi_u(size_data_u), &
                    indj_u(size_data_u), &
                    indi_v(size_data_v), &
                    indj_v(size_data_v), &
                    indi_w(size_data_w), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_B0_u, data_W00_u, data_W10_u, &
                                    data_B0_v, data_W00_v, data_W10_v, &
                                    data_B0_w, data_W00_w, data_W10_w
    dimension ::    data_B0_u(size_data_u), &
                    data_W00_u(size_data_u), &
                    data_W10_u(size_data_u), &
                    data_B0_v(size_data_v), &
                    data_W00_v(size_data_v), &
                    data_W10_v(size_data_v), &
                    data_B0_w(size_data_w), &
                    data_W00_w(size_data_w), &
                    data_W10_w(size_data_w)

    integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w

    integer, parameter :: dimen = 3
    double precision, intent(out) :: data_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen), &
                    indj_result(size_data_I_u*size_data_I_v*size_data_I_w*dimen)

    ! Local data 
    !-------------
    integer :: size_data_result
    double precision, allocatable, dimension(:) :: data_result_temp
    integer, allocatable, dimension(:) :: indi_result_temp_csr, indj_result_temp, indi_result_temp_coo
    integer :: i, k, l, nnz, count
    integer :: init, fin

    size_data_result = size_data_I_u*size_data_I_v*size_data_I_w
    allocate(data_result_temp(size_data_result))
    allocate(indi_result_temp_csr(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w+1), &
            indj_result_temp(size_data_result), indi_result_temp_coo(size_data_result))

    do i = 1, dimen
        call wq_get_advention_3D(nb_qp_total, therstiff_coefs((i-1)*dimen+1:i*dimen, :), &
                            nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B0_u, data_W00_u, data_W10_u, &
                            data_B0_v, data_W00_v, data_W10_v, &
                            data_B0_w, data_W00_w, data_W10_w, &
                            size_data_I_u, size_data_I_v, size_data_I_w, &
                            data_result_temp, indi_result_temp_csr, indj_result_temp)

        count = 1
        indi_result_temp_coo = 0
        do k = 1, nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w
            nnz = indi_result_temp_csr(k+1) - indi_result_temp_csr(k)
            do l = 1, nnz
                indi_result_temp_coo(count) = k-1
                count = count + 1
            end do
        end do

        init = (i-1)*size(data_result_temp) + 1
        fin = i*size(data_result_temp)

        data_result(init : fin) = data_result_temp    
        indi_result(init : fin) = indi_result_temp_coo + (i-1)*nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w
        indj_result(init : fin) = indj_result_temp 

    end do

    deallocate(data_result_temp, indi_result_temp_csr, indi_result_temp_coo, indj_result_temp)

end subroutine wq_get_thermalstiffness_3d

subroutine wq_get_source_3d(nb_qp_total, source_coefs, &
                            nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, &
                            nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W00_u, data_W00_v, data_W00_w, &
                            source_vector)
    !! Computes source vector in 3D case

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: source_coefs
    dimension :: source_coefs(nb_qp_total)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w
    dimension ::    indi_u(size_data_u), &
                    indj_u(size_data_u), &
                    indi_v(size_data_v), &
                    indj_v(size_data_v), &
                    indi_w(size_data_w), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W00_v, data_W00_w
    dimension ::    data_W00_u(size_data_u), &
                    data_W00_v(size_data_v), &
                    data_W00_w(size_data_w)

    double precision, intent(out) :: source_vector
    dimension :: source_vector(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w)
    
    ! Local data
    ! ------------------
    ! Csr format
    integer :: indi_u_csr, indi_v_csr, indi_w_csr
    dimension ::    indi_u_csr(nb_ctrlpts_u+1), &
                    indi_v_csr(nb_ctrlpts_v+1), &
                    indi_w_csr(nb_ctrlpts_w+1)
    integer, dimension(:), allocatable :: indj_u_csr, indj_v_csr, indj_w_csr
    double precision, dimension(:), allocatable :: data_dummy_csr

    ! ====================================================
    ! Initialize
    allocate(data_dummy_csr(size_data_u), indj_u_csr(size_data_u))
    call coo2csr(nb_ctrlpts_u, size_data_u, data_W00_u, indi_u, indj_u, data_dummy_csr, &
                    indj_u_csr, indi_u_csr)
    deallocate(data_dummy_csr, indj_u_csr)

    allocate(data_dummy_csr(size_data_u), indj_v_csr(size_data_v))
    call coo2csr(nb_ctrlpts_v, size_data_v, data_W00_v, indi_v, indj_v, data_dummy_csr, &
                    indj_v_csr, indi_v_csr)
    deallocate(data_dummy_csr, indj_v_csr)

    allocate(data_dummy_csr(size_data_u), indj_w_csr(size_data_w))
    call coo2csr(nb_ctrlpts_w, size_data_w, data_W00_w, indi_w, indj_w, data_dummy_csr, &
                    indj_w_csr, indi_w_csr)
    deallocate(data_dummy_csr, indj_w_csr)
    ! ====================================================

    ! Find vector
    source_vector = 0.0d0
    call tensor3d_sparsedot_vector(nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, indi_u_csr, indj_u, data_W00_u, &
                        size_data_v, indi_v_csr, indj_v, data_W00_v, &
                        size_data_w, indi_w_csr, indj_w, data_W00_w, &
                        source_coefs, source_vector)

end subroutine wq_get_source_3d

subroutine wq_get_force_3d(nb_qp_total, force_coefs, &
                            nb_ctrlpts_u, nb_qp_u, &
                            nb_ctrlpts_v, nb_qp_v, &
                            nb_ctrlpts_w, nb_qp_w, &
                            size_data_u, size_data_v, size_data_w, &
                            indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w, &
                            data_W00_u, data_W00_v, data_W00_w, &
                            force_vector)
    !! Computes force vector in 3D case

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, intent(in) :: nb_qp_total
    integer, intent(in) ::  nb_ctrlpts_u, nb_qp_u, nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w
    double precision, intent(in) :: force_coefs
    dimension :: force_coefs(nb_qp_total, 3)
    integer, intent(in) :: size_data_u, size_data_v, size_data_w
    integer, intent(in) ::  indi_u, indj_u, &
                            indi_v, indj_v, &
                            indi_w, indj_w
    dimension ::    indi_u(size_data_u), &
                    indj_u(size_data_u), &
                    indi_v(size_data_v), &
                    indj_v(size_data_v), &
                    indi_w(size_data_w), &
                    indj_w(size_data_w)
    double precision, intent(in) :: data_W00_u, data_W00_v, data_W00_w
    dimension ::    data_W00_u(size_data_u), &
                    data_W00_v(size_data_v), &
                    data_W00_w(size_data_w)

    integer, parameter :: dimen = 3
    double precision, intent(out) :: force_vector
    dimension :: force_vector(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w*dimen)

    ! Local data
    ! -------------
    double precision, allocatable, dimension(:) :: force_vector_temp
    integer :: i, init, fin

    ! Initialize
    allocate(force_vector_temp(nb_ctrlpts_u*nb_ctrlpts_v*nb_ctrlpts_w))
    do i = 1, dimen
            
        call wq_get_source_3D(nb_qp_total, force_coefs(:, i), nb_ctrlpts_u, nb_qp_u, &
                        nb_ctrlpts_v, nb_qp_v, nb_ctrlpts_w, nb_qp_w, &
                        size_data_u, size_data_v, size_data_w, &
                        indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, data_W00_u, data_W00_v, data_W00_w, &
                        force_vector_temp)

        init = (i-1)*size(force_vector_temp) + 1
        fin = init + size(force_vector_temp) - 1
        force_vector(init : fin) = force_vector_temp   
    end do

    deallocate(force_vector_temp)
end subroutine wq_get_force_3d

