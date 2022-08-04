! ==========================
! module :: Assembly for IGA-WQ 
! author :: Joaquin Cornejo
! ==========================

subroutine eval_conductivity_coefficient(dime, nnz, JJ, nnz_K, KK, Kcoef, info)
    !! Computes coefficient for K
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------    
    integer, intent(in) :: dime, nnz, nnz_K
    double precision, intent(in) :: JJ, KK
    dimension :: JJ(dime, dime, nnz), KK(dime, dime, nnz_K)

    integer, intent(out) :: info
    double precision, intent(out) :: Kcoef
    dimension :: Kcoef(dime, dime, nnz)

    ! Local data
    ! -----------
    integer :: i, nb_tasks
    double precision :: Jt, detJt, invJt
    dimension :: Jt(dime, dime), invJt(dime, dime)
    double precision :: Kt
    dimension :: Kt(dime, dime)  

    info = 1

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

end subroutine eval_conductivity_coefficient

subroutine eval_capacity_coefficient(dime, nnz, JJ, nnz_C, CC, Ccoef, info)
    !! Computes coefficient for C 
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! --------------------    
    integer, intent(in) :: dime, nnz, nnz_C
    double precision, intent(in) :: JJ, CC
    dimension :: JJ(dime, dime, nnz), CC(nnz_C)

    integer, intent(out) :: info
    double precision, intent(out) :: Ccoef
    dimension :: Ccoef(nnz)

    ! Local data
    ! -----------
    integer :: i, nb_tasks
    double precision :: Jt, detJt
    dimension :: Jt(dime, dime)

    info = 1

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

end subroutine eval_capacity_coefficient

subroutine jacobien_physicalposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, ctrlpts, &
                                        jacob, physical_pos, detJ, invJ)
    !! Computes jacobien matrix and physical position in 3D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, parameter :: d = 3
    integer, intent(in) :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(3, nr_u*nr_v*nr_w)

    double precision, intent(out) :: jacob, physical_pos, detJ, invJ
    dimension ::    jacob(3, 3, nc_u*nc_v*nc_w), &
                    physical_pos(3, nc_u*nc_v*nc_w), &
                    detJ(nc_u*nc_v*nc_w), invJ(3, 3, nc_u*nc_v*nc_w)

    ! Local data
    !-----------------
    integer :: i, j, k, nb_tasks, beta
    dimension :: beta(d)
    double precision :: result
    dimension :: result(nc_u*nc_v*nc_w)

    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Compute physical position
    !$OMP PARALLEL PRIVATE(result)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, d/nb_tasks) 
    do i = 1, d
        call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            ctrlpts(i, :), result)
        physical_pos(i, :) = result
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta, result)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            ctrlpts(i, :), result)
            jacob(i, j, :) = result
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
    ! For det J 
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        ! Evaluate determinant
        call MatrixInv(invJ(:, :, k), jacob(:, :, k), detJ(k), d)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_3d

subroutine interpolation_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, ddl_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, intent(in) ::  nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, ddl_u
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(ddl_u, nr_u*nr_v*nr_w)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(ddl_u, nc_u*nc_v*nc_w)

    ! Local data
    !-----------------
    integer :: i
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Interpolation
    do i = 1, ddl_u
        call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                                    nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                    nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                                    u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolation_3d

subroutine jacobien_physicalposition_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                        indi_u, indj_u, indi_v, indj_v, &
                                        data_B_u, data_B_v, ctrlpts, &
                                        jacob, physical_pos, detJ, invJ)
    !! Computes jacobien in 2D case
    !! IN CSR FORMAT
    
    use omp_lib
    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, parameter :: d = 2
    integer, intent(in) :: nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(2, nr_u*nr_v)

    double precision, intent(out) :: jacob, physical_pos, detJ, invJ
    dimension ::    jacob(2, 2, nc_u*nc_v), &
                    physical_pos(2, nc_u*nc_v), &
                    detJ(nc_u*nc_v), invJ(2, 2, nc_u*nc_v)

    ! Local data
    !-----------------
    integer :: i, j, k, nb_tasks, beta
    dimension :: beta(2)
    double precision :: result
    dimension ::  result(nc_u*nc_v)

    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer ::  indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    
    ! Compute physical position
    !$OMP PARALLEL PRIVATE(result)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, d/nb_tasks) 
    do i = 1, d
        call tensor2d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            ctrlpts(i, :), result)
        physical_pos(i, :) = result
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta, result)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call tensor2d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            ctrlpts(i, :), result)
            jacob(i, j, :) = result
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    ! For det J 
    !$OMP PARALLEL PRIVATE(detJt)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        ! Evaluate determinant
        call MatrixInv(invJ(:, :, k), jacob(:, :, k), detJ(k), d)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine jacobien_physicalposition_2d

subroutine interpolation_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, dof_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------   
    integer, intent(in) ::  nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v, dof_u
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(dof_u, nr_u*nr_v)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(dof_u, nc_u*nc_v)

    ! Local data
    !-----------------
    integer :: i
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer :: indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    ! Interpolation
    do i = 1, dof_u
        call tensor2d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, &
                                    nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                                    nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                    u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolation_2d