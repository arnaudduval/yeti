! ==========================
! module :: geometrical properties and behavior of materials
! author :: Joaquin Cornejo
! 
! In this module, one could find functions to compute thermal properties and jacobian given the data at control points
! For mechanical properties, we suggest to go to plasticity module
! It is possible that conductivity and capacity are computed in python, then 
! these functions are deprecated
! ==========================

subroutine eval_jacobien_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, ctrlpts, jacob, detJ, invJ)
    !! Computes jacobien matrix, its determinant and its inverse in 3D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------  
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

    double precision, intent(out) :: jacob, detJ, invJ
    dimension :: jacob(3, 3, nc_u*nc_v*nc_w), detJ(nc_u*nc_v*nc_w), invJ(3, 3, nc_u*nc_v*nc_w)

    ! Local data
    !-----------
    integer :: i, j, k, nb_tasks, beta(d)
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            ctrlpts(i, :), jacob(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
    !$OMP PARALLEL
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        call MatrixInv(invJ(:, :, k), jacob(:, :, k), detJ(k), d)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_3d

subroutine interpolate_fieldphy_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, ddl_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! ------------------- 
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
    !-----------
    integer :: i
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    do i = 1, ddl_u
        call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_fieldphy_3d

subroutine eval_jacobien_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, ctrlpts, &
                            jacob, detJ, invJ)
    !! Computes jacobien matrix, its determinant and its inverse in 2D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 2
    integer, intent(in) :: nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(2, nr_u*nr_v)

    double precision, intent(out) :: jacob, detJ, invJ
    dimension :: jacob(2, 2, nc_u*nc_v), detJ(nc_u*nc_v), invJ(2, 2, nc_u*nc_v)

    ! Local data
    !-----------
    integer :: i, j, k, nb_tasks, beta(2)
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer ::  indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    
    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call sumproduct2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            ctrlpts(i, :), jacob(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJ)/nb_tasks) 
    do k = 1, size(detJ)
        call MatrixInv(invJ(:, :, k), jacob(:, :, k), detJ(k), d)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_2d

subroutine interpolate_fieldphy_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, dof_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------   
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
    !-----------
    integer :: i
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer :: indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    do i = 1, dof_u
        call sumproduct2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_fieldphy_2d

subroutine eval_conductivity_coefficient(dimen, nnz, invJJ, detJJ, nnz_K, Kprop, Kcoefs, info)
    !! Computes conductivity coefficients coef = J^-1 lambda detJ J^-T
    
    use heat_spmf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nnz, nnz_K
    double precision, intent(in) :: invJJ, detJJ
    dimension :: invJJ(dimen, dimen, nnz), detJJ(nnz)
    double precision, target, intent(in) :: Kprop
    dimension :: Kprop(dimen, dimen, nnz_K)

    integer, intent(out) :: info
    double precision, intent(out) :: Kcoefs
    dimension :: Kcoefs(dimen, dimen, nnz)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat

    allocate(mat)
    mat%Kprop => Kprop
    allocate(mat%Kcoefs(dimen, dimen, nnz))
    call update_conductivity_coefs(mat, dimen, nnz, invJJ, detJJ, info)
    Kcoefs = mat%Kcoefs
    
end subroutine eval_conductivity_coefficient

subroutine eval_capacity_coefficient(nnz, detJJ, nnz_C, Cprop, Ccoefs, info)
    !! Computes capacity coefficient coef = sigma * detJ
    
    use heat_spmf
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, intent(in) :: nnz, nnz_C
    double precision, intent(in) :: detJJ
    dimension :: detJJ(nnz)
    double precision, target, intent(in) :: Cprop
    dimension :: Cprop(nnz_C)

    integer, intent(out) :: info
    double precision, intent(out) :: Ccoefs
    dimension :: Ccoefs(nnz)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat

    allocate(mat)
    mat%Cprop => Cprop
    allocate(mat%Ccoefs(nnz))
    call update_capacity_coefs(mat, nnz, detJJ, info)
    Ccoefs = mat%Ccoefs

end subroutine eval_capacity_coefficient