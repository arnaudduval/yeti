! ==========================
! module :: geometrical properties and behavior of materials
! author :: Joaquin Cornejo
! 
! In this module, one could find functions to compute thermal properties and jacobian given the data at control points
! For mechanical properties, we suggest to go to plasticity module
! It is possible that conductivity and capacity are computed in python, then 
! these functions are deprecated
! ==========================

subroutine tensor_n_mode_product_py(ncols, ncX, X, &
                                    nrU, nnzU, dataU, indi, indj, &
                                    mode, isdense, istransp, nrR, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = (nc_u, nc_v, nc_w)
    !! Matrix U = (nr, nc). Since U is in CSR format, nc is not necessary to be declared
    !! Tensor R = (nu, nv, nw) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw
    !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u
    !! This function is adapted to python

    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: ncols, ncX, nrU, nnzU, mode, nrR
    dimension :: ncols(3)
    integer, intent(in) :: indi, indj
    dimension :: indi(nrU+1), indj(nnzU)
    double precision, intent(in) :: X, dataU
    dimension :: X(ncX), dataU(nnzU)
    logical, intent(in) :: isdense, istransp

    double precision, intent(out) :: R
    dimension :: R(nrR)

    ! Local data
    ! ----------
    integer :: nrows(3)
    double precision, allocatable, dimension(:, :) :: UU
    integer, allocatable, dimension(:) :: indi_T, indj_T
    double precision, allocatable, dimension(:, :) :: dataUt, dataU_T

    ! Set conditions
    if (any(indi.le.0)) stop 'Indices in Fortran format are needed'
    if (any(indj.le.0)) stop 'Indices in Fortran format are needed'
    if ((mode.lt.1).or.(mode.gt.3)) stop 'Mode must be greater than 0 and less than 4'
    if (product(ncols).ne.ncX) stop 'Wrong dimensions'
    nrows = ncols; nrows(mode) = nrU
    if (product(nrows).ne.nrR) stop 'Wrong dimensions'

    if (isdense) then
        ! Create dense matrix
        allocate(UU(nrU, ncols(mode)))
        call csr2dense(nnzU, indi, indj, dataU, nrU, ncols(mode), UU)

        if (istransp) then
            ! Compute tensor n mode product of transpose matrix
            call tensor_n_mode_product(ncols(1), ncols(2), ncols(3), X, ncols(mode), &
                                    nrU, transpose(UU), mode, nrR, R)
        else
            ! Compute tensor n mode product
            call tensor_n_mode_product(ncols(1), ncols(2), ncols(3), X, nrU, &
                                    ncols(mode), UU, mode, nrR, R)
        end if

    else
        if (istransp) then
            ! Convert CSR to CSC
            allocate(indi_T(ncols(mode)+1), indj_T(nnzU), dataUt(nnzU, 1), dataU_T(nnzU, 1))
            dataUt(:, 1) = dataU
            call csr2csc(1, nrU, ncols(mode), nnzU, dataUt, indj, indi, dataU_T, indj_T, indi_T)

            ! Compute tensor n mode product
            call tensor_n_mode_product_sp(ncols(1), ncols(2), ncols(3), X, ncols(mode), &
                                        nnzU, dataU_T(:, 1), indi_T, indj_T, mode, nrR, R)
        
        else
            ! Compute tensor n mode product
            call tensor_n_mode_product_sp(ncols(1), ncols(2), ncols(3), X, nrU, &
                                        nnzU, dataU, indi, indj, mode, nrR, R)
        end if

    end if

end subroutine tensor_n_mode_product_py

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
    double precision :: temp
    dimension :: temp(nc_u*nc_v*nc_w)

    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta, temp)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call sumproduct3d_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            ctrlpts(i, :), temp)
            jacob(i, j, :) = temp
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
    ! Compute inverse and determinant 
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

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Interpolation
    do i = 1, ddl_u
        call sumproduct3d_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
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
    double precision :: temp
    dimension :: temp(nc_u*nc_v)

    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer ::  indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    
    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta, temp)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, d*d/nb_tasks) 
    do j = 1, d
        do i = 1, d
            beta = 1; beta(j) = 2
            call sumproduct2d_sp(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            ctrlpts(i, :), temp)
            jacob(i, j, :) = temp
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    ! Compute inverse and determinant
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

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    ! Interpolation
    do i = 1, dof_u
        call sumproduct2d_sp(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_fieldphy_2d

subroutine eval_conductivity_coefficient(dime, nnz, JJ, nnz_K, KK, Kcoef, info)
    !! Computes conductivity coefficients coef = J^-1 lambda detJ J^-T
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dime, nnz, nnz_K
    double precision, intent(in) :: JJ, KK
    dimension :: JJ(dime, dime, nnz), KK(dime, dime, nnz_K)

    integer, intent(out) :: info
    double precision, intent(out) :: Kcoef
    dimension :: Kcoef(dime, dime, nnz)

    ! Local data
    ! ----------
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
            ! Get jacobien
            Jt = JJ(:, :, i)

            ! Evaluate inverse
            call MatrixInv(invJt, Jt, detJt, dime)

            ! Compute K = invJ * prop * detJ * invJ'
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
            ! Get jacobien
            Jt = JJ(:, :, i)

            ! Evaluate inverse
            call MatrixInv(invJt, Jt, detJt, dime)

            ! Compute K = invJ * prop * detJ * invJ'
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
    !! Computes capacity coefficient coef = sigma * detJ
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, intent(in) :: dime, nnz, nnz_C
    double precision, intent(in) :: JJ, CC
    dimension :: JJ(dime, dime, nnz), CC(nnz_C)

    integer, intent(out) :: info
    double precision, intent(out) :: Ccoef
    dimension :: Ccoef(nnz)

    ! Local data
    ! ----------
    integer :: i, nb_tasks
    double precision :: Jt, detJt
    dimension :: Jt(dime, dime)

    info = 1

    if (nnz_C.eq.1) then 

        !$OMP PARALLEL PRIVATE(Jt, detJt)
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
        do i = 1, nnz
            ! Get jacobien
            Jt = JJ(:, :, i)

            ! Compute C = detJ  * prop
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
            ! Get jacobien
            Jt = JJ(:, :, i)

            ! Compute C = detJ  * prop
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