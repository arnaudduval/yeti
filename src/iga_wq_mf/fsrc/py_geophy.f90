! ==============================================
! In this module, one could find functions to compute thermal properties and jacobian given the data at control points
! For mechanical properties, we suggest to go to plasticity module
! ==============================================

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
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
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
        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
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
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
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
                            data_B_u, data_B_v, ddl_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------   
    integer, intent(in) ::  nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v, ddl_u
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(ddl_u, nr_u*nr_v)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(ddl_u, nc_u*nc_v)

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

    do i = 1, ddl_u
        call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_fieldphy_2d

subroutine eval_capacity_coefficient(nnzJ, detJ, nnzP, prop, coefs, info)
    !! Computes capacity coefficient coef = sigma * detJ
    
    use matrixfreeheat
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, intent(in) :: nnzJ, nnzP
    double precision, intent(in) :: detJ
    dimension :: detJ(nnzJ)
    double precision, target, intent(in) :: prop
    dimension :: prop(nnzP)

    integer, intent(out) :: info
    double precision, intent(out) :: coefs
    dimension :: coefs(nnzJ)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    double precision :: invJJ(1, 1, nnzJ)
    
    allocate(mat)
    mat%dimen = 1; invJJ = 0.d0
    call setup_geometry(mat, nnzJ, invJJ, detJ)
    mat%Cprop => prop
    allocate(mat%Ccoefs(nnzJ))
    call update_capacitycoefs(mat, info)
    coefs = mat%Ccoefs

end subroutine eval_capacity_coefficient

subroutine eval_conductivity_coefficient(dimen, nnzJ, invJ, detJ, nnzP, prop, coefs, info)
    !! Computes conductivity coefficients coef = J^-1 lambda detJ J^-T
    
    use matrixfreeheat
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nnzJ, nnzP
    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nnzJ), detJ(nnzJ)
    double precision, target, intent(in) :: prop
    dimension :: prop(dimen, dimen, nnzP)

    integer, intent(out) :: info
    double precision, intent(out) :: coefs
    dimension :: coefs(dimen, dimen, nnzJ)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat

    allocate(mat)
    mat%dimen = dimen
    call setup_geometry(mat, nnzJ, invJ, detJ)
    mat%Kprop => prop
    allocate(mat%Kcoefs(dimen, dimen, nnzJ))
    call update_conductivitycoefs(mat, info)
    coefs = mat%Kcoefs
    
end subroutine eval_conductivity_coefficient

subroutine eval_intforce_coefficient(dimen, nnz, invJ, detJ, stress, coefs)
    
    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nnz
    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nnz), detJ(nnz)
    double precision, target, intent(in) :: stress
    dimension :: stress(dimen, nnz)

    double precision, intent(out) :: coefs
    dimension :: coefs(dimen, dimen, nnz)

    ! Local data
    ! ----------
    type(mecamat), pointer :: mat
    double precision :: Tstress
    dimension :: Tstress(dimen, dimen)
    integer :: i

    allocate(mat)
    mat%dimen  = dimen
    mat%nvoigt = dimen*(dimen+1)/2
    call setup_geometry(mat, nnz, invJ, detJ)

    do i = 1, nnz
        call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
        coefs(:,:,i) = matmul(mat%invJ(:,:,i), Tstress)*mat%detJ(i)
    end do

end subroutine eval_intforce_coefficient

subroutine eigen_decomposition_py(nr, nc, nnz, indi, indj, data_B, data_W, &
                                Mcoefs, Kcoefs, robcond1, robcond2, eigenvalues, eigenvectors)
    !! Eigen decomposition generalized KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B, data_W
    dimension :: data_B(nnz, 2), data_W(nnz, 4)
    double precision, intent(in) :: Mcoefs, Kcoefs
    dimension :: Mcoefs(nc), Kcoefs(nc)
    integer, intent(in) :: robcond1, robcond2
    dimension :: robcond1(2), robcond2(2)
            
    double precision, intent(out) :: eigenvalues, eigenvectors
    dimension :: eigenvalues(nr), eigenvectors(nr, nr)

    ! Local data
    ! ----------
    double precision :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    call eigen_decomposition(nr, nc, Mcoefs, Kcoefs, nnz, indi, indj, data_B, data_W, &
                            robcond1, robcond2, eigenvalues, eigenvectors, Kdiag, Mdiag)

end subroutine eigen_decomposition_py
