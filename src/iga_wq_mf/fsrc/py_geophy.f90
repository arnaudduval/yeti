! ==============================================
! In this module, one could find functions to compute thermal properties and jacobian given the data at control points
! For mechanical properties, we suggest to go to plasticity module
! ==============================================

subroutine eval_jacobien_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, ctrlpts, JJ, detJJ, invJJ)
    !! Computes jacobien matrix, its determinant and its inverse in 3D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(dimen, nr_u*nr_v*nr_w)

    double precision, intent(out) :: JJ, detJJ, invJJ
    dimension :: JJ(dimen, dimen, nc_u*nc_v*nc_w), detJJ(nc_u*nc_v*nc_w), invJJ(dimen, dimen, nc_u*nc_v*nc_w)

    ! Local data
    !-----------
    integer :: i, j, k, nb_tasks, beta(dimen)
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
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, dimen*dimen/nb_tasks) 
    do j = 1, dimen
        do i = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            ctrlpts(i, :), JJ(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
    !$OMP PARALLEL
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJJ)/nb_tasks) 
    do k = 1, size(detJJ)
        call MatrixInv(invJJ(:, :, k), JJ(:, :, k), detJJ(k), dimen)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_3d

subroutine interpolate_meshgrid_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, nm, u_ctrlpts, u_interp)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) ::  nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, nm
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v*nr_w)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nm, nc_u*nc_v*nc_w)

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

    do i = 1, nm
        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_meshgrid_3d

subroutine l2projection_ctrlpts_3d(coefs, nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use matrixfreeheat
    use solverheat3
    use datastructure

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nm, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nm, nr_total), resPCG(nm, nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    type(structure), allocatable :: struct
    double precision :: ones(3)
    integer :: i

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    ones = 1.d0
    allocate(struct)
    call init_3datastructure(struct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)
    call getcsr2csc(struct)
    call eigendecomposition(struct, ones)

    ! Set material and solver
    allocate(mat, solv)
    call setup_capacitycoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 1

    do i = 1, nm
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    struct%indi_T(1, 1:nc_u+1), struct%indj_T(1, 1:nnz_u), struct%indi_T(2, 1:nc_v+1), &
                    struct%indj_T(2, 1:nnz_v), struct%indi_T(3, 1:nc_w+1), struct%indj_T(3, 1:nnz_w), &
                    struct%bw_T(1, 1:nnz_u, :2), struct%bw_T(2, 1:nnz_v, :2), struct%bw_T(3, 1:nnz_w, :2), &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                    struct%eigvec(1, 1:nr_u, 1:nr_u), struct%eigvec(2, 1:nr_v, 1:nr_v), &
                    struct%eigvec(3, 1:nr_w, 1:nr_w), nbIterPCG, threshold, b(i, :), x(i, :), resPCG(i, :))
    end do

end subroutine l2projection_ctrlpts_3d

subroutine eval_jacobien_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, ctrlpts, &
                            JJ, detJJ, invJJ)
    !! Computes jacobien matrix, its determinant and its inverse in 2D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: ctrlpts
    dimension :: ctrlpts(dimen, nr_u*nr_v)

    double precision, intent(out) :: JJ, detJJ, invJJ
    dimension :: JJ(dimen, dimen, nc_u*nc_v), detJJ(nc_u*nc_v), invJJ(dimen, dimen, nc_u*nc_v)

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
    !$OMP DO SCHEDULE(STATIC, dimen*dimen/nb_tasks) 
    do j = 1, dimen
        do i = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            ctrlpts(i, :), JJ(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detJJ)/nb_tasks) 
    do k = 1, size(detJJ)
        call MatrixInv(invJJ(:, :, k), JJ(:, :, k), detJJ(k), dimen)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_2d

subroutine interpolate_meshgrid_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, nm, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------   
    integer, intent(in) ::  nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v, nm
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nm, nc_u*nc_v)

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

    do i = 1, nm
        call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_meshgrid_2d

subroutine l2projection_ctrlpts_2d(coefs, nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use matrixfreeheat
    use solverheat2
    use datastructure

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nm, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nm, nr_total), resPCG(nm, nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    type(structure), allocatable :: struct
    double precision :: ones(2)
    integer :: i

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    ones = 1.d0
    allocate(struct)
    call init_2datastructure(struct, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v)
    call getcsr2csc(struct)
    call eigendecomposition(struct, ones)

    ! Set material and solver
    allocate(mat, solv)
    call setup_capacitycoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 1

    do i = 1, nm
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    struct%indi_T(1, 1:nc_u+1), struct%indj_T(1, 1:nnz_u), struct%indi_T(2, 1:nc_v+1), &
                    struct%indj_T(2, 1:nnz_v), struct%bw_T(1, 1:nnz_u, :2), struct%bw_T(2, 1:nnz_v, :2), &
                    indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, &
                    struct%eigvec(1, 1:nr_u, 1:nr_u), struct%eigvec(2, 1:nr_v, 1:nr_v), &
                    nbIterPCG, threshold, b(i, :), x(i, :), resPCG(i, :))
    end do

end subroutine l2projection_ctrlpts_2d


!! --------- TO ERASE EVENTUALLY

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
