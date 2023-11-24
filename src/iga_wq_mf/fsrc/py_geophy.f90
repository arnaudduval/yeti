! ==============================================
! In this module, one could find functions to compute thermal properties and jacobian given the data at control points
! For mechanical properties, we suggest to go to plasticity module
! ==============================================

subroutine stiffmass_eigendecomposition_dense(nr, A, B, eigval, eigvec)
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double precision, intent(in) :: A, B 
    dimension :: A(nr, nr), B(nr, nr)
            
    double precision, intent(out) :: eigval, eigvec
    dimension :: eigval(nr), eigvec(nr, nr)

    ! Local data
    ! ----------
    call compute_eigdecomp_pdr(nr, A, B, eigval, eigvec)

end subroutine stiffmass_eigendecomposition_dense

subroutine stiffmass_eigendecomposition_sp(nr, nc, univMcoefs, univKcoefs, nnz, indi, indj, &
                                data_B, data_W, eigval, eigvec)
    !! Generalized eigen decomposition KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    double precision, intent(in) :: univMcoefs, univKcoefs
    dimension :: univMcoefs(nc), univKcoefs(nc)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B, data_W
    dimension :: data_B(nnz, 2), data_W(nnz, 4)
            
    double precision, intent(out) :: eigval, eigvec
    dimension :: eigval(nr), eigvec(nr, nr)

    ! Local data
    ! ----------
    double precision :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    call stiffmass_eigendecomposition(nr, nc, univMcoefs, univKcoefs, nnz, indi, indj, &
                            data_B, data_W, eigval, eigvec, Kdiag, Mdiag)

end subroutine stiffmass_eigendecomposition_sp

subroutine eval_normal(nr, nc_total, JJ, normal)
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc_total
    double precision, intent(in) :: JJ
    dimension :: JJ(nr, nr-1, nc_total)

    double precision, intent(out) :: normal 
    dimension :: normal(nr, nc_total)

    ! Local data
    ! ----------
    integer :: i
    double precision, dimension(3) :: t1, t2, t3
    double precision :: norm

    if ((nr.lt.2).or.(nr.gt.3)) stop 'Size problem'
    t1 = 0.d0; t2 = 0.d0; normal = 0.d0
    if (nr.eq.3) then 
        do i = 1, nc_total
            t1 = JJ(:, 1, i); t2 = JJ(:, 2, i)
            call crossproduct(t1, t2, t3)
            norm = norm2(t3)
            normal(:, i) = t3/norm
        end do
    else
        t2(3) = 1.d0
        do i = 1, nc_total
            t1(:nr) = JJ(:, 1, i)
            call crossproduct(t1, t2, t3)
            norm = norm2(t3)
            normal(:, i) = t3(:nr)/norm
        end do
    end if
end subroutine eval_normal

subroutine eval_inverse_det(nr, nnz, M, detM, invM)
    !! Computes the inverse and determinant of a matrix in a format (:, :, nnz)
    use omp_lib
    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nr, nnz
    double precision, intent(in) :: M
    dimension :: M(nr, nr, nnz)

    double precision, intent(out) :: invM, detM
    dimension :: invM(nr, nr, nnz), detM(nnz)

    ! Local data
    !-----------
    integer :: i, nb_tasks

    !$OMP PARALLEL
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(detM)/nb_tasks) 
    do i = 1, size(detM)
        call MatrixInv(invM(:, :, i), M(:, :, i), detM(i), nr)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_inverse_det

subroutine eval_jacobien_3d(nm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, u_ctrlpts, JJ)
    !! Computes jacobien matrix, its determinant and its inverse in 3D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: dimen = 3
    integer, intent(in) :: nm, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v*nr_w)

    double precision, intent(out) :: JJ
    dimension :: JJ(nm, dimen, nc_u*nc_v*nc_w)

    ! Local data
    !-----------
    integer :: i, j, nb_tasks, beta(dimen)
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
        do i = 1, nm
            beta = 1; beta(j) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            u_ctrlpts(i, :), JJ(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_3d

subroutine interpolate_meshgrid_3d(nm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, u_ctrlpts, u_interp)
    !! Computes interpolation in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nm, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
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

subroutine l2projection_ctrlpts_3d(nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, detJ, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use matrixfreeheat
    use heatsolver3
    use datastructure

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, target, intent(in) :: detJ
    dimension :: detJ(nc_total)
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nm, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nm, nr_total), resPCG(nm, nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: i 

    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    logical :: table(dimen, 2)   
    double precision :: mean(dimen) = 1.d0
    double precision :: ones(nc_total)

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set material and solver
    mat%dimen = dimen; ones = 1.d0
    mat%detJ => detJ; mat%ncols_sp = nc_total
    call setup_capacityprop(mat, nc_total, ones)
    solv%matrixfreetype = 1

    table = .false.; solv%withdiag = .false.
    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                            data_W_u, data_W_v, data_W_w, table, mean)

    do i = 1, nm
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b(i, :), x(i, :), resPCG(i, :))
    end do
                
end subroutine l2projection_ctrlpts_3d

subroutine eval_jacobien_2d(nm, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, u_ctrlpts, JJ)
    !! Computes jacobien matrix, its determinant and its inverse in 2D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nm, nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v)

    double precision, intent(out) :: JJ
    dimension :: JJ(nm, dimen, nc_u*nc_v)

    ! Local data
    !-----------
    integer :: i, j, nb_tasks, beta(2)
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
        do i = 1, nm
            beta = 1; beta(j) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            u_ctrlpts(i, :), JJ(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_2d

subroutine interpolate_meshgrid_2d(nm, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------   
    integer, intent(in) :: nm, nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
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

subroutine l2projection_ctrlpts_2d(nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, detJ, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use matrixfreeheat
    use heatsolver2
    use datastructure

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nm, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    double precision, target, intent(in) :: detJ
    dimension :: detJ(nc_total)
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nm, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nm, nr_total), resPCG(nm, nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: i

    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)
    logical :: table(dimen, 2)   
    double precision :: mean(dimen) = 1.d0
    double precision :: ones(nc_total)

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    mat%dimen = dimen; ones = 1.d0
    mat%detJ => detJ; mat%ncols_sp = nc_total
    call setup_capacityprop(mat, nc_total, ones)
    solv%matrixfreetype = 1

    table = .false.; solv%withdiag = .false.
    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                            data_W_u, data_W_v, table, mean)

    do i = 1, nm
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, nbIterPCG, threshold, b(i, :), x(i, :), resPCG(i, :))
    end do

end subroutine l2projection_ctrlpts_2d

subroutine eval_jacobien_1d(nm, nr_u, nc_u, nnz_u, indi_u, indj_u, &
                            data_B_u, u_ctrlpts, JJ)
    !! Computes jacobien matrix, its determinant and its inverse in 2D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 1
    integer, intent(in) :: nm, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_B_u
    dimension :: data_B_u(nnz_u, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u)

    double precision, intent(out) :: JJ
    dimension :: JJ(nm, dimen, nc_u)

    ! Local data
    !-----------
    integer :: i
    integer :: indi_T_u, indj_T_u
    dimension :: indi_T_u(nc_u+1), indj_T_u(nnz_u)
    double precision :: data_BT_u
    dimension :: data_BT_u(nnz_u, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    
    do i = 1, nm
        call spmat_dot_dvec(nc_u, nr_u, nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 2), &
                        u_ctrlpts(i, :), JJ(i, 1, :))
    end do

end subroutine eval_jacobien_1d

subroutine interpolate_meshgrid_1d(nm, nr_u, nc_u, nnz_u, indi_u, indj_u, &
                            data_B_u, u_ctrlpts, u_interp)
    !! Computes interpolation in 2D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------   
    integer, intent(in) :: nm, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_B_u
    dimension :: data_B_u(nnz_u, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nm, nc_u)

    ! Local data
    !-----------
    integer :: i
    integer :: indi_T_u, indj_T_u
    dimension :: indi_T_u(nc_u+1), indj_T_u(nnz_u)
    double precision :: data_BT_u
    dimension :: data_BT_u(nnz_u, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)

    do i = 1, nm
        call spmat_dot_dvec(nc_u, nr_u, nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_meshgrid_1d

subroutine l2projection_ctrlpts_1d(nm, nr_total, nc_total, nr_u, nc_u, &
                            nnz_u, indi_u, indj_u, data_B_u, data_W_u, detJ, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nr_total, nc_total, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_B_u, data_W_u
    dimension :: data_B_u(nnz_u, 2), data_W_u(nnz_u, 4)

    double precision, intent(in) :: detJ
    dimension :: detJ(nc_total)
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nm, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nm, nr_total), resPCG(nm, nbIterPCG+1)

    ! Local data
    ! ----------
    integer :: i 
    double precision :: data_Wt
    dimension :: data_Wt(nnz_u)
    double precision :: WW, BB, A
    dimension :: WW(nr_u, nc_u), BB(nr_u, nc_u), A(nr_u, nr_u)

    if (nr_total.ne.nr_u) stop 'Size problem'
    resPCG = threshold
    do i = 1, nnz_u
        data_Wt(i) = data_W_u(i, 1) * detJ(indj_u(i)) 
    end do
    call csr2dense(nnz_u, indi_u, indj_u, data_Wt, nr_u, nc_u, WW)
    call csr2dense(nnz_u, indi_u, indj_u, data_B_u(:, 1), nr_u, nc_u, BB)
    A = matmul(WW, transpose(BB))
    do i = 1, nm
        call solve_linear_system(nr_u, nr_u, A, b(i, :), x(i, :))
    end do

end subroutine l2projection_ctrlpts_1d

subroutine get_forcevol_2d(nm, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, detJ, prop, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v)

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp(nc_total)

    do i = 1, nm
        tmp = prop(i, :)*detJ
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            tmp, array_out(i, :))
    end do

end subroutine get_forcevol_2d

subroutine get_forcevol_3d(nm, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, detJ, prop, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp(nc_total)

    do i = 1, nm
        tmp = prop(i, :)*detJ
        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            tmp, array_out(i, :))
    end do

end subroutine get_forcevol_3d

subroutine get_forcesurf_2d(nm, nc_total, nr_u, nc_u, nnz_u, indi_u, indj_u, data_W_u, JJ, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nm, nc_total, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_W_u
    dimension :: data_W_u(nnz_u, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, W00
    dimension :: coefs(nm, nc_total), v1(dimen), W00(nr_u, nc_total)
    integer :: i

    if (nc_total.ne.nc_u) stop 'Size problem'

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        dsurf = sqrt(dot_product(v1, v1))
        coefs(:, i) = prop(:, i)*dsurf
    end do

    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W00)
    do i = 1, nm
        array_out(i, :) = matmul(W00, coefs(i, :))
    end do
    
end subroutine get_forcesurf_2d

subroutine get_forcesurf_3d(nm, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, JJ, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nm, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(nm, nc_total), v1(dimen), v2(dimen), v3(dimen)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))
        coefs(:, i) = prop(:, i)*dsurf
    end do

    do i = 1, nm
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs(i, :), array_out(i, :))
    end do
    
end subroutine get_forcesurf_3d

!! SPACE TIME APPROACH

subroutine eval_jacobien_4d(nm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, &
                            data_B_u, data_B_v, data_B_w, data_B_t, u_ctrlpts, JJ)
    !! Computes jacobien matrix, its determinant and its inverse in 3D
    !! IN CSR FORMAT
    
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: dimen = 4
    integer, intent(in) :: nm, nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_B_t
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), data_B_t(nnz_t, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v*nr_w*nr_t)

    double precision, intent(out) :: JJ
    dimension :: JJ(nm, dimen, nc_u*nc_v*nc_w*nc_t)

    ! Local data
    !-----------
    integer :: i, j, nb_tasks, beta(dimen)
    integer :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
    integer :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
    double precision :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    call csr2csc(2, nr_t, nc_t, nnz_t, data_B_t, indj_t, indi_t, data_BT_t, indj_T_t, indi_T_t)

    ! Compute jacobien matrix
    !$OMP PARALLEL PRIVATE(beta)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, dimen*dimen/nb_tasks) 
    do j = 1, dimen
        do i = 1, nm
            beta = 1; beta(j) = 2
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, beta(4)), &
                            u_ctrlpts(i, :), JJ(i, j, :))
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine eval_jacobien_4d

subroutine interpolate_meshgrid_4d(nm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, &
                            data_B_u, data_B_v, data_B_w, data_B_t, u_ctrlpts, u_interp)
    !! Computes interpolation in 4D case (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nm, nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_B_t
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), data_B_t(nnz_t, 2)
    double precision, intent(in) :: u_ctrlpts
    dimension :: u_ctrlpts(nm, nr_u*nr_v*nr_w*nr_t)

    double precision, intent(out) :: u_interp
    dimension :: u_interp(nm, nc_u*nc_v*nc_w*nc_t)

    ! Local data
    !-----------
    integer :: i
    integer :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
    integer :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
    double precision :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    call csr2csc(2, nr_t, nc_t, nnz_t, data_B_t, indj_t, indi_t, data_BT_t, indj_T_t, indi_T_t)

    do i = 1, nm
        call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), &
                            u_ctrlpts(i, :), u_interp(i, :))
    end do

end subroutine interpolate_meshgrid_4d

subroutine get_forcevol_st2d(nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, &
                            data_W_u, data_W_v, data_W_t, detJ, detG, prop, array_out)
    !! Computes volumetric force vector in 4D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_t
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)
    double precision, intent(in) :: detJ, detG, prop
    dimension :: detJ(nc_sp), detG(nc_tm), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v*nr_t)

    ! Local data
    ! ----------
    integer :: l, i, j, k
    double precision :: tmp(nc_total)

    if ((nc_total.ne.nc_sp*nc_tm).or.(nc_total.ne.nc_u*nc_v*nc_t)) stop 'Size problem'
    do l = 1, nm

        do j = 1, nc_t
            do i = 1, nc_u*nc_v
                k = i + (j-1)*nc_u*nc_v
                tmp(k) = prop(l, k)*detJ(i)*detG(j)
            end do
        end do
        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                            tmp, array_out(l, :))
    end do

end subroutine get_forcevol_st2d

subroutine get_forcevol_st3d(nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, &
                            nc_t, nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, &
                            indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, detJ, detG, prop, array_out)
    !! Computes volumetric force vector in 4D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)
    double precision, intent(in) :: detJ, detG, prop
    dimension :: detJ(nc_sp), detG(nc_tm), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v*nr_w*nr_t)

    ! Local data
    ! ----------
    integer :: l, i, j, k
    double precision :: tmp(nc_total)
    if ((nc_total.ne.nc_sp*nc_tm).or.(nc_total.ne.nc_u*nc_v*nc_w*nc_t)) stop 'Size problem'

    do l = 1, nm

        do j = 1, nc_t
            do i = 1, nc_u*nc_v*nc_w
                k = i + (j-1)*nc_u*nc_v*nc_w
                tmp(k) = prop(l, k)*detJ(i)*detG(j)
            end do
        end do

        call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                            tmp, array_out(l, :))
    end do

end subroutine get_forcevol_st3d

subroutine get_forcesurf_st2d(nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_t, nc_t, nnz_u, nnz_t, &
                            indi_u, indj_u, indi_t, indj_t, data_W_u, data_W_t, &
                            JJ, detG, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 2
    integer, intent(in) :: nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_t, nc_t, nnz_u, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_W_u, data_W_t
    dimension :: data_W_u(nnz_u, 4), data_W_t(nnz_t, 4)
    double precision, intent(in) :: JJ, detG, prop
    dimension :: JJ(dimen_sp, dimen_sp-1, nc_sp), detG(nc_tm), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_t)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1
    dimension :: coefs(nm, nc_total), v1(dimen_sp)
    integer :: i, j, k
    if (nc_total.ne.nc_sp*nc_tm) stop 'Size problem'

    ! Compute coefficients
    do i = 1, nc_u
        v1 = JJ(:, 1, i)
        dsurf = sqrt(dot_product(v1, v1))
        do j = 1, nc_t
            k = i + (j-1)*nc_u
            coefs(:, k) = prop(:, k)*dsurf*detG(j)
        end do
    end do

    do i = 1, nm
        call sumfacto2d_spM(nr_u, nc_u, nr_t, nc_t, nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 1), coefs(i, :), array_out(i, :))
    end do
    
end subroutine get_forcesurf_st2d

subroutine get_forcesurf_st3d(nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, data_W_t, &
                            JJ, detG, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nm, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_t
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)
    double precision, intent(in) :: JJ, detG, prop
    dimension :: JJ(dimen_sp, dimen_sp-1, nc_sp), detG(nc_tm), prop(nm, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nm, nr_u*nr_v*nr_t)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(nm, nc_total), v1(dimen_sp), v2(dimen_sp), v3(dimen_sp)
    integer :: i, j, k

    if (nc_total.ne.nc_sp*nc_tm) stop 'Size problem'

    ! Compute coefficients
    do i = 1, nc_u*nc_v
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))
        do j = 1, nc_t
            k = i + (j-1)*nc_u*nc_v
            coefs(:, k) = prop(:, k)*dsurf*detG(j)
        end do
    end do

    do i = 1, nm
        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                            coefs(i, :), array_out(i, :))
    end do
    
end subroutine get_forcesurf_st3d