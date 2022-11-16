! ==========================
! module: Matrix free solver
! author: Joaquin Cornejo
! ==========================

subroutine mf_wq_interpolate_cp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    
    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, Deigen
    double precision :: U_u, U_v, U_w
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w)

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
    ! Eigen decomposition
    allocate(   Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w), &
                Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    call eigen_decomposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, (/1.d0, 1.d0, 1.d0/), .false., &
                                U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
    deallocate( Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, &
                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Set material and solver
    allocate(mat, solv)
    call setupCcoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 1

    call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

end subroutine mf_wq_interpolate_cp_3d

subroutine mf_iga_interpolate_cp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIterPCG, threshold, x, resPCG)
    !! Preconditioned conjugate gradient to solve interpolation problem
    !! IN CSR FORMAT
    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, Deigen
    double precision :: U_u, U_v, U_w
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w)
    
    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call iga2wq3d(nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, indj_u, indj_v, indj_w, &
                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, data_W_u, data_W_v, data_W_w)
    
    ! Eigen decomposition
    allocate(   Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w), &
                Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    call eigen_decomposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, (/1.d0, 1.d0, 1.d0/), .false., &
                                U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
    deallocate( Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, &
                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
    
    ! Set material and solver
    allocate(mat, solv)
    call setupCcoefs(mat, nc_total, coefs)
    solv%matrixfreetype = 1

    call PCG(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)
            
end subroutine mf_iga_interpolate_cp_3d

subroutine mf_wq_steady_heat_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    integer :: i
    double precision :: kmean(3), kappa

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_total)
    double precision, dimension(:), allocatable :: Dphysical, Dparametric

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    allocate(mat, solv)
    call setupKcoefs(mat, 3, nc_total, coefs)
    solv%matrixfreetype = 2

    if (methodPCG.eq.'WP') then ! CG algorithm
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)
        
    else  

        ! Customize method
        allocate(Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        kmean = 1.d0

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'TDC')) then 
            
            ! Separation of variables
            do i = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((methodPCG.eq.'JMS').or.(methodPCG.eq.'JMC')) then 

            ! Jacobian mean
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(1, 1, :), kmean(1)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(2, 2, :), kmean(2)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(3, 3, :), kmean(3)) 
        
        end if

        ! Condition number 
        if ((methodPCG.eq.'C').or.(methodPCG.eq.'JMC')) then
            
            call compute_steady_condition_number(nc_total, coefs, kmean, kappa)
            print*, 'Method: ', methodPCG, ', condition number: ', kappa

        end if

        ! Eigen decomposition
        allocate(Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
        call eigen_decomposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, kmean, .true., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        call setup_eigendiag(solv, nr_total, Deigen)

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then

            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, kmean(1), kmean(2), kmean(3), Dparametric)
            
            ! Find diagonal of real matrix
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)

            call setup_FDscaling(solv, nr_total, Dparametric, Dphysical)
        end if
        deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

        ! Set solver
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

    end if

end subroutine mf_wq_steady_heat_3d

subroutine mf_iga_steady_heat_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, &
                                nbIterPCG, threshold, methodPCG, x, resPCG)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python)
    !! IN CSR FORMAT

    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(3, 3, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2) 
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    character(len = 10) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    integer :: i
    double precision :: kmean(3)

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_total)
    double precision, dimension(:), allocatable :: Dphysical, Dparametric

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call iga2wq3d(nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, indj_u, indj_v, indj_w, &
                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, data_W_u, data_W_v, data_W_w)

    allocate(mat, solv)
    call setupKcoefs(mat, 3, nc_total, coefs)
    solv%matrixfreetype = 2

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call CG(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

    else 
        ! Customize method
        allocate(Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        kmean = 1.d0

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'TDC')) then 
            
            ! Separation of variables
            do i = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((methodPCG.eq.'JMS').or.(methodPCG.eq.'JMC')) then 

            ! Jacobian mean
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(1, 1, :), kmean(1)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(2, 2, :), kmean(2)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(3, 3, :), kmean(3)) 
        
        end if

        ! Eigen decomposition
        allocate(Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
        call eigen_decomposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, kmean, .true., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
        call setup_eigendiag(solv, nr_total, Deigen)

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then

            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, kmean(1), kmean(2), kmean(3), Dparametric)
            
            ! Find diagonal of real matrix
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)

            call setup_FDscaling(solv, nr_total, Dparametric, Dphysical)
        end if
        deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

        ! Set solver
        call PCG(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

    end if

end subroutine mf_iga_steady_heat_3d

subroutine mf_wq_transient_linear_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                b, thetadt, nbIterPCG, threshold, methodPCG, x, resPCG)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(3, 3, nc_total), Ccoefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG    
    double precision, intent(in) :: thetadt, threshold, b
    dimension :: b(nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! ----------
    type(thermomat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision :: kmean(3), cmean, kappa

    ! Fast diagonalization
    double precision, dimension(:), allocatable ::  Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, &
                                                    Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_total)
    double precision, dimension(:), allocatable :: Dphysical, Dparametric, Dtemp

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    allocate(mat, solv)
    call setupCcoefs(mat, nc_total, Ccoefs)
    call setupKcoefs(mat, 3, nc_total, Kcoefs)
    mat%scalars = (/1.d0, thetadt/)
    solv%matrixfreetype = 3

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)
        
    else  

        ! Customize method  
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        kmean = 1.d0; cmean = 1.d0

        if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'JMS')) then 
            
            ! Jacobian mean 
            call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(1, 1, :), kmean(1)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(2, 2, :), kmean(2)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(3, 3, :), kmean(3)) 
            call compute_mean_3d(nc_u, nc_v, nc_w, Ccoefs(:), cmean)
        
        end if

        ! Eigen decomposition
        allocate(Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w))            
        call eigen_decomposition_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, kmean, .true., &
                            U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        Deigen = cmean + thetadt*Deigen
        call setup_eigendiag(solv, nr_total, Deigen)

        if (methodPCG.eq.'JMS') then
            allocate(Dparametric(nr_total), Dphysical(nr_total), Dtemp(nr_total))

            ! Find diagonal of preconditioner
            Dtemp = 0.d0
            call kronvec3d(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, cmean)
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, kmean(1), kmean(2), kmean(3), Dparametric)
            Dparametric = Dtemp + thetadt*Dparametric

            ! Find diagonal of real matrix 
            Dtemp = 0.d0
            call wq_find_capacity_diagonal_3D(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dtemp)
            call wq_find_conductivity_diagonal_3D(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
            Dphysical = Dtemp + thetadt*Dphysical

            call setup_FDscaling(solv, nr_total, Dparametric, Dphysical)
        end if
        deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

        ! Condition number P^-1 A
        call compute_transient_condition_number(nc_total, Kcoefs, Ccoefs, kmean, cmean, kappa)
        if (methodPCG.eq.'JMS') then
            Dtemp = Dparametric/Dphysical
            kappa = kappa*maxval(Dtemp)/minval(Dtemp)
        end if
        print*, 'Method: ', methodPCG, ' condition number: ', kappa

        ! Set solver
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

    end if

end subroutine mf_wq_transient_linear_3d

subroutine mf_wq_transient_nonlinear_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    !&
                                    nr_total_t, nr_u_t, nr_v_t, nr_w_t, nnz_u_t, nnz_v_t, nnz_w_t, & 
                                    indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                    data_B_u_t, data_B_v_t, data_B_w_t, data_W_u_t, data_W_v_t, data_W_w_t, &
                                    !&
                                    nbsteps, time_list, FF, ndod, dod, GG, nbpts, table_cond, table_cap, invJJ, detJJ, &
                                    theta, methodPCG, solution, resPCG)
    !! Time-integration scheme (with Newton-Raphson method) to solve non linear transient heat problems
    !! For the moment, only for isotropic materials.
    !! It assumes that initial temperature equals 0
    !! IN CSR FORMAT
    
    use heat_spmf
    use heat_solver

    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: thresholdPCG = 1.d-11, thresholdNL=1.d-11
    integer, parameter :: nbIterNL = 20, nbIterPCG = 100
    integer, intent(in) :: nr_total_t, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    integer :: nr_u_t, nr_v_t, nr_w_t, nnz_u_t, nnz_v_t, nnz_w_t
    integer, intent(in) :: indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t
    dimension ::    indi_u_t(nr_u_t+1), indj_u_t(nnz_u_t), &
                    indi_v_t(nr_v_t+1), indj_v_t(nnz_v_t), &
                    indi_w_t(nr_w_t+1), indj_w_t(nnz_w_t)
    double precision, intent(in) :: data_B_u_t, data_W_u_t, data_B_v_t, data_W_v_t, data_B_w_t, data_W_w_t
    dimension ::    data_B_u_t(nnz_u_t, 2), data_W_u_t(nnz_u_t, 4), &
                    data_B_v_t(nnz_v_t, 2), data_W_v_t(nnz_v_t, 4), &
                    data_B_w_t(nnz_w_t, 2), data_W_w_t(nnz_w_t, 4)

    integer, intent(in) :: nbsteps, ndod, nbpts
    double precision, intent(in) :: time_list, FF, GG
    dimension :: time_list(nbsteps), FF(nr_total_t, nbsteps), GG(ndod, nbsteps)
    integer, intent(in) :: dod
    dimension :: dod(ndod)

    double precision, intent(in) :: table_cond, table_cap, invJJ, detJJ
    dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2), invJJ(3, 3, nc_total), detJJ(nc_total)
    double precision, intent(in) :: theta
    character(len=10), intent(in) :: methodPCG

    double precision, intent(out) :: solution, resPCG
    dimension :: solution(nr_total_t, nbsteps), resPCG(nbIterPCG+3, nbIterNL*nbsteps)

    ! Local data
    ! ----------  
    type(thermomat), pointer :: mat
    integer :: i, j, k, ndof
    double precision :: dt, dt2, factor, resNL
    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT, ddGG, ddFFdof, ddVVdof
    double precision, dimension(nr_total_t) ::  TTn0, TTn1i0, VVn0, VVn1, CdT, KT, Fstep, KTCdT, ddFF, ddVV, TTn1
    double precision :: TTinterp, Kprop, Cprop
    dimension :: TTinterp(nc_total), Kprop(nc_total), Cprop(nc_total)
    double precision, target :: Ccoefs, Kcoefs
    dimension :: Ccoefs(nc_total), Kcoefs(3, 3, nc_total)

    integer :: indi_T_u_t, indi_T_v_t, indi_T_w_t, indj_T_u_t, indj_T_v_t, indj_T_w_t
    dimension ::    indi_T_u_t(nc_u+1), indi_T_v_t(nc_v+1), indi_T_w_t(nc_w+1), &
                    indj_T_u_t(nnz_u_t), indj_T_v_t(nnz_v_t), indj_T_w_t(nnz_w_t)
    double precision :: data_BT_u_t, data_BT_v_t, data_BT_w_t
    dimension :: data_BT_u_t(nnz_u_t, 2), data_BT_v_t(nnz_v_t, 2), data_BT_w_t(nnz_w_t, 2)
        
    if (nr_total_t.ne.nr_u_t*nr_v_t*nr_w_t) stop 'Size problem'
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    if (any(dod.gt.nr_total_t)) stop 'Indices must be less than nb of rows'

    call csr2csc(2, nr_u_t, nc_u, nnz_u_t, data_B_u_t, indj_u_t, indi_u_t, data_BT_u_t, indj_T_u_t, indi_T_u_t)
    call csr2csc(2, nr_v_t, nc_v, nnz_v_t, data_B_v_t, indj_v_t, indi_v_t, data_BT_v_t, indj_T_v_t, indi_T_v_t)
    call csr2csc(2, nr_w_t, nc_w, nnz_w_t, data_B_w_t, indj_w_t, indi_w_t, data_BT_w_t, indj_T_w_t, indi_T_w_t)

    solution(dod, :) = GG

    allocate(mat)
    Ccoefs = 0.d0; Kcoefs = 0.d0
    call setupCcoefs(mat, nc_total, Ccoefs)
    call setupKcoefs(mat, 3, nc_total, Kcoefs)
    
    ! Compute initial velocity from boundary conditions
    allocate(ddGG(ndod))
    if (nbsteps.eq.2) then 
        dt   = time_list(2) - time_list(1)
        ddGG = 1.d0/dt*(solution(dod, 2) - solution(dod, 1))
    else if (nbsteps.gt.2) then
        dt     = time_list(2) - time_list(1)
        dt2    = time_list(3) - time_list(1)
        factor = dt2/dt
        ddGG   = 1.d0/(dt*(factor - factor**2))*(solution(dod, 3) &
                - (factor**2)*solution(dod, 2) - (1 - factor**2)*solution(dod, 1))
    else
        stop 'This solver needs at least 2 steps'
    end if

    VVn0 = 0.d0; VVn0(dod) = ddGG
    deallocate(ddGG)

    ! Create block L
    ndof = nr_total_t - ndod
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total_t+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)

    ! Solve
    if (ndof.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    allocate(ddFFdof(ndof), ddVVdof(ndof))
    resPCG = 0.d0; k = 1
    do i = 2, nbsteps

        dt = time_list(i) - time_list(i-1)

        ! Get values of last step 
        TTn0 = solution(:, i-1)

        ! Get values of new step
        TTn1 = TTn0 + dt*(1 - theta)*VVn0; TTn1(dod) = solution(dod, i)
        TTn1i0 = TTn1; VVn1 = 0.d0
        VVn1(dod) = 1.d0/theta*(1.0d0/dt*(solution(dod, i) - solution(dod, i-1)) - (1 - theta)*VVn0(dod))
        Fstep = FF(:, i)
        
        print*, 'Step: ', i - 1
        do j = 1, nbIterNL ! Newton-Raphson algorithm
            ! Compute temperature and properties at each quadrature point
            call interp_temperature(nr_u_t, nc_u, nr_v_t, nc_v, nr_w_t, nc_w, nnz_u_t, nnz_v_t, nnz_w_t, &
                                        indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                        data_B_u_t, data_B_v_t, data_B_w_t, TTn1, TTinterp)
            call interp_isotropic_prop(nbpts, table_cap, nc_total, TTinterp, Cprop)
            call interp_isotropic_prop(nbpts, table_cond, nc_total, TTinterp, Kprop)
            
            ! Compute coefficients to compute tangent matrix
            call eval_isotropic_coefs(nc_total, Kprop, Cprop, invJJ, detJJ, Kcoefs, Ccoefs)

            ! Compute internal force = C dT + K T
            call mf_wq_get_cu_3d(mat, nr_total_t, nc_total, nr_u_t, nc_u, nr_v_t, nc_v, nr_w_t, nc_w, &
                            nnz_u_t, nnz_v_t, nnz_w_t, indi_T_u_t, indj_T_u_t, indi_T_v_t, indj_T_v_t, &
                            indi_T_w_t, indj_T_w_t, data_BT_u_t, data_BT_v_t, data_BT_w_t, &
                            indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                            data_W_u_t, data_W_v_t, data_W_w_t, VVn1, CdT)

            call mf_wq_get_ku_3d(mat, nr_total_t, nc_total, nr_u_t, nc_u, nr_v_t, nc_v, nr_w_t, nc_w, &
                            nnz_u_t, nnz_v_t, nnz_w_t, indi_T_u_t, indj_T_u_t, indi_T_v_t, indj_T_v_t, &
                            indi_T_w_t, indj_T_w_t, data_BT_u_t, data_BT_v_t, data_BT_w_t, &
                            indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                            data_W_u_t, data_W_v_t, data_W_w_t, TTn1, KT)

            KTCdT = KT + CdT

            ! Compute residue
            ddFF = Fstep - KTCdT
            call spMdotdV(ndof, nr_total_t, ndof, indi_L, indj_L, L, ddFF, ddFFdof)
            resNL = maxval(abs(ddFFdof))
            print*, 'Raphson error:', resNL
            if (resNL.le.thresholdNL) exit

            ! Iterative solver 
            resPCG(1, k) = dble(i-1); resPCG(2, k) = dble(j)
            call mf_wq_transient_linear_3d(Ccoefs, Kcoefs, ndof, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                        ddFFdof, theta*dt, nbIterPCG, thresholdPCG, methodPCG, ddVVdof, resPCG(3:, k))
            ! Update values
            call spMdotdV(nr_total_t, ndof, ndof, indi_LT, indj_LT, LT, ddVVdof, ddVV)
            VVn1 = VVn1 + ddVV
            TTn1 = TTn1i0 + theta*dt*VVn1
            k = k + 1
                
        end do
        
        solution(:, i) = TTn1
        VVn0 = VVn1

    end do

end subroutine mf_wq_transient_nonlinear_3d