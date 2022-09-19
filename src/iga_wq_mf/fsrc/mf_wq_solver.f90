! ==========================
! module: Matrix free solver
! author: Joaquin Cornejo
! ==========================

subroutine mf_iga_steady_heat_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIterPCG, threshold, & 
                                method, directsol, x, RelRes, RelError)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! IN CSR FORMAT
                        
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    character(len = 10) :: method
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b, directsol
    dimension :: b(nr_u*nr_v*nr_w), directsol(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_u*nr_v*nr_w), RelRes(nbIterPCG+1), RelError(nbIterPCG+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, normb
    double precision, allocatable, dimension(:) :: r, p, Ap, dummy, z
    integer :: nr_total, iter, i

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: k_u, k_v, k_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set initial variables
    nr_total = nr_u*nr_v*nr_w
    allocate(r(nr_total), p(nr_total), Ap(nr_total), dummy(nr_total), z(nr_total))
    x = 0.d0; r = b
    RelRes = 0.d0; RelRes(1) = 1.d0
    RelError = 0.d0; RelError(1) = 1.d0
    normb = maxval(abs(r))
    if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 

    if (method.eq.'WP') then 
        if (nbIterPCG.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            p = r
            rsold = dot_product(r, r)
            
            do iter = 1, nbIterPCG
                call mf_iga_get_ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                    data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, p, Ap)
                
                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                RelRes(iter+1) = maxval(abs(r))/normb
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.threshold) exit
            
                rsnew = dot_product(r, r)
                p = r + rsnew/rsold * p
                rsold = rsnew
            end do
        end if

    else 
        ! Custumize method
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        k_u = 1.d0; k_v = 1.d0; k_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! Compute mean
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(1, 1, :), k_u) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(2, 2, :), k_v) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(3, 3, :), k_w) 
        
        end if

        ! --------------------
        ! Eigen decomposition
        ! --------------------
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
        allocate(data_W_u(nnz_u, 2), Kdiag_u(nr_u), Mdiag_u(nr_u))
        do i = 1, nnz_u
            data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
        end do
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 2), (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)
        deallocate(data_W_u)
        
        allocate(data_W_v(nnz_v, 2), Kdiag_v(nr_v), Mdiag_v(nr_v))
        do i = 1, nnz_v
            data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
        end do
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 2), (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)    
        deallocate(data_W_v)

        allocate(data_W_w(nnz_w, 2), Kdiag_w(nr_w), Mdiag_w(nr_w))
        do i = 1, nnz_w
            data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
        end do
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 2), (/0, 0/), D_w, U_w, Kdiag_w, Mdiag_w)  
        deallocate(data_W_w)
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        ! Find diagonal of eigen values
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(nr_total))
        I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_u, k_v, k_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, k_u, k_v, k_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call iga_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, &
                                nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, Dphysical)
        end if
        
        if (nbIterPCG.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            dummy = r 
            if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
            end if
            call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
            if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                call fd_sqr_scaling(nr_total, Dparametric, Dphysical, z) 
            end if
            p = z
            rsold = dot_product(r, z)

            do iter = 1, nbIterPCG
                call mf_iga_get_ku_3D(coefs, nc_total, nr_u, nc_u, &
                            nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, p, Ap)

                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                RelRes(iter+1) = maxval(abs(r))/normb
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.threshold) exit       

                dummy = r
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, z) 
                end if
                rsnew = dot_product(r, z)

                p = z + rsnew/rsold * p
                rsold = rsnew
            end do
        end if
    end if

end subroutine mf_iga_steady_heat_3d

subroutine mf_wq_steady_heat_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            b, nbIterPCG, threshold, method, directsol, x, RelRes, RelError)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    character(len=10), intent(in) :: method
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b, directsol
    dimension :: b(nr_u*nr_v*nr_w), directsol(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_u*nr_v*nr_w), RelRes(nbIterPCG+1), RelError(nbIterPCG+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, normb
    double precision, allocatable, dimension(:) :: r, rhat, p, Ap, s, As, dummy, ptilde, Aptilde, stilde, Astilde
    integer :: nr_total, iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: k_u, k_v, k_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set initial values
    nr_total = nr_u*nr_v*nr_w
    allocate(r(nr_total), rhat(nr_total), p(nr_total), Ap(nr_total), As(nr_total), s(nr_total), &
            dummy(nr_total), ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total))
    x = 0.d0; r = b; rhat = r; p = r
    RelRes = 0.d0; RelRes(1) = 1.d0
    RelError = 0.d0; RelError(1) = 1.d0
    rsold = dot_product(r, rhat); normb = maxval(abs(r))
    if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 

    if (method.eq.'WP') then 
        if (nbIterPCG.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            do iter = 1, nbIterPCG
                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, Ap)
                alpha = rsold/dot_product(Ap, rhat)
                s = r - alpha*Ap

                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, s, As)
                omega = dot_product(As, s)/dot_product(As, As)
                x = x + alpha*p + omega*s
                r = s - omega*As
    
                RelRes(iter+1) = maxval(abs(r))/normb
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.threshold) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew
            end do
        end if

    else  
        ! Custumize method
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0
        k_u = 1.d0; k_v = 1.d0; k_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! Compute mean
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(1, 1, :), k_u) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(2, 2, :), k_v) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(3, 3, :), k_w) 
        
        end if

        ! --------------------
        ! Eigen decomposition
        ! -------------------- 
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w), Kdiag_u(nr_u), Mdiag_u(nr_u))

        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)

        allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)    

        allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), (/0, 0/), D_w, U_w, Kdiag_w, Mdiag_w)   
        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        ! Find diagonal of eigen values
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(nr_total))
        I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_u, k_v, k_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, k_u, k_v, k_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        end if

        if (nbIterPCG.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            do iter = 1, nbIterPCG

                dummy = p
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if 
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, ptilde)  
                end if

                call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)

                alpha = rsold/dot_product(Aptilde, rhat)
                s = r - alpha*Aptilde

                dummy = s
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy)  
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
                if ((method.eq.'TDS').or.(method.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, stilde) 
                end if

                call mf_wq_get_Ku_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)

                omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
                x = x + alpha*ptilde + omega*stilde
                r = s - omega*Astilde    

                RelRes(iter+1) = maxval(abs(r))/normb
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.threshold) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                rsold = rsnew
            end do
        end if
    end if

end subroutine mf_wq_steady_heat_3d

subroutine fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, pardiag, phydiag, method, &
                    ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, array_in, array_out)
    !! Solves MM.s = r (MM is the preconditioner) in steady (or transient) heat 3D  with substitution method, 
    !! where MM is an approximation of K or K+C.
    !! IN CSR FORMAT

    implicit none
    ! Input / output data 
    !--------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w, ndof
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag, pardiag, phydiag
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigendiag(nr_total), pardiag(nr_total), phydiag(nr_total)

    integer, intent(in) :: indi_L, indj_L, indi_LT, indj_LT
    dimension :: indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof)
    double precision :: L, LT
    dimension :: L(ndof), LT(ndof)
    character(len=10), intent(in) :: method

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data
    ! ----------
    double precision :: array_temp_0, array_temp_1
    dimension :: array_temp_0(nr_total), array_temp_1(nr_total)

    ! Compute L'.array                 
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, array_temp_0)

    ! Do scaling depending on the method used
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, array_temp_0)
    end if

    ! By now, we test this approximation
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, array_temp_0, array_temp_1)

    ! Do scaling depending on the method used
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, array_temp_1)
    end if

    ! Compute L.array             
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, array_temp_1, array_out)

end subroutine fd_tshs_3d

subroutine mf_wq_transient_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            newmarkdt, b, nbIterPCG, threshold, method, x, residue)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    character(len=10), intent(in) :: method
    double precision, intent(in) :: newmarkdt
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x, residue
    dimension :: x(nr_u*nr_v*nr_w), residue(nbIterPCG+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, normb
    double precision, dimension(size(b)) :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, dummy
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen, Dtemp
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: csigma, k_u, k_v, k_w

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set variables
    x = 0.d0; r = b; rhat = r; p = r
    rsold = dot_product(r, rhat); normb = maxval(abs(r))
    residue = 0.d0; residue(1) = 1.d0
    if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 

    ! Custumize method  
    k_u = 1.d0; k_v = 1.d0; k_w = 1.d0; csigma = 1.d0
    if ((method.eq.'JMC').or.(method.eq.'JMS')) then 
        ! Compute mean 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(1, 1, :), k_u) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(2, 2, :), k_v) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(3, 3, :), k_w) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Ccoefs(:), csigma)
    end if

    ! --------------------
    ! Eigen decomposition
    ! -------------------- 
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w), Kdiag_u(nr_u), Mdiag_u(nr_u))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), (/0, 0/), D_w, U_w, Kdiag_w, Mdiag_w)   
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(size(b)))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_u, k_v, k_w, Deigen)
    Deigen = csigma + newmarkdt*Deigen
    deallocate(I_u, I_v, I_w)

    ! Scaling using diagonals
    allocate(Dparametric(size(b)), Dphysical(size(b)))
    Dparametric = 1.d0; Dphysical = 1.d0
    if (method.eq.'JMS') then
        ! Find diagonal of preconditioner
        allocate(Dtemp(size(b)))
        Dtemp = 0.d0
        call kron_product_3vec(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, csigma)
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, k_u, k_v, k_w, Dparametric)
        Dparametric = newmarkdt*Dparametric + Dtemp 

        ! Find diagonal of real matrix (K + C in this case)
        call wq_find_capacity_diagonal_3D(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dtemp)
        call wq_find_conductivity_diagonal_3D(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        Dphysical = newmarkdt*Dphysical + Dtemp
    end if

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    do iter = 1, nbIterPCG

        dummy = p
        if (method.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, dummy) 
        call fd_steady_heat_3d(size(b), nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
        if (method.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, ptilde)  

        call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, ptilde, Aptilde)

        alpha = rsold/dot_product(Aptilde, rhat)
        s = r - alpha*Aptilde

        dummy = s
        if (method.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, dummy) 
        call fd_steady_heat_3d(size(b), nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
        if (method.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, stilde)  

        call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, stilde, Astilde)

        omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
        x = x + alpha*ptilde + omega*stilde
        r = s - omega*Astilde    
        
        residue(iter+1) = maxval(abs(r))/normb
        if (residue(iter+1).le.threshold) exit

        rsnew = dot_product(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        rsold = rsnew

    end do

end subroutine mf_wq_transient_linear_3d

subroutine mf_wq_transient_nonlinear_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, nbpts, table_cond, table_cap, &
                        newmark, invJJ, detJJ, sizeF, time_list, FF, temperature)
    
    use heat_transfer
    implicit none 
    ! Input / output data
    ! -------------------
    character(len=10) :: method = 'FDC'
    double precision, parameter :: threshold = 1.d-8
    integer, parameter :: nbIterNL = 20, nbIterPCG = 100, d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF, nbpts
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: table_cond, table_cap, newmark
    dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2)
    double precision, intent(in) :: time_list, invJJ, detJJ, FF
    dimension ::    time_list(sizeF), invJJ(d, d, nc_total), &
                    detJJ(nc_total), FF(nr_u*nr_v*nr_w, sizeF)
    
    double precision, intent(out) :: temperature
    dimension :: temperature(nr_u*nr_v*nr_w, sizeF)

    ! Local data
    ! ----------  
    integer :: i, j, nr_total
    double precision, allocatable, dimension(:) :: TTn0, TTn1i0, VVn0, VVn1, CdT, KT, Fstep, KTCdT, ddFF, ddVV, TTn1
    double precision :: TTinterp, KK, CC, Kcoefs, CCoefs
    dimension :: TTinterp(nc_total), KK(nc_total), CC(nc_total), Kcoefs(dimen, dimen, nc_total), Ccoefs(nc_total)

    double precision :: resPCG, resNL, dt
    dimension :: resPCG(nbIterPCG+1)

    ! --------------------------------------------
    ! SOLVE
    ! -------------------------------------------- 
    nr_total = nr_u*nr_v*nr_w
    allocate(TTn0(nr_total), TTn1i0(nr_total), VVn0(nr_total), VVn1(nr_total), CdT(nr_total), KT(nr_total), &
            Fstep(nr_total), KTCdT(nr_total), ddFF(nr_total), ddVV(nr_total), TTn1(nr_total))

    ! Initialize
    temperature = 0.d0; VVn0 = 0.d0
    
    do i = 2, sizeF
        ! Get delta time
        dt = time_list(i) - time_list(i-1)

        ! Get values of last simulation 
        TTn0 = temperature(:, i-1)

        ! Prediction of new step
        TTn1 = TTn0 + dt*(1-newmark)*VVn0
        TTn1i0 = TTn1; VVn1 = 0.d0

        ! Get force of new step
        Fstep = FF(:, i)
        
        ! Solver Newton-Raphson
        print*, 'Step: ', i - 1
        do j = 1, nbIterNL

            ! Compute temperature (at each quadrature point) 
            call interpolate_temperature_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, TTn1, TTinterp)

            ! Interpolate capacity and conductivity at each quadrature point 
            call compute_heat_properties(nbpts, table_cond, table_cap, nc_total, TTinterp, KK, CC)

            ! Compute coefficients to compute tangent matrix
            call compute_heat_coefficients(nc_total, KK, CC, invJJ, detJJ, Kcoefs, Ccoefs)
            
            ! Compute Fint = C dT + K T 
            call mf_wq_get_cu_3d_csr(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, VVn1, CdT)

            call mf_wq_get_ku_3d_csr(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, TTn1, KT)
            KTCdT = KT + CdT

            ! Compute residue
            ddFF = Fstep - KTCdT
            resNL = sqrt(dot_product(ddFF, ddFF))
            print*, " with Raphson error: ", resNL
            if (isnan(resNL)) stop 'Error is NAN'
            if (resNL.le.1.d-6) exit

            ! Solve by iterations 
            call mf_wq_transient_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    newmark*dt, ddFF, nbIterPCG, threshold, method, ddVV, resPCG)

            ! Update values
            VVn1 = VVn1 + ddVV
            TTn1 = TTn1i0 + newmark*dt*VVn1
                
        end do
        
        ! Set values
        temperature(:, i) = TTn1
        VVn0 = VVn1
                
    end do

end subroutine mf_wq_transient_nonlinear_3d

! ============================
! 'COMPLETE' METHODS
! ============================

! STEADY HEAT TRANSFER CASE:
! --------------------------

subroutine mf_wq_get_Au_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_L, indj_L, indi_LT, indj_LT, L, LT, &
                                ndof, array_in, array_out)
    !! Computes Ann.u in steady heat 3D with substitution method, where 
    !! But A is given as [Ann, And; Adn, Add]. So Ann u =  L A L' u, where L is a zeros and ones matrix.
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    integer, intent(in) :: indi_L, indj_L, indi_LT, indj_LT
    dimension :: indi_L(ndof+1), indj_L(ndof), indi_LT(nr_u*nr_v*nr_w+1), indj_LT(ndof)
    double precision :: L, LT
    dimension :: L(ndof), LT(ndof)

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data 
    ! ----------
    integer :: nr_total
    double precision, allocatable, dimension(:) :: array_temp_0, array_temp_1

    ! Set total number of rows
    nr_total = nr_u*nr_v*nr_w

    ! Compute L'.u    
    allocate(array_temp_0(nr_total))               
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, array_temp_0)

    ! Compute K.u
    allocate(array_temp_1(nr_total))
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_temp_0, array_temp_1)
    deallocate(array_temp_0)

    ! Compute L.u                   
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, array_temp_1, array_out)

end subroutine mf_wq_get_Au_shs_3d

subroutine mf_wq_solve_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIterPCG, threshold, method, x, residue)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve K x = F, with B.x = g (Dirichlet condition). Using substitution method, this is
    !! Knn xn = Fn - Knd xd and xd = g 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: table, dod
    dimension :: table(d, 2), dod(ndod)
    character(len=10), intent(in) :: method
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, f, g
    dimension :: f(nr_u*nr_v*nr_w), g(ndod)
    
    double precision, intent(out) :: x, residue
    dimension :: x(nr_u*nr_v*nr_w), residue(nbIterPCG+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, normb
    double precision, allocatable, dimension(:) :: r, rhat, p, Ap, s, As, ptilde, Aptilde, stilde, Astilde, xn, Kx, b
    integer :: nr_total, iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: k_u, k_v, k_w

    ! Block L
    integer :: ndof
    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set initial variables
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    nr_total = nr_u*nr_v*nr_w
    ndof = nr_total - ndod
    allocate(r(ndof), rhat(ndof), p(ndof), Ap(ndof), As(ndof), s(ndof), ptilde(ndof), Aptilde(ndof), &
            Astilde(ndof), stilde(ndof), xn(ndof), Kx(nr_total), b(nr_total))

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! Compute K.x where x = [0, xd], then K.x = [Knd xd, Kdd xd]
    x = 0.d0; x(dod) = g
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Kx)
    b = f - Kx ! This is the real b in Ax = b equation
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, b, r)

    ! Set variables
    xn = 0.d0; rhat = r; p = r
    rsold = dot_product(r, rhat); normb = maxval(abs(r))
    residue = 0.d0; residue(1) = 1.d0
    if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 
    
    if (method.eq.'WP') then 
        ! ----------------------------
        ! Conjugate Gradient algorithm
        ! ----------------------------
        do iter = 1, nbIterPCG
            call mf_wq_get_Au_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, p, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call mf_wq_get_Au_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, s, As)
            omega = dot_product(As, s)/dot_product(As, As)
            xn = xn + alpha*p + omega*s
            r = s - omega*As

            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
            x(dod) = g
            
            residue(iter+1) = maxval(abs(r))/normb
            if (residue(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    else  
        ! Custumize method
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0 
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0

        k_u = 1.d0; k_v = 1.d0; k_w = 1.d0
        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! Compute mean 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(1, 1, :), k_u) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(2, 2, :), k_v) 
            call compute_mean_3d(nc_u, nc_v, nc_w, coefs(3, 3, :), k_w) 
    
        end if

        ! --------------------
        ! Eigen decomposition
        ! --------------------
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w), Kdiag_u(nr_u), Mdiag_u(nr_u))

        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), table(1, :), D_u, U_u, Kdiag_u, Mdiag_u)

        allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), table(2, :), D_v, U_v, Kdiag_v, Mdiag_v)    

        allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), table(3, :), D_w, U_w, Kdiag_w, Mdiag_w) 

        deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

        ! Find diagonal of eigen values
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(nr_total))
        I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_u, k_v, k_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, k_u, k_v, k_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        end if

        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        do iter = 1, nbIterPCG

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, p, ptilde)

            call mf_wq_get_Au_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, s, stilde)

            call mf_wq_get_Au_shs_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            xn = xn + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
            x(dod) = g

            residue(iter+1) = dot_product(r, r)/normb
            if (residue(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew

        end do

    end if

end subroutine mf_wq_solve_shs_3d

! TRANSIENT HEAT TRANSFER CASE: 
! -----------------------------

subroutine mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_L, indj_L, indi_LT, indj_LT, L, LT, &
                                ndof, newmarkdt, array_in, array_out)
    !! Computes Ann.u in transient heat 3D with substitution method, where 
    !! But A is given as [Ann, And; Adn, Add]. So Ann u =  L A L' u, where L is a zeros and ones matrix.
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)

    integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) :: indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    integer, intent(in) :: indi_L, indj_L, indi_LT, indj_LT
    dimension :: indi_L(ndof+1), indj_L(ndof), indi_LT(nr_u*nr_v*nr_w+1), indj_LT(ndof)
    double precision :: L, LT, newmarkdt
    dimension :: L(ndof), LT(ndof)

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data 
    ! ---------- 
    integer :: nr_total
    double precision, allocatable, dimension(:) :: array_temp_0, array_temp_1

    ! Set total number of rows
    nr_total = nr_u*nr_v*nr_w

    ! Compute L'.array   
    allocate(array_temp_0(nr_total))                
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, array_temp_0)

    ! Compute w = (C + alpha dt K).v
    allocate(array_temp_1(nr_total))
    call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, array_temp_0, array_temp_1)
    deallocate(array_temp_0)

    ! Compute L.w                 
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, array_temp_1, array_out)

end subroutine mf_wq_get_Au_ths_3d

subroutine mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            newmarkdt, table, ndod, dod, f, g, nbIterPCG, threshold, method, x, residue)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    character(len=10), intent(in) :: method
    double precision, intent(in) :: newmarkdt
    integer, intent(in) :: table, dod
    dimension :: table(d, 2), dod(ndod)
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, f, g
    dimension :: f(nr_u*nr_v*nr_w), g(ndod)
    
    double precision, intent(out) :: x, residue
    dimension :: x(nr_u*nr_v*nr_w), residue(nbIterPCG+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, normb
    double precision, allocatable, dimension(:) :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, xn, Ax, b
    integer :: nr_total, iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen, Dtemp
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: csigma, k_u, k_v, k_w

    ! Block L
    integer :: ndof
    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Set initial variables
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    nr_total = nr_u*nr_v*nr_w
    ndof = nr_total - ndod
    allocate(r(ndof), rhat(ndof), p(ndof), s(ndof), ptilde(ndof), Aptilde(ndof), &
            Astilde(ndof), stilde(ndof), xn(ndof), Ax(nr_total), b(nr_total))

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! Compute A.x where x = [0, xd] and A = (C + alpha dt K)
    x = 0.d0; x(dod) = g
    call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, 1.d0, newmarkdt, x, Ax)

    b = f - Ax ! This is the real b in Ax = b equation
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, b, r)

    ! Set variables
    xn = 0.d0; rhat = r; p = r
    rsold = dot_product(r, rhat); normb = maxval(abs(r))
    residue = 0.d0; residue(1) = 1.d0
    if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 

    ! Custumize method  
    k_u = 1.d0; k_v = 1.d0; k_w = 1.d0; csigma = 1.d0
    if ((method.eq.'JMC').or.(method.eq.'JMS')) then 
        ! Compute mean 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(1, 1, :), k_u) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(2, 2, :), k_v) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(3, 3, :), k_w) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Ccoefs(:), csigma)
    end if

    ! --------------------
    ! Eigen decomposition
    ! -------------------- 
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w), Kdiag_u(nr_u), Mdiag_u(nr_u))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), table(1, :), D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), table(2, :), D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), table(3, :), D_w, U_w, Kdiag_w, Mdiag_w)   
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(nr_total))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_u, k_v, k_w, Deigen)
    Deigen = csigma + newmarkdt*Deigen
    deallocate(I_u, I_v, I_w)

    ! Scaling using diagonals
    allocate(Dparametric(nr_total), Dphysical(nr_total))
    Dparametric = 1.d0; Dphysical = 1.d0
    if (method.eq.'JMS') then
        ! Find diagonal of preconditioner
        allocate(Dtemp(nr_total))
        Dtemp = 0.d0
        call kron_product_3vec(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, csigma)
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, k_u, k_v, k_w, Dparametric)
        Dparametric = newmarkdt*Dparametric + Dtemp 

        ! Find diagonal of real matrix (C + alpha dt K in this case)
        call wq_find_capacity_diagonal_3D(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dtemp)
        call wq_find_conductivity_diagonal_3D(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        Dphysical = newmarkdt*Dphysical + Dtemp
    end if

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    do iter = 1, nbIterPCG

        call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                        ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, p, ptilde)

        call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, newmarkdt, ptilde, Aptilde)

        alpha = rsold/dot_product(Aptilde, rhat)
        s = r - alpha*Aptilde

        call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                        ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, s, stilde)

        call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, newmarkdt, stilde, Astilde)

        omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
        xn = xn + alpha*ptilde + omega*stilde
        r = s - omega*Astilde    
        
        call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
        x(dod) = g
        residue(iter+1) = maxval(abs(r))/normb
        if (residue(iter+1).le.threshold) exit

        rsnew = dot_product(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        rsold = rsnew

    end do

end subroutine mf_wq_solve_ths_linear_3d

subroutine mf_wq_ths_nonlinear_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, nbpts, table_cond, table_cap, &
                        newmark, table_dir, ndod, dod, invJJ, detJJ, sizeF, time_list, FF, GG, temperature)
    
    use heat_transfer
    implicit none 
    ! Input / output data
    ! -------------------
    character(len=10) :: method = 'FDC'
    double precision, parameter :: threshold = 1.d-8
    integer, parameter :: nbIterNL = 20, nbIterPCG = 100, d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF, nbpts
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: table_cond, table_cap, newmark
    dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2)
    integer, intent(in) :: ndod, table_dir, dod
    dimension :: table_dir(d, 2), dod(ndod)
    double precision, intent(in) :: time_list, invJJ, detJJ, FF, GG
    dimension ::    time_list(sizeF), invJJ(d, d, nc_total), &
                    detJJ(nc_total), FF(nr_u*nr_v*nr_w, sizeF), GG(ndod, sizeF)
    
    double precision, intent(out) :: temperature
    dimension :: temperature(nr_u*nr_v*nr_w, sizeF)

    ! Local data
    ! ----------  
    integer :: nr_total, ndof, i, j
    double precision, allocatable, dimension(:) :: TTn0, TTn1i0, VVn0, VVn1, CdT, KT, Fstep, &
                                                    KTCdT, ddFF, ddVV, TTn1, TTinterp, GGtmp, ddGG
    double precision :: resPCG, resNL, dt, dt2, factor
    dimension :: resPCG(nbIterPCG+1)

    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT

    double precision :: KK, CC, Kcoefs, CCoefs
    dimension :: KK(nc_total), CC(nc_total), Kcoefs(dimen, dimen, nc_total), Ccoefs(nc_total)

    ! Set initial values
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    nr_total = nr_u*nr_v*nr_w
    ndof = nr_total - ndod
    allocate(GGtmp(ndod), ddGG(ndod))
    GGtmp = 0.d0; ddGG = 0.d0

    ! Compute initial velocity from boundary conditions
    if (sizeF.eq.2) then 
        dt = time_list(2) - time_list(1)
        ddGG = 1.d0/dt*(GG(:, 2) - GG(:, 1))
    else if (sizeF.gt.2) then
        dt = time_list(2) - time_list(1)
        dt2 = time_list(3) - time_list(1)
        factor = dt2/dt
        ddGG = 1.d0/(dt*(factor-factor**2))*(GG(:, 3) - (factor**2)*GG(:, 2) - (1 - factor**2)*GG(:, 1))
    else
        stop 'This solver needs at least 2 steps'
    end if

    ! Save data
    allocate(VVn0(nr_total))
    VVn0 = 0.d0; VVn0(dod) = ddGG
    deallocate(ddGG)

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! --------------------------------------------
    ! SOLVE
    ! -------------------------------------------- 
    ! Initialize
    allocate(TTn0(nr_total), TTn1i0(nr_total), VVn1(nr_total), CdT(nr_total), KT(nr_total), Fstep(nr_total), &
            KTCdT(nr_total), ddFF(nr_total), ddVV(nr_total), TTn1(nr_total), TTinterp(nc_total))
    temperature = 0.d0; 
    
    do i = 2, sizeF
        ! Get delta time
        dt = time_list(i) - time_list(i-1)

        ! Get values of last simulation 
        TTn0 = temperature(:, i-1)

        ! Prediction of new step
        TTn1 = TTn0 + dt*(1-newmark)*VVn0; TTn1(dod) = GG(:, i)
        TTn1i0 = TTn1; VVn1 = 0.d0
        VVn1(dod) = 1.d0/newmark*(1.0d0/dt*(GG(:, i) - GG(:, i-1)) - (1 - newmark)*VVn0(dod))

        ! Get force of new step
        Fstep = FF(:, i)
        
        ! Solver Newton-Raphson
        print*, 'Step: ', i - 1
        do j = 1, nbIterNL

            ! Compute temperature (at each quadrature point) 
            call interpolate_temperature_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, TTn1, TTinterp)

            ! Interpolate capacity and conductivity at each quadrature point 
            call compute_heat_properties(nbpts, table_cond, table_cap, nc_total, TTinterp, KK, CC)

            ! Compute coefficients to compute tangent matrix
            call compute_heat_coefficients(nc_total, KK, CC, invJJ, detJJ, Kcoefs, Ccoefs)
            
            ! Compute Fint = C dT + K T 
            call mf_wq_get_cu_3d_csr(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, VVn1, CdT)

            call mf_wq_get_ku_3d_csr(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, TTn1, KT)
            KTCdT = KT + CdT

            ! Compute residue
            ddFF = Fstep - KTCdT
            call clean_dirichlet_1dim(nr_total, ddFF, ndod, dod)
            resNL = sqrt(dot_product(ddFF, ddFF))
            print*, " with Raphson error: ", resNL
            if (isnan(resNL)) stop 'Error is NAN'
            if (resNL.le.1.d-6) exit

            ! Solve by iterations 
            call mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    newmark*dt, table_dir, ndod, dod, ddFF, GGtmp, nbIterPCG, threshold, method, ddVV, resPCG)

            ! Update values
            VVn1 = VVn1 + ddVV
            TTn1 = TTn1i0 + newmark*dt*VVn1
                
        end do
        
        ! Set values
        temperature(:, i) = TTn1
        VVn0 = VVn1
                
    end do

end subroutine mf_wq_ths_nonlinear_3d