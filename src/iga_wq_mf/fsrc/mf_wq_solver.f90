! ==========================
! module: Matrix free solver
! author: Joaquin Cornejo
! ==========================

subroutine mf_iga_steady_heat_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIterPCG, threshold, & 
                                methodPCG, directsol, x, RelRes, RelError)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python)
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

    character(len = 10) :: methodPCG
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

    if (methodPCG.eq.'WP') then 
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

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((methodPCG.eq.'JMS').or.(methodPCG.eq.'JMC')) then 
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

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
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
            if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
                call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
            end if
            call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
            if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
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
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, z)
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
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
                            b, nbIterPCG, threshold, methodPCG, directsol, x, RelRes, RelError)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
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

    character(len=10), intent(in) :: methodPCG
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

    if (methodPCG.eq.'WP') then 
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

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((methodPCG.eq.'JMS').or.(methodPCG.eq.'JMC')) then 
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

        if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
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
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
                end if 
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, ptilde)  
                end if

                call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)

                alpha = rsold/dot_product(Aptilde, rhat)
                s = r - alpha*Aptilde

                dummy = s
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
                    call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy)  
                end if
                call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
                if ((methodPCG.eq.'TDS').or.(methodPCG.eq.'JMS')) then
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

subroutine mf_wq_transient_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                newmarkdt, b, nbIterPCG, threshold, methodPCG, solution, residue)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (newmarkdt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

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
    character(len=10), intent(in) :: methodPCG
    double precision, intent(in) :: newmarkdt
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold, b
    dimension :: b(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: solution, residue
    dimension :: solution(nr_u*nr_v*nr_w), residue(nbIterPCG+1)

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
    double precision :: cm, k_um, k_vm, k_wm, kmeanvector(d), kappa

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
    solution = 0.d0; r = b; rhat = r; p = r
    rsold = dot_product(r, rhat); normb = maxval(abs(r))
    residue = 0.d0; residue(1) = 1.d0
    ! if (normb.lt.threshold) stop 'Force is almost zero, then it is a trivial solution' 

    ! Custumize method  
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    k_um = 1.d0; k_vm = 1.d0; k_wm = 1.d0; cm = 1.d0

    if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'JMS')) then 
        ! Compute mean 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(1, 1, :), k_um) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(2, 2, :), k_vm) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Kcoefs(3, 3, :), k_wm) 
        call compute_mean_3d(nc_u, nc_v, nc_w, Ccoefs(:), cm)
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
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w), Deigen(size(b)))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, k_um, k_vm, k_wm, Deigen)
    Deigen = cm + newmarkdt*Deigen
    deallocate(I_u, I_v, I_w)

    if (methodPCG.eq.'JMS') then
        allocate(Dparametric(size(b)), Dphysical(size(b)), Dtemp(size(b)))

        ! Find diagonal of preconditioner
        Dtemp = 0.d0
        call kronvec3d(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, cm)
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, k_um, k_vm, k_wm, Dparametric)
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

    ! Condition number P^-1 A
    kmeanvector = (/k_um, k_vm, k_wm/)
    call compute_transient_condition_number(nc_total, Kcoefs, Ccoefs, kmeanvector, cm, kappa)
    if (methodPCG.eq.'JMS') then
        Dtemp = Dparametric/Dphysical
        kappa = kappa*maxval(Dtemp)/minval(Dtemp)
    end if
    print*, 'Condition number: ', kappa

    ! -------------------------------------------
    ! Preconditioned Conjugate Gradient algorithm
    ! -------------------------------------------
    do iter = 1, nbIterPCG

        dummy = p
        if (methodPCG.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, dummy) 
        call fd_steady_heat_3d(size(b), nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
        if (methodPCG.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, ptilde)  

        call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, ptilde, Aptilde)

        alpha = rsold/dot_product(Aptilde, rhat)
        s = r - alpha*Aptilde

        dummy = s
        if (methodPCG.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, dummy) 
        call fd_steady_heat_3d(size(b), nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
        if (methodPCG.eq.'JMS') call fd_sqr_scaling(size(b), Dparametric, Dphysical, stilde)  

        call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, stilde, Astilde)

        omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
        solution = solution + alpha*ptilde + omega*stilde
        r = s - omega*Astilde    
        
        residue(iter+1) = maxval(abs(r))/normb
        if (residue(iter+1).le.threshold) exit

        rsnew = dot_product(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        rsold = rsnew

    end do

end subroutine mf_wq_transient_linear_3d

subroutine mf_wq_transient_nonlinear_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    nr_u_t, nc_u_t, nr_v_t, nc_v_t, nr_w_t, nc_w_t, nnz_u_t, nnz_v_t, nnz_w_t, & 
                                    indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                    data_B_u_t, data_B_v_t, data_B_w_t, data_W_u_t, data_W_v_t, data_W_w_t, &
                                    nbpts, table_cond, table_cap, newmark, invJJ, detJJ, sizeF, time_list, FF, &
                                    ndod, dod, GG, methodPCG, solution, resPCG)
    !! Time-integration scheme (with Newton-Raphson method) to solve non linear transient heat problems
    !! For the moment, only for isotropic materials.
    !! It assumes that initial temperature equals 0
    !! IN CSR FORMAT
    
    use heat_transfer
    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: thresholdPCG = 1.d-14, thresholdNL=1.d-14
    integer, parameter :: nbIterNL = 20, nbIterPCG = 100, d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    integer :: nr_u_t, nc_u_t, nr_v_t, nc_v_t, nr_w_t, nc_w_t, nnz_u_t, nnz_v_t, nnz_w_t
    integer, intent(in) :: indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t
    dimension ::    indi_u_t(nr_u_t+1), indj_u_t(nnz_u_t), &
                    indi_v_t(nr_v_t+1), indj_v_t(nnz_v_t), &
                    indi_w_t(nr_w_t+1), indj_w_t(nnz_w_t)
    double precision, intent(in) :: data_B_u_t, data_W_u_t, data_B_v_t, data_W_v_t, data_B_w_t, data_W_w_t
    dimension ::    data_B_u_t(nnz_u_t, 2), data_W_u_t(nnz_u_t, 4), &
                    data_B_v_t(nnz_v_t, 2), data_W_v_t(nnz_v_t, 4), &
                    data_B_w_t(nnz_w_t, 2), data_W_w_t(nnz_w_t, 4)

    integer, intent(in) :: nr_total, nbpts, sizeF, ndod, dod
    dimension :: dod(ndod)
    double precision, intent(in) :: table_cond, table_cap, newmark
    dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2)
    double precision, intent(in) :: time_list, invJJ, detJJ, FF, GG
    dimension :: time_list(sizeF), invJJ(d, d, nc_total), detJJ(nc_total), FF(nr_total, sizeF), GG(ndod, sizeF)
    character(len=10), intent(in) :: methodPCG

    double precision, intent(out) :: solution, resPCG
    dimension :: solution(nr_total, sizeF), resPCG(nbIterPCG+3, 8*sizeF)

    ! Local data
    ! ----------  
    integer :: i, j, c, ndof
    double precision :: resNL, dt, dt2, factor
    double precision, allocatable, dimension(:) ::  TTn0, TTn1i0, VVn0, VVn1, CdT, KT, Fstep, &
                                                    KTCdT, ddFF, ddVV, TTn1, ddGG, ddFFdof, ddVVdof
    double precision :: TTinterp, KK, CC, Kcoefs, CCoefs
    dimension ::    TTinterp(nc_total), KK(nc_total), CC(nc_total), Kcoefs(dimen, dimen, nc_total), &
                    Ccoefs(nc_total)
    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT

    ! Set initial values
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    ndof = nr_total - ndod
    allocate(ddGG(ndod))
    solution(dod, :) = GG
    resPCG = 0.0d0

    ! Compute initial velocity from boundary conditions
    if (sizeF.eq.2) then 
        dt = time_list(2) - time_list(1)
        ddGG = 1.d0/dt*(solution(dod, 2) - solution(dod, 1))
    else if (sizeF.gt.2) then
        dt = time_list(2) - time_list(1)
        dt2 = time_list(3) - time_list(1)
        factor = dt2/dt
        ddGG = 1.d0/(dt*(factor - factor**2))*(solution(dod, 3) - (factor**2)*solution(dod, 2) - (1 - factor**2)*solution(dod, 1))
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
    allocate(TTn0(nr_total), TTn1i0(nr_total), VVn1(nr_total), CdT(nr_total), KT(nr_total), ddVV(nr_total), Fstep(nr_total), &
            KTCdT(nr_total), ddFF(nr_total), ddVVdof(nr_u*nr_v*nr_w), TTn1(nr_total), ddFFdof(nr_u*nr_v*nr_w))
    
    c = 1
    do i = 2, sizeF
        ! Get delta time
        dt = time_list(i) - time_list(i-1)

        ! Get values of last simulation 
        TTn0 = solution(:, i-1)

        ! Prediction of new step
        TTn1 = TTn0 + dt*(1-newmark)*VVn0; TTn1(dod) = solution(dod, i)
        TTn1i0 = TTn1; VVn1 = 0.d0
        VVn1(dod) = 1.d0/newmark*(1.0d0/dt*(solution(dod, i) - solution(dod, i-1)) - (1 - newmark)*VVn0(dod))

        ! Get force of new step
        Fstep = FF(:, i)
        
        ! Solver Newton-Raphson
        print*, 'Step: ', i - 1
        do j = 1, nbIterNL

            ! Compute temperature (at each quadrature point) 
            call interpolate_temperature_3d(nr_u_t, nc_u_t, nr_v_t, nc_v_t, nr_w_t, nc_w_t, nnz_u_t, nnz_v_t, nnz_w_t, &
                                        indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                        data_B_u_t, data_B_v_t, data_B_w_t, TTn1, TTinterp)

            ! Interpolate capacity and conductivity at each quadrature point 
            call compute_heat_properties(nbpts, table_cond, table_cap, nc_total, TTinterp, KK, CC)

            ! Compute coefficients to compute tangent matrix
            call compute_heat_coefficients(nc_total, KK, CC, invJJ, detJJ, Kcoefs, Ccoefs)
            
            ! Compute Fint = C dT + K T 
            call mf_wq_get_cu_3d_csr(Ccoefs, nc_total, nr_u_t, nc_u_t, nr_v_t, nc_v_t, nr_w_t, nc_w_t, &
                                nnz_u_t, nnz_v_t, nnz_w_t, indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                data_B_u_t, data_B_v_t, data_B_w_t, data_W_u_t, data_W_v_t, data_W_w_t, VVn1, CdT)

            call mf_wq_get_ku_3d_csr(Kcoefs, nc_total, nr_u_t, nc_u_t, nr_v_t, nc_v_t, nr_w_t, nc_w_t, &
                                nnz_u_t, nnz_v_t, nnz_w_t, indi_u_t, indj_u_t, indi_v_t, indj_v_t, indi_w_t, indj_w_t, &
                                data_B_u_t, data_B_v_t, data_B_w_t, data_W_u_t, data_W_v_t, data_W_w_t, TTn1, KT)
            KTCdT = KT + CdT

            ! Compute residue
            ddFF = Fstep - KTCdT
            call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, ddFF, ddFFdof)
            resNL = maxval(abs(ddFFdof))
            print*, " with Raphson error: ", resNL
            if (resNL.le.thresholdNL) exit

            ! Solve by iterations 
            resPCG(1, c) = dble(i-1); resPCG(2, c) = dble(j)
            call mf_wq_transient_linear_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    newmark*dt, ddFFdof, nbIterPCG, thresholdPCG, methodPCG, ddVVdof, resPCG(3:, c))

            ! Update values
            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, ddVVdof, ddVV)
            VVn1 = VVn1 + ddVV
            TTn1 = TTn1i0 + newmark*dt*VVn1
            c = c + 1
                
        end do
        
        ! Set values
        solution(:, i) = TTn1
        VVn0 = VVn1

    end do

end subroutine mf_wq_transient_nonlinear_3d
