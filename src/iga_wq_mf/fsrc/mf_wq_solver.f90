! ==========================
! module: Matrix free solver
! author: Joaquin Cornejo
! ==========================

subroutine mf_iga_steady_heat_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, b, nbIter, epsilon, & 
                                method, nnz_cond, cond, JJ, directsol, x, RelRes, RelError)
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
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, b, JJ, directsol
    dimension :: b(nr_u*nr_v*nr_w), JJ(3, 3, nc_total), directsol(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_u*nr_v*nr_w), RelRes(nbIter+1), RelError(nbIter+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, norm2b
    double precision, allocatable, dimension(:) :: r, p, Ap, dummy, z
    integer :: nr_total, iter, i

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: c_u, c_v, c_w

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
    norm2b = norm2(b)
    if (norm2b.lt.epsilon) stop 'Force is almost zero, then it is a trivial solution' 

    if (method.eq.'WP') then 
        if (nbIter.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            p = r
            rsold = dot_product(r, r)
            
            do iter = 1, nbIter
                call mf_iga_get_Ku_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                    data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, p, Ap)
                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                RelRes(iter+1) = norm2(r)/norm2b
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.epsilon) exit
            
                rsnew = dot_product(r, r)
                p = r + rsnew/rsold * p
                rsold = rsnew
            end do
        end if
    else 
        ! Initialize
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))
        Mcoef_u = 1.d0; Kcoef_u = 1.d0; c_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0; c_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0; c_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! Jacobian mean
            call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, Lu, Lv, Lw)
            call conductivity_mean_3d(nc_u, nc_v, nc_w, nnz_cond, cond, lamb_u, lamb_v, lamb_w)
                    
            c_u = lamb_u*Lv*Lw/Lu
            c_v = lamb_v*Lu*Lw/Lv
            c_w = lamb_w*Lu*Lv/Lw

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
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, c_u, c_v, c_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, c_u, c_v, c_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call iga_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, &
                                nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, Dphysical)
        end if
        
        if (nbIter.gt.0) then
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

            do iter = 1, nbIter
                call mf_iga_get_Ku_3D(coefs, nr_total, nc_total, nr_u, nc_u, &
                            nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, W_u, W_v, W_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, p, Ap)

                alpha = rsold/dot_product(p, Ap)
                x = x + alpha * p
                r = r - alpha * Ap

                RelRes(iter+1) = norm2(r)/norm2b
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.epsilon) exit       

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
                            b, nbIter, epsilon, method, nnz_cond, cond, JJ, directsol, x, RelRes, RelError)
    
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
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, b, JJ, directsol
    dimension :: b(nr_u*nr_v*nr_w), JJ(3, 3, nc_total), directsol(nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x, RelRes, RelError
    dimension :: x(nr_u*nr_v*nr_w), RelRes(nbIter+1), RelError(nbIter+1)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, norm2b
    double precision, allocatable, dimension(:) :: r, rhat, p, Ap, s, As, dummy, ptilde, Aptilde, stilde, Astilde
    integer :: nr_total, iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: c_u, c_v, c_w

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
    rsold = dot_product(r, rhat); norm2b = norm2(b)
    if (norm2b.lt.epsilon) stop 'Force is almost zero, then it is a trivial solution' 

    if (method.eq.'WP') then 
        if (nbIter.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            do iter = 1, nbIter
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
    
                RelRes(iter+1) = norm2(r)/norm2b
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.epsilon) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew
            end do
        end if

    else  
        ! Initialize
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0; c_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0; c_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0; c_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! Diagonal decomposition
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! Jacobian mean
            call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, Lu, Lv, Lw)
            call conductivity_mean_3d(nc_u, nc_v, nc_w, nnz_cond, cond, lamb_u, lamb_v, lamb_w)

            c_u = lamb_u*Lv*Lw/Lu
            c_v = lamb_v*Lu*Lw/Lv
            c_w = lamb_w*Lu*Lv/Lw

        end if

        ! --------------------------------------------
        ! Eigen decomposition
        ! -------------------------------------------- 
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
        allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
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
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, c_u, c_v, c_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, c_u, c_v, c_w, Dparametric)
            deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

            ! Find diagonal of real matrix (K in this case)
            allocate(Dphysical(nr_total))
            call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
        end if

        if (nbIter.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            do iter = 1, nbIter

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

                RelRes(iter+1) = norm2(r)/norm2b
                RelError(iter+1) = maxval(abs(directsol - x))/maxval(abs(directsol))
                if (RelRes(iter+1).le.epsilon) exit

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
    !! Solves MM.s = r (MM is the preconditioner) in steady heat 3D (or transient) case with substitution method, where 
    !! MM is an approximation of M = K.
    !! Indices must be in CSR format 

    implicit none
    ! Input / output  data 
    !---------------------
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
    double precision :: LTu, KLTu
    dimension :: LTu(nr_total), KLTu(nr_total)

    ! Compute L'.u                   
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, LTu)

    ! Scaling
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, LTu)
    end if

    ! By now, we test this approximation. It could change later
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, LTu, KLTu)

    ! Scaling
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, KLTu)
    end if

    ! Compute L.u                   
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, KLTu, array_out)

end subroutine fd_tshs_3d

! STEADY HEAT TRANSFER CASE:
! -------------------------------------------------------------

subroutine mf_wq_get_Au_shs_3d(Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_L, indj_L, indi_LT, indj_LT, L, LT, &
                                ndof, array_in, array_out)
    !! Computes Ann.u in steady heat 3D case with substitution method, where 
    !! But A is given as [Ann, And; Adn, Add]. So Ann u =  L A L' u, where C is a zeros and ones matrix.
    !! Indices must be in CSR format

    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
    double precision, intent(in) :: Kcoefs
    dimension :: Kcoefs(d, d, nc_total)

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
    dimension :: indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof)
    double precision :: L, LT
    dimension :: L(ndof), LT(ndof)

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data 
    ! ------------------   
    double precision :: KLTu, LTu
    dimension :: KLTu(nr_total), LTu(nr_total)

    ! Compute L'.u                   
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, LTu)

    ! Compute K.u
    call mf_wq_get_ku_3d(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, LTu, KLTu)

    ! Compute L.u                   
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, KLTu, array_out)

end subroutine mf_wq_get_Au_shs_3d

subroutine mf_wq_solve_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, method, nnz_cond, cond, JJ, x, residue)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition). Then, we use substitution method, so
    !! Mnn xn = Fn - Mnd xd and xd = g 
    !! CSR FORMAT

    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
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
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, f, g, JJ
    dimension :: f(nr_total), g(ndod), JJ(3, 3, nc_total)
    
    double precision, intent(out) :: x, residue
    dimension :: x(nr_total), residue(nbIter+1)

    ! Local data
    ! ------------------
    ! Pre / Conjugate gradient algoritm
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: rsold, rsnew, alpha, omega, beta, norm2b
    double precision :: r, rhat, p, Ap, s, As, ptilde, Aptilde, stilde, Astilde, xn, Kx, fAx
    dimension ::    r(nr_total-ndod), rhat(nr_total-ndod), p(nr_total-ndod), Ap(nr_total-ndod), &
                    As(nr_total-ndod), s(nr_total-ndod), ptilde(nr_total-ndod), Aptilde(nr_total-ndod), &
                    Astilde(nr_total-ndod), stilde(nr_total-ndod), xn(nr_total-ndod), Kx(nr_total), fAx(nr_total)
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: c_u, c_v, c_w

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

    ! Initiate variables
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    ndof = nr_total - ndod
    x = 0.d0; x(dod) = g

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! Compute K . x where x = [0, xd], then K.x = [Knd xd, Kdd xd]
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Kx)
    fAx = f - Kx ! This is the real b in Ax = b equation
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, fAx, r)

    ! Set variables
    xn = 0.d0; rhat = r; p = r
    rsold = dot_product(r, rhat); norm2b = rsold
    residue = 0.d0; residue(1) = 1.d0
    if (norm2b.lt.epsilon) stop 'Fext is almost zero, then it is a trivial solution' 
    
    if (method.eq.'WP') then 
        ! ----------------------------
        ! Conjugate Gradient algorithm
        ! ----------------------------
        do iter = 1, nbIter
            call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, p, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, s, As)
            omega = dot_product(As, s)/dot_product(As, As)
            xn = xn + alpha*p + omega*s
            r = s - omega*As

            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
            x(dod) = g
            
            residue(iter+1) = dot_product(r, r)/norm2b
            if (residue(iter+1).le.epsilon) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    else  
        ! Initialize
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0; c_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0; c_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0; c_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! DIAGONAL DECOMPOSITION
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! MY METHOD
            ! Find dimensions and conductivity
            call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, Lu, Lv, Lw)
            call conductivity_mean_3d(nc_u, nc_v, nc_w, nnz_cond, cond, lamb_u, lamb_v, lamb_w)
            
            c_u = lamb_u*Lv*Lw/Lu
            c_v = lamb_v*Lu*Lw/Lv
            c_w = lamb_w*Lu*Lv/Lw

        end if

        ! --------------------------------------------
        ! EIGEN DECOMPOSITION
        ! -------------------------------------------- 
        allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
        allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
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
        allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
        allocate(Deigen(nr_total))
        I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, c_u, c_v, c_w, Deigen)
        deallocate(I_u, I_v, I_w)

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
            ! Find diagonal of preconditioner
            allocate(Dparametric(nr_total))
            call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, c_u, c_v, c_w, Dparametric)
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
        do iter = 1, nbIter

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, p, ptilde)

            call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, s, stilde)

            call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            xn = xn + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
            x(dod) = g

            residue(iter+1) = dot_product(r, r)/norm2b
            if (residue(iter+1).le.epsilon) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew

        end do

    end if

end subroutine mf_wq_solve_shs_3d

! TRANSIENT HEAT TRANSFER CASE: 
! -------------------------------------------------------------

subroutine mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_L, indj_L, indi_LT, indj_LT, L, LT, &
                                ndof, newmarkdt, array_in, array_out)
    !! Computes Ann.u in steady heat 3D case with substitution method, where 
    !! But A is given as [Ann, And; Adn, Add]. So Ann u =  L A L' u, where L is a zeros and ones matrix.
    !! Indices must be in CSR format

    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
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
    dimension :: indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof)
    double precision :: L, LT, newmarkdt
    dimension :: L(ndof), LT(ndof)

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data 
    ! ------------------   
    double precision :: KCLTu, LTu
    dimension :: KCLTu(nr_total), LTu(nr_total)

    ! Compute v = L'.u                   
    call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, array_in, LTu)

    ! Compute w = (C + alpha dt K).v
    call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            1.d0, newmarkdt, LTu, KCLTu)

    ! Compute L.w                 
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, KCLTu, array_out)

end subroutine mf_wq_get_Au_ths_3d

subroutine mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            newmarkdt, table, ndod, dod, f, g, nbIter, epsilon, x, residue)

    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d = 3
    character (len=10) :: method = 'FDC'
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
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

    double precision, intent(in) :: newmarkdt
    integer, intent(in) :: table, dod
    dimension :: table(d, 2), dod(ndod)
    integer, intent(in) :: nbIter
    double precision, intent(in) :: epsilon, f, g
    dimension :: f(nr_total), g(ndod)
    
    double precision, intent(out) :: x, residue
    dimension :: x(nr_total), residue(nbIter+1)

    ! Local data
    ! ------------------
    ! Pre / Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta, bdotb
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, xn, Ax, fAx
    dimension ::    r(nr_total-ndod), rhat(nr_total-ndod), p(nr_total-ndod), &
                    s(nr_total-ndod), ptilde(nr_total-ndod), Aptilde(nr_total-ndod), &
                    Astilde(nr_total-ndod), stilde(nr_total-ndod), xn(nr_total-ndod), Ax(nr_total), fAx(nr_total)
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

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

    ! Initiate variables
    if (any(dod.le.0)) stop 'Indices must be greater than 0'
    ndof = nr_total - ndod
    x = 0.d0; x(dod) = g

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! Compute A . x where x = [0, xd] and A = (C + alpha dt K)
    call mf_wq_get_kcu_3d(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, 1.d0, newmarkdt, x, Ax)

    fAx = f - Ax ! This is the real b in Ax = b equation
    call spMdotdV(ndof, nr_total, ndof, indi_L, indj_L, L, fAx, r)

    ! Set variables
    xn = 0.d0; rhat = r; p = r
    rsold = dot_product(r, rhat); bdotb = rsold
    residue = 0.d0; residue(1) = 1.d0
    if (bdotb.lt.epsilon) stop 'Fext is almost zero, then it is a trivial solution' 

    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
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
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    allocate(Deigen(nr_total))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, 1.d0, 1.d0, 1.d0, Deigen)

    ! Update eigen diagonal
    Deigen = 1.d0 + newmarkdt*Deigen
    deallocate(I_u, I_v, I_w)

    ! --------------------------------------------
    ! SCALING
    ! --------------------------------------------
    ! Find diagonal of preconditioner
    allocate(Dparametric(nr_total))
    Dparametric = 1.d0

    ! Find diagonal of real matrix (K in this case)
    allocate(Dphysical(nr_total))
    Dphysical = 1.d0

    if (nbIter.gt.0) then
        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        ! The system to solve is Knn xn = bn
        ! Where bn = fn - Knd xd 
        ! So, to solve this system we initialize r = bn - Knn xn, if xn = 0, r = bn
        do iter = 1, nbIter

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, p, ptilde)

            call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, & 
                                nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, newmarkdt, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                            ndof, indi_L, indj_L, indi_LT, indj_LT, L, LT, s, stilde)

            call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, &
                                nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_L, indj_L, indi_LT, indj_LT, L, LT, ndof, newmarkdt, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            xn = xn + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            call spMdotdV(nr_total, ndof, ndof, indi_LT, indj_LT, LT, xn, x)
            x(dod) = g
            
            residue(iter+1) = dot_product(r, r)/bdotb
            if (residue(iter+1).le.epsilon) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew

        end do

    end if

end subroutine mf_wq_solve_ths_linear_3d

subroutine mf_wq_ths_nonlinear_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, nbpts, table_cond, table_cap, &
                        newmark, table_dir, ndod, dod, invJ, detJ, sizeF, time_list, FF, GG, temperature)
    
    use heat_transfer
    implicit none 
    ! Input / output data
    ! ---------------------
    double precision, parameter :: tol = 1.d-8
    ! Geometry
    integer, parameter :: nbIterRaphson = 30, nbIterSolver = 100, d = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF, nbpts
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    ! Physics
    double precision, intent(in) :: table_cond, table_cap, newmark
    dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2)
    integer, intent(in) :: ndod, table_dir, dod
    dimension :: table_dir(d, 2), dod(ndod)
    double precision, intent(in) :: time_list, invJ, detJ, FF, GG
    dimension :: time_list(sizeF), invJ(d, d, nc_total), detJ(nc_total), FF(nr_total, sizeF), GG(ndod, sizeF)
    
    double precision, intent(out) :: temperature
    dimension :: temperature(nr_total, sizeF)

    ! Local data
    ! -----------   
    integer :: ndof, i, j
    double precision :: TTn0, TTn1, TTn1i0, TTinterp, ddTTn0, ddTTn1, CdT, KT, Fit, KTCdT, ddFF, sol
    dimension :: TTn0(nr_total), TTn1(1, nr_total), TTn1i0(nr_total), TTinterp(1, nc_total), &
                ddTTn0(nr_total), ddTTn1(nr_total), CdT(nr_total), KT(nr_total), Fit(nr_total), &
                KTCdT(nr_total), ddFF(nr_total), sol(nr_total)
    double precision, allocatable, dimension(:) :: GGtmp, ddGG, VVn0
    double precision :: resSolver, resRaphson, prod1, prod2, dt, dt2, factor
    dimension :: resSolver(nbIterSolver+1)

    integer, allocatable, dimension(:) :: indi_L, indj_L, indi_LT, indj_LT
    double precision, allocatable, dimension(:) :: L, LT

    double precision :: KK, CC, Kcoefs, CCoefs
    dimension :: KK(nc_total), CC(nc_total), Kcoefs(dimen, dimen, nc_total), Ccoefs(nc_total)

    ! Initialize
    ndof = nr_total - ndod
    allocate(GGtmp(ndod), ddGG(ndod), VVn0(ndof))
    GGtmp = 0.d0; ddGG = 0.d0; VVn0 = 0.d0

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
    ddTTn0 = 0.d0; ddTTn0(dod) = ddGG
    deallocate(ddGG)

    ! Create block L
    allocate(indi_L(ndof+1), indj_L(ndof), indi_LT(nr_total+1), indj_LT(ndof), L(ndof), LT(ndof))
    call create_block_L(ndof, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)     

    ! --------------------------------------------
    ! SOLVE
    ! -------------------------------------------- 
    ! Initialize
    temperature = 0.d0; 
    
    do i = 2, sizeF+1
        ! Get delta time
        dt = time_list(i) - time_list(i-1)

        ! Get values of last simulation 
        TTn0 = temperature(:, i-1)

        ! Prediction of new step
        TTn1(1, :) = TTn0 + dt*(1-newmark)*ddTTn0; TTn1(1, dod) = GG(:, i)
        TTn1i0 = TTn1(1, :); ddTTn1 = 0.d0
        ddTTn1(dod) = 1.d0/newmark*(1.0d0/dt*(GG(:, i) - GG(:, i-1)) - (1 - newmark)*ddTTn0(dod))

        ! Get force of new step
        Fit = FF(:, i)
        prod2 = dot_product(Fit, Fit)

        ! Newton Raphson
        do j = 1, nbIterRaphson
            print*, 'Step: ', i-1, ' Iteration: ', j-1

            ! Compute temperature (at each quadrature point) 
            call interpolate_temperature_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, TTn1(1, :), TTinterp)

            ! Interpolate capacity and conductivity at each quadrature point 
            call compute_heat_properties(nbpts, table_cond, table_cap, nc_total, TTinterp(1, :), KK, CC)

            ! Compute coefficients to compute tangent matrix
            call compute_heat_coefficients(nc_total, KK, CC, invJ, detJ, Kcoefs, Ccoefs)
            
            ! Compute Fint = C dT + K T 
            call mf_wq_get_cu_3d_csr(Ccoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, ddTTn1, CdT)

            call mf_wq_get_ku_3d_csr(Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, TTn1(1, :), KT)

            KTCdT = KT + CdT

            ! Compute residue
            ddFF = Fit - KTCdT
            prod1 = dot_product(ddFF, ddFF)
            resRaphson = sqrt(prod1/prod2)
            print*, "Raphson with error: ", resRaphson
            if (isnan(resRaphson)) stop
            
            ! Verify
            if (resRaphson.le.1e-6) then 
                exit
            else

                ! Solve by iterations 
                call mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                newmark*dt, table_dir, ndod, dod, ddFF, GGtmp, nbIterSolver, tol, sol, resSolver)

                ! Update values
                ddTTn1 = ddTTn1 + sol
                TTn1(1, :) = TTn1i0 + newmark*dt*sol
                
            end if                
        end do
        
        ! Set values
        temperature(:, i) = TTn1(1, :)
        ddTTn0 = ddTTn1
                
    end do

end subroutine mf_wq_ths_nonlinear_3d