! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
! ==========================

subroutine wq_get_grad_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes gradient (?) !!!
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, nc_total)
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(nr_u*nr_v*nr_w)

    ! Local data
    ! ---------------
    integer :: i, alpha
    dimension :: alpha(d)
    double precision :: result_temp
    dimension :: result_temp(nr_u*nr_v*nr_w)

    ! Initialize
    result = 0.d0

    do i = 1, d
        alpha = 1; alpha(i) = 3
        ! Find vector temp
        call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, alpha(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, alpha(2)), &
                            nnz_w, indi_w, indj_w, data_W_w(:, alpha(3)), &
                            coefs(i, :), result_temp)
        result = result + result_temp
    end do

end subroutine wq_get_grad_3d

subroutine wq_get_forceVol_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes force vector in 3D case
    !! coefs defined in python

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, nc_total)
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(d, nr_u*nr_v*nr_w)

    ! Local data
    ! -------------
    double precision :: result_temp
    dimension :: result_temp(nr_u*nr_v*nr_w)
    integer :: i

    do i = 1, d
        call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                    coefs(i, :), result_temp)
        result(i, :) = result_temp   
    end do

end subroutine wq_get_forceVol_3d

subroutine wq_get_forceSurf_3d(force, JJ, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, result)
    !! Computes force vector in 3D case
    !! Coefs defined in python

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) ::  nc_total, nr_u, nc_u, nr_v, nc_v
    double precision, intent(in) :: force, JJ
    dimension :: force(d), JJ(d, d-1, nc_total)
    integer, intent(in) :: nnz_u, nnz_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: result
    dimension :: result(d, nr_u*nr_v)

    ! Local data
    ! -------------
    double precision :: result_temp, coefs, norm, v1, v2, v3
    dimension :: result_temp(nr_u*nr_v), coefs(d, nc_total), v1(d), v2(d), v3(d)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total
        ! Define vectors
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        
        ! Compute norm
        call crossproduct(v1, v2, v3)
        norm = sqrt(dot_product(v3, v3))

        ! Compute coefs
        coefs(:, i) = force * norm

    end do

    do i = 1, d
        call tensor2d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                                data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                coefs(i, :), result_temp)
        result(i, :) = result_temp   
    end do
    
end subroutine wq_get_forceSurf_3d

subroutine wq_get_forceInt_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes force vector in 3D case
    !! By now, we assume that is correct

    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, nc_total)
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(d, nr_u*nr_v*nr_w)

    ! Local data
    ! -------------
    double precision :: result_temp
    dimension :: result_temp(nr_u*nr_v*nr_w)
    integer :: i
    
    ! Compute vector
    do i = 1, d
        call wq_get_grad_3d(coefs((i-1)*d+1:i*d, :), nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, result_temp)

        result(i, :) = result_temp   
    end do

end subroutine wq_get_forceInt_3d

subroutine mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_input, array_output)
    !! Computes S.u in 3D case (Stiffness)
    !! Indices must be in CSR format

    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, d*d, nc_total)

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

    double precision, intent(in) :: array_input
    dimension :: array_input(d, nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(d, nr_total)

    ! Local data 
    ! ------------------
    double precision :: array_temp
    dimension :: array_temp(nr_total)
    integer :: i, j
    
    ! Initiliaze
    array_output = 0.d0

    do i = 1, d
        do j = 1, d
            call mf_wq_get_ku_3d(coefs((i-1)*d+1:i*d, (j-1)*d+1:j*d, :), nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_input(j, :), array_temp)

            array_output(i, :) = array_output(i, :) + array_temp    

        end do 
    end do

end subroutine mf_wq_get_su_3d

subroutine mf_wq_get_su_3d_csr(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, array_input, array_output)
    !! Computes S.u in 3D case (Stiffness)
    !! Indices must be in CSR format

    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, d*d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_input
    dimension :: array_input(d, nr_total)

    double precision, intent(out) :: array_output
    dimension :: array_output(d, nr_total)

    ! Local data 
    ! ------------------
    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_input, array_output)

end subroutine mf_wq_get_su_3d_csr

subroutine wq_mf_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            U_u, U_v, U_w, Deigen, nbIter, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    logical :: isPrecond = .True.
    double precision, parameter :: epsilon = 1e-10
    integer, parameter :: d = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, d*d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u, d), U_v(nr_v, nr_v, d), U_w(nr_w, nr_w, d), Deigen(d, nr_total)
    integer, intent(in) :: ndu, ndv, ndw, nbIter
    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    double precision, intent(in) :: b
    dimension :: b(d, nr_total)
    
    double precision, intent(out) :: x 
    dimension :: x(d, nr_total)

    ! Local data
    ! ------------------
    ! Pre / Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta, RelRes, prod, prod2
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, Ap, As
    dimension ::    r(d, nr_total), rhat(d, nr_total), p(d, nr_total), s(d, nr_total), &
                    ptilde(d, nr_total), Aptilde(d, nr_total), Astilde(d, nr_total), stilde(d, nr_total), &
                    Ap(d, nr_total), As(d, nr_total)
    integer :: iter

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Initiate variables
    x = 0.d0

    if (.not.isPrecond) then
        ! -------------------------------------------
        ! Conjugate gradient algorithm
        ! -------------------------------------------
        ! In this algorithm we assume that A = Snn and b = Fn
        r = b; 
        call update_dirichlet_3d(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call dot_prod_plasticity(d, nr_total, r, rhat, rsold)

        do iter = 1, nbIter
            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call update_dirichlet_3d(nr_total, Ap, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call dot_prod_plasticity(d, nr_total, Ap, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Ap ! Normally s is alrady Dirichlet updated

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            call update_dirichlet_3d(nr_total, As, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call dot_prod_plasticity(d, nr_total, As, s, prod)
            call dot_prod_plasticity(d, nr_total, As, As, prod2)
            omega = prod/prod2
            x = x + alpha*p + omega*s ! Normally x is alrady Dirichlet updated
            r = s - omega*As ! Normally r is alrady Dirichlet updated

            RelRes = maxval(abs(r))/maxval(abs(b))
            ! print*, RelRes
            if (RelRes.le.epsilon) exit
            call dot_prod_plasticity(d, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do
    else
        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        ! In this algorithm we assume that A = Snn and b = Fn
        r = b; 
        call update_dirichlet_3d(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call dot_prod_plasticity(d, nr_total, r, rhat, rsold)

        do iter = 1, nbIter
            call fast_diag_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, p, ptilde)
            call update_dirichlet_3d(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call update_dirichlet_3d(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            
            call dot_prod_plasticity(d, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
            
            call fast_diag_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, s, stilde)
            call update_dirichlet_3d(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call update_dirichlet_3d(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call dot_prod_plasticity(d, nr_total, Astilde, s, prod)
            call dot_prod_plasticity(d, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde ! Normally x is alrady Dirichlet updated
            r = s - omega*Astilde ! Normally r is alrady Dirichlet updated
            
            RelRes = maxval(abs(r))/maxval(abs(b))
            if (RelRes.le.epsilon) exit
            call dot_prod_plasticity(d, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do
    end if

end subroutine wq_mf_elasticity_3d

subroutine wq_mf_elasticity_3d_csr(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, &
                            ndu, ndv, ndw, dod_u, dod_v, dod_w, JJ, b, x)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, d*d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    integer, intent(in) :: ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    double precision, intent(in) :: JJ, b
    dimension :: JJ(3, 3, nc_total), b(d, nr_total)
    
    double precision, intent(out) :: x
    dimension :: x(d, nr_total)

    ! Local data
    ! ----------
    integer :: i
    character(len = 10) :: Method = 'FDC'
    double precision, dimension(:), allocatable :: Mcoef, Kcoef
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    double precision :: dum1, dum2, dum3, dum
    dimension :: dum(3, 3, 1)
    double precision :: c_u, c_v, c_w, Lu, Lv, Lw

    ! Jacobian Mean
    dum = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3, 1/))
    call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, &
                            1, dum, Lu, Lv, Lw, dum1, dum2, dum3)
                                
    c_u = Lv*Lw/Lu
    c_v = Lw*Lu/Lv
    c_w = Lu*Lv/Lw

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    ! Initialize 
    allocate(U_u(nr_u, nr_u, d), D_u(nr_u), U_v(nr_v, nr_v, d), D_v(nr_v), U_w(nr_w, nr_w, d), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Deigen(d, nr_total))

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, d

        call eigen_decomposition(nr_u, nc_u, Mcoef, Kcoef, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), method, table(1, :, i), D_u, U_u(:, :, i), Kdiag_u, Mdiag_u)

        call eigen_decomposition(nr_v, nc_v, Mcoef, Kcoef, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), method, table(2, :, i), D_v, U_v(:, :, i), Kdiag_v, Mdiag_v)

        call eigen_decomposition(nr_w, nc_w, Mcoef, Kcoef, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), method, table(3, :, i), D_w, U_w(:, :, i), Kdiag_w, Mdiag_w) 

        call find_parametric_diag_3d(nr_u, nr_v, nr_w, c_u, c_v, c_w, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen(i, :))

    end do
    deallocate(I_u, I_v, I_w, D_u, D_v, D_w, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    call wq_mf_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        U_u, U_v, U_w, Deigen, 100, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)

end subroutine wq_mf_elasticity_3d_csr

subroutine solver_plasticity(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        table, ndu, ndv, ndw, dod_u, dod_v, dod_w, JJ, E, nu, sigma_Y, sizeF, Fext, disp)

    use elastoplasticity
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    ! Geometry
    integer, parameter :: d = 3, dof = d*(d+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    ! Physics
    double precision, intent(in) :: E, nu, sigma_Y
    integer, intent(in) :: ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    double precision, intent(in) :: JJ, Fext
    dimension :: JJ(d, d, nc_total), Fext(d, nr_total, sizeF+1)
    
    double precision, intent(out) :: disp
    dimension :: disp(d, nr_total, sizeF+1)

    ! Local data
    ! -----------
    double precision :: invJ, detJ, invJtemp, detJtemp
    dimension :: invJ(d, d, nc_total), detJ(nc_total), invJtemp(d, d)

    character(len = 10) :: Method = 'FDC'
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    type(material), pointer :: mat
    double precision :: ep_n1, ep_n0, e_n1, sigma_n1, dSdE
    dimension ::    ep_n1(dof, nc_total), ep_n0(dof, nc_total), e_n1(dof, nc_total), &
                    sigma_n1(dof, nc_total), dSdE(dof, dof, nc_total)
    double precision :: disp_t, Fext_t, coef_fint, coef_S, Fint, delta_F, delta_disp
    dimension ::    disp_t(d, nr_total), Fext_t(d, nr_total), coef_fint(d*d, nc_total), &
                    coef_S(d*d, d*d, nc_total), Fint(d, nr_total), delta_F(d, nr_total), delta_disp(d, nr_total)
    double precision :: c_u, c_v, c_w, error, prod1, prod2
    integer :: i, j, k, nbIterRaphson, nbIterSolver, dorobin(2)

    ! --------------------------------------------
    ! GEOMETRY
    ! -------------------------------------------- 
    do i = 1, nc_total
        call MatrixInv(invJtemp, JJ(:, :, i), detJtemp, d)
        invJ(:, :, i) = invJtemp
        detJ(i) = detJtemp
    end do

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    ! Initialize 
    dorobin = (/0, 0/)
    c_u = 1.d0; c_v = 1.d0; c_w = 1.d0
    allocate(U_u(nr_u, nr_u, d), D_u(nr_u), U_v(nr_v, nr_v, d), D_v(nr_v), U_w(nr_w, nr_w, d), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Deigen(d, nr_total))

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, d
        dorobin = table(1, :, i)
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), method, dorobin, D_u, U_u(:, :, i), Kdiag_u, Mdiag_u)

        dorobin = table(2, :, i)
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), method, dorobin, D_v, U_v(:, :, i), Kdiag_v, Mdiag_v)

        dorobin = table(3, :, i)
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), method, dorobin, D_w, U_w(:, :, i), Kdiag_w, Mdiag_w) 

        call find_parametric_diag_3d(nr_u, nr_v, nr_w, c_u, c_v, c_w, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen(i, :))
    end do
    deallocate(I_u, I_v, I_w, D_u, D_v, D_w)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    ! --------------------------------------------
    ! SOLVE
    ! -------------------------------------------- 
    ! Initialize
    call initialize_mat(mat, E, nu)
    disp = 0.d0; ep_n1 = 0.d0; ep_n0 = 0.d0
    nbIterRaphson = 30; nbIterSolver = 100
    
    do i = 2, sizeF+1

        ! Initialize 
        disp_t = disp(:, :, i-1)
        Fext_t = Fext(:, :, i)        
        call dot_prod_plasticity(d, nr_total, Fext_t, Fext_t, prod2)

        ! Newton Raphson
        do j = 1, nbIterRaphson
            print*, 'Step: ', i-1, ' Iteration: ', j-1

            ! Compute strain as a function of displacement (at each quadrature point)
            call interpolate_strain(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, invJ, disp_t, e_n1)

            ! Closest point projection in perfect plasticity 
            do k = 1, nc_total
                call cpp_perfplasticity(mat%Ctensor, mat%Stensor, sigma_Y, e_n1(:, k), ep_n0(:, k), &
                                        ep_n1(:, k), sigma_n1(:, k), dSdE(:, :, k))
            end do

            ! Compute coefficients to compute Fint and Stiffness
            call compute_coefficients(nc_total, sigma_n1, dSdE, invJ, detJ, coef_fint, coef_S)

            ! Compute Fint
            call wq_get_forceInt_3d(coef_fint, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, Fint)

            ! Update F
            delta_F = Fext_t - Fint
            call update_dirichlet_3d(nr_total, delta_F, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call dot_prod_plasticity(d, nr_total, delta_F, delta_F, prod1)
            error = sqrt(prod1/prod2)
            print*, "with error: ", error

            ! Verify
            if ((error.le.1e-6).or.(j.eq.nbIterRaphson)) then 
                disp(:, :, i) = disp_t
                ep_n0 = ep_n1
                exit
            else
                ! Solve by iterations 
                call wq_mf_elasticity_3d(coef_S, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                    U_u, U_v, U_w, Deigen, nbIterSolver, ndu, ndv, ndw, dod_u, dod_v, dod_w, &
                                    delta_F, delta_disp)

                ! Update displacement
                disp_t = disp_t + delta_disp
            end if
            
        end do
        print*, "------------"
    end do

end subroutine solver_plasticity
