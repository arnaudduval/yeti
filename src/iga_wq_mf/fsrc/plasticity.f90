! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
! ==========================

subroutine interpolate_strain_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, invJ, u, eps)
    !! Computes strain in 3D case (from parametric space to physical space)
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input/ output
    ! --------------------  
    integer, parameter :: d = 3, ddl = d*(d+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: invJ, u
    dimension :: invJ(d, d, nc_total), u(d, nr_total)

    double precision, intent(out) :: eps
    dimension :: eps(ddl, nc_total)

    ! Local data
    !-----------------
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer :: i, j, k, beta
    dimension :: beta(d)
    double precision :: MM, invJext, result, temp
    dimension :: MM(ddl, d*d), invJext(d*d, d*d), result(d*d, nc_total), temp(d*d)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Construct MM
    call create_der2sym(d, ddl, MM)

    ! Compute derivatives of the displacement in physical space
    do j = 1, d
        do i = 1, d
            k = i + (j-1)*d
            beta = 1; beta(i) = 2
            call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            u(j, :), result(k, :))
        end do
    end do

    ! Compute true strain 
    do i = 1, nc_total

        ! Compute inverse of jacobian matrix extended
        invJext = 0.d0
        do j = 1, d
            invJext((j-1)*d+1:j*d, (j-1)*d+1:j*d) = invJ(:, :, i)
        end do
        
        ! Compute invJext dot result
        temp = matmul(invJext, result(:, i)) 

        ! Evaluate MM dot temp
        eps(:, i) = matmul(MM, temp)

    end do

end subroutine interpolate_strain_3d

subroutine wq_get_forcevol_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes volumetric force vector in 3D case

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(d, nr_u*nr_v*nr_w)

    ! Local data
    ! -------------
    integer :: i

    do i = 1, d
        call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                    coefs(i, :), result(i, :))
    end do

end subroutine wq_get_forcevol_3d

subroutine wq_get_forcesurf_3d(vforce, JJ, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, result)
    !! Computes boundary force vector in 3D case

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) ::  nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: vforce, JJ
    dimension :: vforce(d), JJ(d, d-1, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: result
    dimension :: result(d, nr_u*nr_v)

    ! Local data
    ! -------------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(d, nc_total), v1(d), v2(d), v3(d)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total
        ! Define vectors
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        
        ! Compute the surface 
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))

        ! Compute coefs
        coefs(:, i) = vforce * dsurf
    end do

    ! Compute force
    do i = 1, d
        call tensor2d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                                data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                coefs(i, :), result(i, :))
    end do
    
end subroutine wq_get_forcesurf_3d

subroutine wq_get_forceint_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes internal force vector in 3D case
    !! --------- IT IS NOT GOOD

    use tensor_methods
    implicit none 
    ! Input / output 
    ! --------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(d*d, nc_total)
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
    integer :: i, j, k, alpha(d)
    
    ! Initialize
    result = 0.d0

    ! Compute vector
    do j = 1, d
        do i = 1, d
        
            k = i + (j-1)*d
            alpha = 1; alpha(i) = 3
            call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, alpha(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, alpha(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, alpha(3)), &
                                coefs(k, :), result_temp)
            result(j, :) = result(j, :) + result_temp
            
        end do
    end do

end subroutine wq_get_forceint_3d

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

! ------------------------------

subroutine mf_wq_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isPrecond, nbIter, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)
    
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    logical, intent(in) :: isPrecond 
    double precision, parameter :: epsilon = 1.d-7
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
    double precision :: rsold, rsnew, alpha, omega, beta, RelRes, prod, prod2, norm2b
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
        call clean_dirichlet_3d(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call block_dot_product(d, nr_total, r, rhat, rsold)
        norm2b = norm2(r)

        do iter = 1, nbIter
            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call clean_dirichlet_3d(nr_total, Ap, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, Ap, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Ap ! Normally s is alrady Dirichlet updated

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            call clean_dirichlet_3d(nr_total, As, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, As, s, prod)
            call block_dot_product(d, nr_total, As, As, prod2)
            omega = prod/prod2
            x = x + alpha*p + omega*s ! Normally x is alrady Dirichlet updated
            r = s - omega*As ! Normally r is alrady Dirichlet updated

            RelRes = sqrt(norm2(r)/norm2b)
            ! print*, RelRes
            if (RelRes.le.epsilon) exit
            call block_dot_product(d, nr_total, r, rhat, rsnew)
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
        call clean_dirichlet_3d(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call block_dot_product(d, nr_total, r, rhat, rsold)
        norm2b = norm2(r)

        do iter = 1, nbIter
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, p, ptilde)
            call clean_dirichlet_3d(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call clean_dirichlet_3d(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            
            call block_dot_product(d, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
            
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, s, stilde)
            call clean_dirichlet_3d(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3D(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call clean_dirichlet_3d(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, Astilde, s, prod)
            call block_dot_product(d, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde ! Normally x is alrady Dirichlet updated
            r = s - omega*Astilde ! Normally r is alrady Dirichlet updated
            
            RelRes = sqrt(norm2(r)/norm2b)
            ! print*, RelRes
            if (RelRes.le.epsilon) exit
            call block_dot_product(d, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do
    end if
    print*, 'Bi CGstab with error: ', RelRes, 'after iterations', iter
end subroutine mf_wq_elasticity_3d

subroutine mf_wq_plasticity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, properties, &
                        table, ndu, ndv, ndw, dod_u, dod_v, dod_w, invJ, detJ, sizeF, Fext, disp)

    use elastoplasticity
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    ! Geometry
    integer, parameter :: nbIterRaphson = 30, nbIterSolver = 300
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
    double precision, intent(in) :: properties(5)
    integer, intent(in) :: ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    double precision, intent(in) :: invJ, detJ, Fext
    dimension :: invJ(d, d, nc_total), detJ(nc_total), Fext(d, nr_total, sizeF+1)
    
    double precision, intent(out) :: disp
    dimension :: disp(d, nr_total, sizeF+1)

    ! Local data
    ! -----------
    type(material), pointer :: mat
    double precision :: E, H, beta, nu, sigma_Y

    character(len = 10) :: Method = 'FDC'
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    double precision :: alpha_n0, alpha_n1, ep_n1, ep_n0, deps, sigma_n0, sigma_n1, Dalg
    dimension ::    alpha_n0(nc_total), alpha_n1(nc_total), ep_n0(nc_total), ep_n1(nc_total), deps(dof, nc_total), &
                    sigma_n0(dof, nc_total), sigma_n1(dof, nc_total), Dalg(dof, dof, nc_total)
    
    double precision :: coef_fint, coef_S, Fext_t, Fint, dF, ddisp
    dimension ::    coef_fint(d*d, nc_total), coef_S(d*d, d*d, nc_total), &
                    Fext_t(d, nr_total), Fint(d, nr_total), dF(d, nr_total), ddisp(d, nr_total)
    
    double precision :: c_u, c_v, c_w, relerror, prod1, prod2
    integer :: i, j, k

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    ! Initialize 
    c_u = 1.d0; c_v = 1.d0; c_w = 1.d0
    allocate(U_u(nr_u, nr_u, d), D_u(nr_u), U_v(nr_v, nr_v, d), D_v(nr_v), U_w(nr_w, nr_w, d), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Deigen(d, nr_total))

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, d
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), method, table(1, :, i), D_u, U_u(:, :, i), Kdiag_u, Mdiag_u)

        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), method, table(2, :, i), D_v, U_v(:, :, i), Kdiag_v, Mdiag_v)

        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), method, table(3, :, i), D_w, U_w(:, :, i), Kdiag_w, Mdiag_w) 

        call find_parametric_diag_3d(nr_u, nr_v, nr_w, c_u, c_v, c_w, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen(i, :))
    end do
    deallocate(I_u, I_v, I_w, D_u, D_v, D_w)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    ! --------------------------------------------
    ! SOLVE
    ! -------------------------------------------- 
    ! Initialize
    E = properties(1); H = properties(2);  beta = properties(3)
    nu = properties(4); sigma_Y = properties(5)
    call initialize_mat(mat, E, H, beta, nu, sigma_Y)
    disp = 0.d0; ep_n0 = 0.d0; sigma_n0 = 0.d0
    
    do i = 2, sizeF+1

        ! Initialize 
        ddisp = 0.d0
        Fext_t = Fext(:, :, i)        
        call block_dot_product(d, nr_total, Fext_t, Fext_t, prod2)

        ! Newton Raphson
        do j = 1, 2
            print*, 'Step: ', i-1, ' Iteration: ', j-1

            ! Compute strain as a function of displacement (at each quadrature point) 
            call interpolate_strain_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, invJ, ddisp, deps)

            ! Closest point projection in perfect plasticity 
            do k = 1, nc_total
                call cpp_combined_hardening(mat, deps(:, k), alpha_n0(k), ep_n0(k), sigma_n0(:, k), &
                                    alpha_n1(k), ep_n1(k), sigma_n1(:, k), Dalg(:, :, k))
            end do

            ! Compute coefficients to compute Fint and Stiffness
            call compute_coefficients(nc_total, sigma_n1, Dalg, invJ, detJ, coef_fint, coef_S)
            
            ! Compute Fint
            call wq_get_forceint_3d(coef_fint, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, Fint)

            dF = Fext_t - Fint
            call clean_dirichlet_3d(nr_total, dF, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(d, nr_total, dF, dF, prod1)
            relerror = sqrt(prod1/prod2)
            print*, "Raphson with error: ", relerror
            if (isnan(relerror)) stop
            
            ! Verify
            if (relerror.le.sigma_Y*1e-6) then 
                exit
            else

                ! Solve by iterations 
                call mf_wq_elasticity_3d(coef_S, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, .true., &
                nbIterSolver, U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, &
                dF, ddisp)

            end if            
        end do
        
        ! Set values
        disp(:, :, i) = disp(:, :, i-1) + ddisp
        ep_n0 = ep_n1
        sigma_n0 = sigma_n1
                
    end do

end subroutine mf_wq_plasticity_3d

subroutine mf_wq_elasticity_3d_py(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            isPrecond, nbIter, table, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)

    use elastoplasticity
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    logical, intent(in) :: isPrecond 
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

    integer, intent(in) :: nbIter, ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    
    double precision, intent(in) :: b
    dimension :: b(d, nr_total)
    
    double precision, intent(out) :: x 
    dimension :: x(d, nr_total)

    ! Local data
    ! -----------
    character(len = 10) :: Method = 'FDC'
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    integer :: i

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
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), method, table(1, :, i), D_u, U_u(:, :, i), Kdiag_u, Mdiag_u)

        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), method, table(2, :, i), D_v, U_v(:, :, i), Kdiag_v, Mdiag_v)

        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), method, table(3, :, i), D_w, U_w(:, :, i), Kdiag_w, Mdiag_w) 

        call find_parametric_diag_3d(nr_u, nr_v, nr_w, 1.d0, 1.d0, 1.d0, &
                                I_u, I_v, I_w, D_u, D_v, D_w, Deigen(i, :))
    end do
    deallocate(I_u, I_v, I_w, D_u, D_v, D_w)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)


    ! Call elasticity
    call mf_wq_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isPrecond, nbIter, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)

end subroutine mf_wq_elasticity_3d_py