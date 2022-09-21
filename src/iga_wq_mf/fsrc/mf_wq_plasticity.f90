! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
! 
! This module solves elasto-plasticity problems
! By the moment, it uses a return mapping algorithm based on 
! combined isotropic/kinematic hardening theory. Further information 
! in "Introduction to Nonlinear Finite Element Analysis" by Nam-Ho kim
! Moreover, it uses Newton-Raphson algorithm to solve non-linear equation
! Further information in "Computational Inelasticity" by Simo and Hughes.
! ==========================

subroutine interpolate_strain_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, invJ, u, eps)
    !! Computes strain in 3D (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: d = 3, ddl = d*(d+1)/2
    integer, intent(in) :: nc_total, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: invJ, u
    dimension :: invJ(d, d, nc_total), u(d, nr_u*nr_v*nr_w)

    double precision, intent(out) :: eps
    dimension :: eps(ddl, nc_total)

    ! Local data
    !-----------
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer :: i, j, k, beta(d)
    double precision :: EE, ders, eps_temp
    dimension :: EE(ddl, d, d), ders(d*d, nc_total), eps_temp(ddl)

    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Construct passage matrix
    call create_incidence_matrix(d, ddl, EE)

    ! Compute derivatives of displacement field with respect to u (in parametric space)
    do j = 1, d
        do i = 1, d
            beta = 1; beta(i) = 2
            call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            u(j, :), ders(i+(j-1)*d, :))
        end do
    end do

    ! Compute symetric derivatives (strain) 
    do k = 1, nc_total

        ! Initialize 
        eps_temp = 0.d0

        ! Compute E invJ result
        do i = 1, d
            eps_temp = eps_temp + matmul(EE(:, :, i), matmul(transpose(invJ(:, :, k)), ders((i-1)*d+1:i*d, k)))
        end do

        ! Save data
        eps(:, k) = eps_temp

    end do

end subroutine interpolate_strain_3d

subroutine wq_get_forcevol_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
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

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i

    do i = 1, d
        call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            coefs(i, :), array_out(i, :))
    end do

end subroutine wq_get_forcevol_3d

subroutine wq_get_forcesurf_3d(vforce, JJ, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) ::  nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: vforce, JJ
    dimension :: vforce(d), JJ(d, d-1, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(d, nc_total), v1(d), v2(d), v3(d)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total
        ! Define vectors
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        
        ! Compute the delta surface 
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))

        ! Compute coefficients
        coefs(:, i) = vforce * dsurf
    end do

    ! Compute force
    do i = 1, d
        call sumproduct2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs(i, :), array_out(i, :))
    end do
    
end subroutine wq_get_forcesurf_3d

subroutine wq_get_forceint_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_out)
    !! Computes internal force vector in 3D 
    !! Probably correct (?)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
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

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    double precision :: array_temp
    dimension :: array_temp(nr_u*nr_v*nr_w)
    integer :: i, j, k, alpha(d)
    
    ! Initialize
    array_out = 0.d0

    ! Compute internal force
    do j = 1, d
        do i = 1, d
        
            k = i + (j-1)*d
            alpha = 1; alpha(i) = 4
            call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, alpha(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, alpha(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, alpha(3)), &
                                coefs(k, :), array_temp)
            array_out(j, :) = array_out(j, :) + array_temp
            
        end do
    end do
    
end subroutine wq_get_forceint_3d

subroutine mf_wq_get_su_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    double precision, intent(in) :: array_in
    dimension :: array_in(d, nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    double precision :: array_temp
    dimension :: array_temp(nr_u*nr_v*nr_w)
    integer :: i, j
    
    ! Initialize
    array_out = 0.d0

    do i = 1, d
        do j = 1, d
            
            call mf_wq_get_ku_3d(coefs((i-1)*d+1:i*d, (j-1)*d+1:j*d, :), nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in(j, :), array_temp)

            array_out(i, :) = array_out(i, :) + array_temp    

        end do 
    end do

end subroutine mf_wq_get_su_3d

subroutine mf_wq_get_su_3d_csr(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix. 
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    double precision, intent(in) :: array_in
    dimension :: array_in(d, nr_u*nr_v*nr_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    ! Convert CSR to CSC
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_wq_get_su_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_su_3d_csr

subroutine mf_wq_elasticity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isPrecond, nbIterPCG, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)
    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! IN CSR FORMAT                
    
    implicit none 
    ! Input / output data
    ! -------------------
    logical, intent(in) :: isPrecond 
    double precision, parameter :: threshold = 1.d-7
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    dimension :: U_u(nr_u, nr_u, d), U_v(nr_v, nr_v, d), U_w(nr_w, nr_w, d), Deigen(d, nr_u*nr_v*nr_w)
    integer, intent(in) :: ndu, ndv, ndw, nbIterPCG
    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    double precision, intent(in) :: b
    dimension :: b(d, nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x 
    dimension :: x(d, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    ! Conjugate gradient algorithm
    double precision :: rsold, rsnew, alpha, omega, beta, prod, prod2, normb, RelRes
    double precision, allocatable, dimension(:, :) :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, Ap, As
    integer :: nr_total, iter

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
    allocate(r(d, nr_total), rhat(d, nr_total), p(d, nr_total), s(d, nr_total), &
            ptilde(d, nr_total), Aptilde(d, nr_total), Astilde(d, nr_total), & 
            stilde(d, nr_total), Ap(d, nr_total), As(d, nr_total))
    x = 0.d0; r = b; 
    call clean_dirichlet_3dim(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    rhat = r; p = r
    call block_dot_product(d, nr_total, r, rhat, rsold)
    normb = maxval(abs(r))

    if (.not.isPrecond) then
        ! -------------------------------------------
        ! Conjugate gradient algorithm
        ! -------------------------------------------
        do iter = 1, nbIterPCG
            call mf_wq_get_su_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call clean_dirichlet_3dim(nr_total, Ap, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, Ap, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Ap ! Normally s is alrady Dirichlet updated

            call mf_wq_get_su_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            call clean_dirichlet_3dim(nr_total, As, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, As, s, prod)
            call block_dot_product(d, nr_total, As, As, prod2)
            omega = prod/prod2
            x = x + alpha*p + omega*s ! Normally x is alrady Dirichlet updated
            r = s - omega*As ! Normally r is alrady Dirichlet updated

            RelRes = maxval(abs(r))/normb
            if (RelRes.le.threshold) exit
            call block_dot_product(d, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do
    else
        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        do iter = 1, nbIterPCG
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, p, ptilde)
            call clean_dirichlet_3dim(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call clean_dirichlet_3dim(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            
            call block_dot_product(d, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
            
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, s, stilde)
            call clean_dirichlet_3dim(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call clean_dirichlet_3dim(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(d, nr_total, Astilde, s, prod)
            call block_dot_product(d, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde ! Normally x is alrady Dirichlet updated
            r = s - omega*Astilde ! Normally r is alrady Dirichlet updated
            
            RelRes = maxval(abs(r))/normb
            if (RelRes.le.threshold) exit
            call block_dot_product(d, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do
    end if
    print*, 'Bi CGstab with error: ', RelRes, 'after iterations', iter
end subroutine mf_wq_elasticity_3d

subroutine mf_wq_plasticity_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, properties, &
                        table, ndu, ndv, ndw, dod_u, dod_v, dod_w, invJ, detJ, sizeF, Fext, disp)
    !! Solves elasto-plasticity problems using combined isotropic/kinematic hardening theory
    !! and Newton-Raphson method for non-linear equations
    !! IN CSR FORMAT
    !! IT HAS NOT BEEN TESTED

    use elastoplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: nbIterNL = 30, nbIterPCG = 200
    integer, parameter :: d = 3, dof = d*(d+1)/2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: properties(5)
    integer, intent(in) :: ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    double precision, intent(in) :: invJ, detJ, Fext
    dimension :: invJ(d, d, nc_total), detJ(nc_total), Fext(d, nr_u*nr_v*nr_w, sizeF)
    
    double precision, intent(out) :: disp
    dimension :: disp(d, nr_u*nr_v*nr_w, sizeF)

    ! Local data
    ! ----------
    type(mecamat), pointer :: mat
    double precision :: E, H, beta, nu, sigma_Y

    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:, :), allocatable :: D_u, D_v, D_w
    double precision, dimension(:), allocatable :: I_u, I_v, I_w
    double precision :: s_u, s_v, s_w

    double precision :: alpha_n0, alpha_n1, ep_n1, ep_n0, deps, sigma_n0, sigma_n1, Dalg, coef_fint, coef_S
    dimension ::    alpha_n0(dof, nc_total), alpha_n1(dof, nc_total), ep_n0(nc_total), ep_n1(nc_total), &
                    deps(dof, nc_total), sigma_n0(dof, nc_total), sigma_n1(dof, nc_total), &
                    Dalg(dof, dof, nc_total), coef_fint(d*d, nc_total), coef_S(d*d, d*d, nc_total)
    
    double precision, allocatable, dimension(:, :) :: Fstep, Fint, dF, ddisp, delta_disp 
    double precision :: relerror, prod
    integer :: i, j, k, nr_total

    ! Set total number of rows
    nr_total = nr_u*nr_v*nr_w
    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_w.le.0)) stop 'Indices must be greater than 0'

    ! --------------------
    ! Eigen decomposition
    ! -------------------- 
    ! Initialize 
    allocate(U_u(nr_u, nr_u, d), D_u(nr_u, d), U_v(nr_v, nr_v, d), D_v(nr_v, d), &
            U_w(nr_w, nr_w, d), D_w(nr_w, d), Deigen(d, nr_total))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, d
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), table(1, :, i), D_u(:, i), U_u(:, :, i), Kdiag_u, Mdiag_u)
        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), table(2, :, i), D_v(:, i), U_v(:, :, i), Kdiag_v, Mdiag_v)
        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), table(3, :, i), D_w(:, i), U_w(:, :, i), Kdiag_w, Mdiag_w) 
    end do
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    ! ------------------------
    ! Solver non linear system
    ! ------------------------
    ! Initialize
    allocate(Fstep(d, nr_total), Fint(d, nr_total), dF(d, nr_total), ddisp(d, nr_total), delta_disp(d, nr_total))
    E = properties(1); H = properties(2);  beta = properties(3); nu = properties(4); sigma_Y = properties(5)
    call initialize_mecamat(mat, E, H, beta, nu, sigma_Y)
    disp = 0.d0; ep_n0 = 0.d0; sigma_n0 = 0.d0
    
    do i = 2, sizeF

        ! Initialize displacement
        ddisp = 0.d0

        ! Get force of new step
        Fstep = Fext(:, :, i)        

        ! Solver Newton-Raphson
        print*, 'Step: ', i - 1
        do j = 1, nbIterNL

            ! Compute strain as a function of displacement (at each quadrature point) 
            call interpolate_strain_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                            invJ, ddisp, deps)

            ! Closest point projection in perfect plasticity 
            do k = 1, nc_total
                call cpp_combined_hardening(mat, deps(:, k), alpha_n0(:, k), ep_n0(k), sigma_n0(:, k), &
                                    alpha_n1(:, k), ep_n1(k), sigma_n1(:, k), Dalg(:, :, k))
            end do

            ! Compute coefficients to compute Fint and Stiffness
            call compute_meca_coefficients(nc_total, sigma_n1, Dalg, invJ, detJ, coef_fint, coef_S)

            ! Compute Deigen
            do k = 1, d
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*d+1, (k-1)*d+1, :), s_u)
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*d+2, (k-1)*d+2, :), s_v)
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*d+3, (k-1)*d+3, :), s_w)
                call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, k), D_v(:, k), D_w(:, k), &
                                            s_u, s_v, s_w, Deigen(k, :))
            end do
            
            ! Compute Fint
            call wq_get_forceint_3d(coef_fint, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, Fint)

            dF = Fstep - Fint
            call clean_dirichlet_3dim(nr_total, dF, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(d, nr_total, dF, dF, prod)
            relerror = sqrt(prod)
            print*, "Raphson with error: ", relerror
            if (relerror.le.1e-8) exit
    
            ! Solver Bi-conjugate gradient
            call mf_wq_elasticity_3d(coef_S, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, .true., &
                                nbIterPCG, U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, &
                                dF, delta_disp)
            ddisp = ddisp + delta_disp
        end do
        
        ! Save values
        disp(:, :, i) = disp(:, :, i-1) + ddisp
        ep_n0 = ep_n1
        sigma_n0 = sigma_n1
        alpha_n0 = alpha_n1
    
    end do

end subroutine mf_wq_plasticity_3d

subroutine mf_wq_elasticity_3d_py(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            isPrecond, nbIterPCG, table, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)
    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! This function is adapted to python 
    !! IN CSR FORMAT 

    use elastoplasticity
    implicit none 
    ! Input / output data
    ! ---------------------
    logical, intent(in) :: isPrecond 
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    integer, intent(in) :: nbIterPCG, ndu, ndv, ndw
    integer, intent(in) :: table, dod_u, dod_v, dod_w
    dimension :: table(d, 2, d), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    
    double precision, intent(in) :: b
    dimension :: b(d, nr_u*nr_v*nr_w)
    
    double precision, intent(out) :: x 
    dimension :: x(d, nr_u*nr_v*nr_w)

    ! Local data
    ! -----------
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w
    integer :: nr_total, i

    ! Set total number of rows
    nr_total = nr_u*nr_v*nr_w

    ! --------------------
    ! Eigen decomposition
    ! --------------------
    ! Initialize 
    allocate(U_u(nr_u, nr_u, d), D_u(nr_u), U_v(nr_v, nr_v, d), D_v(nr_v), U_w(nr_w, nr_w, d), D_w(nr_w), Deigen(d, nr_total))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, d
        call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                                data_W_u(:, 4), table(1, :, i), D_u, U_u(:, :, i), Kdiag_u, Mdiag_u)

        call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                                data_W_v(:, 4), table(2, :, i), D_v, U_v(:, :, i), Kdiag_v, Mdiag_v)

        call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                                data_W_w(:, 4), table(3, :, i), D_w, U_w(:, :, i), Kdiag_w, Mdiag_w) 

        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, 1.d0, 1.d0, 1.d0, Deigen(i, :))
    end do
    deallocate(I_u, I_v, I_w, D_u, D_v, D_w)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    call mf_wq_elasticity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isPrecond, nbIterPCG, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, x)

end subroutine mf_wq_elasticity_3d_py