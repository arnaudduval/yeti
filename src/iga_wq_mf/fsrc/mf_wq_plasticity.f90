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
    integer, parameter :: dimen = 3, ddl = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: invJ, u
    dimension :: invJ(dimen, dimen, nc_total), u(dimen, nr_u*nr_v*nr_w)

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

    integer :: i, j, k, beta(dimen)
    double precision :: EE, ders, eps_temp
    dimension :: EE(ddl, dimen, dimen), ders(dimen*dimen, nc_total), eps_temp(ddl)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    call create_incidence_matrix(dimen, ddl, EE)

    ! Compute derivatives of displacement field with respect to u (in parametric space)
    do j = 1, dimen
        do i = 1, dimen
            beta = 1; beta(i) = 2
            call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            u(j, :), ders(i+(j-1)*dimen, :))
        end do
    end do

    ! Compute symetric derivatives (strain) 
    do k = 1, nc_total

        eps_temp = 0.d0
        do i = 1, dimen
            eps_temp = eps_temp + matmul(EE(:, :, i), matmul(transpose(invJ(:, :, k)), ders((i-1)*dimen+1:i*dimen, k)))
        end do
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
    integer, parameter :: dimen = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i

    do i = 1, dimen
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
    integer, parameter :: dimen = 3
    integer, intent(in) ::  nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: vforce, JJ
    dimension :: vforce(dimen), JJ(dimen, dimen-1, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, v2, v3
    dimension :: coefs(dimen, nc_total), v1(dimen), v2(dimen), v3(dimen)
    integer :: i

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        v2 = JJ(:, 2, i)
        call crossproduct(v1, v2, v3)
        dsurf = sqrt(dot_product(v3, v3))
        coefs(:, i) = vforce * dsurf
    end do

    ! Compute force
    do i = 1, dimen
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
    integer, parameter :: dimen = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen*dimen, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    double precision :: array_temp
    dimension :: array_temp(nr_u*nr_v*nr_w)
    integer :: i, j, k, alpha(dimen)
    
    array_out = 0.d0
    do j = 1, dimen
        do i = 1, dimen
        
            k = i + (j-1)*dimen
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

subroutine mf_wq_get_su_3d_py(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix. 
    !! This function is adapted to python
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen*dimen, dimen*dimen, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data 
    ! ----------
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_su_3d_py

subroutine mf_wq_elasticity_3d_py(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, b, &
                            ndu, ndv, ndw, dod_u, dod_v, dod_w, table, nbIterPCG, threshold, isPrecond, x, resPCG)

    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! This function is adapted to python 
    !! IN CSR FORMAT 

    use elastoplasticity
    implicit none 
    ! Input / output data
    ! --------------------- 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen*dimen, dimen*dimen, nc_total)
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
    dimension :: table(dimen, 2, dimen), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold
    logical, intent(in) :: isPrecond 

    double precision, intent(in) :: b
    dimension :: b(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: D_u, D_v, D_w, Deigen
    double precision, dimension(:), allocatable :: I_u, I_v, I_w
    double precision :: s_u, s_v, s_w
    integer :: i

    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_w.le.0)) stop 'Indices must be greater than 0'

    ! Eigen decomposition
    allocate(U_u(nr_u, nr_u, dimen), D_u(nr_u, dimen), U_v(nr_v, nr_v, dimen), D_v(nr_v, dimen), &
            U_w(nr_w, nr_w, dimen), D_w(nr_w, dimen), Deigen(dimen, nr_total))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, dimen
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

    do i = 1, dimen
        call compute_mean_3d(nc_u, nc_v, nc_w, coefs((i-1)*dimen+1, (i-1)*dimen+1, :), s_u)
        call compute_mean_3d(nc_u, nc_v, nc_w, coefs((i-1)*dimen+2, (i-1)*dimen+2, :), s_v)
        call compute_mean_3d(nc_u, nc_v, nc_w, coefs((i-1)*dimen+3, (i-1)*dimen+3, :), s_w)
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, i), D_v(:, i), D_w(:, i), &
                                    s_u, s_v, s_w, Deigen(i, :))
    end do

    call mf_wq_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, &
                            nbIterPCG, threshold, isPrecond, x, resPCG)

end subroutine mf_wq_elasticity_3d_py

subroutine mf_wq_plasticity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, sizeF, Fext, &
                        ndu, ndv, ndw, dod_u, dod_v, dod_w, table, invJ, detJ, properties, disp, sigma_vm)
    !! Solves elasto-plasticity problems using combined isotropic/kinematic hardening theory
    !! and Newton-Raphson method for non-linear equations
    !! IN CSR FORMAT

    use elastoplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: nbIterNL = 30, nbIterPCG = 200
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF
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
    dimension :: table(dimen, 2, dimen), dod_u(ndu), dod_v(ndv), dod_w(ndw)
    double precision, intent(in) :: invJ, detJ, Fext
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Fext(dimen, nr_total, sizeF)
    
    double precision, intent(out) :: disp, sigma_vm
    dimension :: disp(dimen, nr_total, sizeF), sigma_vm(nc_total, sizeF)

    ! Local data
    ! ----------
    type(mecamat), pointer :: mat
    double precision :: E, H, beta, nu, sigma_Y

    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: D_u, D_v, D_w, Deigen
    double precision, dimension(:), allocatable :: I_u, I_v, I_w
    double precision :: s_u, s_v, s_w

    double precision :: alpha_n0, alpha_n1, ep_n1, ep_n0, deps, sigma_n0, sigma_n1, Dalg, coef_fint, coef_S
    dimension ::    alpha_n0(ddl, nc_total), alpha_n1(ddl, nc_total), ep_n0(nc_total), ep_n1(nc_total), &
                    deps(ddl, nc_total), sigma_n0(ddl, nc_total), sigma_n1(ddl, nc_total), &
                    Dalg(ddl, ddl, nc_total), coef_fint(dimen*dimen, nc_total), coef_S(dimen*dimen, dimen*dimen, nc_total)
    
    double precision, allocatable, dimension(:, :) :: Fstep, Fint, dF, ddisp, delta_disp 
    double precision :: resNL, resPCG, prod
    dimension :: resPCG(nbIterPCG+1)
    integer :: i, j, k

    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_w.le.0)) stop 'Indices must be greater than 0'

    ! --------------------
    ! Eigen decomposition
    ! -------------------- 
    allocate(U_u(nr_u, nr_u, dimen), D_u(nr_u, dimen), U_v(nr_v, nr_v, dimen), D_v(nr_v, dimen), &
            U_w(nr_w, nr_w, dimen), D_w(nr_w, dimen), Deigen(dimen, nr_total))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0

    do i = 1, dimen
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
    allocate(Fstep(dimen, nr_total), Fint(dimen, nr_total), dF(dimen, nr_total), ddisp(dimen, nr_total), delta_disp(dimen, nr_total))
    E = properties(1); H = properties(2);  beta = properties(3); nu = properties(4); sigma_Y = properties(5)
    call initialize_mecamat(mat, E, H, beta, nu, sigma_Y)
    disp = 0.d0; ep_n0 = 0.d0; sigma_n0 = 0.d0; sigma_vm = 0.d0
    
    do i = 2, sizeF

        ddisp = 0.d0
        Fstep = Fext(:, :, i)        

        print*, 'Step: ', i - 1
        do j = 1, nbIterNL ! Newton-Raphson solver

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
            do k = 1, dimen
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*dimen+1, (k-1)*dimen+1, :), s_u)
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*dimen+2, (k-1)*dimen+2, :), s_v)
                call compute_mean_3d(nc_u, nc_v, nc_w, coef_S((k-1)*dimen+3, (k-1)*dimen+3, :), s_w)
                call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, k), D_v(:, k), D_w(:, k), &
                                            s_u, s_v, s_w, Deigen(k, :))
            end do
            
            ! Compute Fint
            call wq_get_forceint_3d(coef_fint, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, Fint)

            dF = Fstep - Fint
            call clean_dirichlet_3d(nr_total, dF, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(dimen, nr_total, dF, dF, prod)
            resNL = sqrt(prod)
            print*, "Raphson with error: ", resNL
            if (resNL.le.1d-8) exit
    
            ! Solve
            call mf_wq_elasticity_3d(coef_S, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, data_W_w, &
                                U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, dF, &
                                nbIterPCG, 1d-8, .true., delta_disp, resPCG)
            ddisp = ddisp + delta_disp
        end do
        
        ! Save values
        disp(:, :, i) = disp(:, :, i-1) + ddisp
        ep_n0 = ep_n1
        sigma_n0 = sigma_n1
        alpha_n0 = alpha_n1

        do k = 1, nc_total
            call compute_stress_vonmises(dimen, ddl, sigma_n1(:, k), sigma_vm(k, i))
        end do

    end do

end subroutine mf_wq_plasticity_3d