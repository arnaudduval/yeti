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
                                data_B_u, data_B_v, data_B_w, invJ, u, isvoigt, strain)
    !! Computes strain in 3D (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: invJ, u
    dimension :: invJ(dimen, dimen, nc_total), u(dimen, nr_u*nr_v*nr_w)
    logical, intent(in) :: isvoigt

    double precision, intent(out) :: strain
    dimension :: strain(nvoigt, nc_total)

    ! Local data
    !-----------
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    integer :: i, j, k, beta(dimen)
    double precision :: dersu, dersx
    dimension :: dersu(dimen, dimen, nc_total), dersx(dimen, dimen)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Compute jacobian of displacement field with respect to u (in parametric space)
    do i = 1, dimen
        do k = 1, dimen
            beta = 1; beta(k) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                            u(i, :), dersu(i, k, :))
        end do
    end do

    ! Compute jacobian with respect to x (in physical space)
    do j = 1, nc_total
        dersx = matmul(dersu(:, :, j), invJ(:, :, j))
        dersx = 0.5*(dersx + transpose(dersx))
        call symtensor2array(dimen, nvoigt, dersx, strain(:, j))
    end do

    if (.not.isvoigt) return
    do i = dimen+1, nvoigt
        strain(i, :) = 2*strain(i, :)
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
        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
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
    dimension :: vforce(dimen, nc_total), JJ(dimen, dimen-1, nc_total)
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
        coefs(:, i) = vforce(:, i) * dsurf
    end do

    ! Compute force
    do i = 1, dimen
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs(i, :), array_out(i, :))
    end do
    
end subroutine wq_get_forcesurf_3d

subroutine wq_get_forceint_3d_py(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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
    integer :: i, j, r, zeta(dimen)
    
    array_out = 0.d0
    do i = 1, dimen
        do j = 1, dimen
        
            r = j + (i-1)*dimen
            zeta = 1; zeta(j) = 3
            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                coefs(r, :), array_temp)
            array_out(i, :) = array_out(i, :) + array_temp
            
        end do
    end do

end subroutine wq_get_forceint_3d_py

subroutine fd_elasticity_3d_py(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u, d), U_v(nr_v, nr_v, d), U_w(nr_w, nr_w, d), &
                    eigen_diag(d, nr_total), array_in(d, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_total)

    call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    
end subroutine fd_elasticity_3d_py

subroutine mf_wq_get_su_3d_py(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            invJ, detJ, properties, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix. 
    !! This function is adapted to python and ONLY for elastric materials
    !! IN CSR FORMAT

    use linearelastoplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision :: invJ, detJ, properties(2)
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total)  

    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data 
    ! ----------
    type(mecamat), pointer :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    double precision :: kwargs(nvoigt+3, nc_total)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call initialize_mecamat(mat, dimen, properties(1), properties(1), 1.d0, properties(2), huge(0.0d0))
    call setup_geo(mat, nc_total, invJ, detJ)
    mat%isElastic = .true.

    kwargs = 0.d0; kwargs(1, :) = mat%lambda; kwargs(2, :) = mat%mu
    call setup_kwargs(mat, nc_total, kwargs)
    
    call mf_wq_get_su_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_su_3d_py

subroutine mf_wq_elasticity_3d_py(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, b, &
                            ndu, ndv, ndw, dod_u, dod_v, dod_w, table, invJ, detJ, properties, &
                            nbIterPCG, threshold, methodPCG, x, resPCG)

    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! This function is adapted to python 
    !! IN CSR FORMAT 

    use linearelastoplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    double precision :: invJ, detJ, properties(2)
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total)    
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold 

    double precision, intent(in) :: b
    dimension :: b(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(mecamat), pointer :: mat
    logical :: isPrecond
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    double precision, dimension(:), allocatable :: I_u, I_v, I_w
    double precision :: U_u, U_v, U_w, D_u, D_v, D_w, Deigen, kwargs
    dimension ::    U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), U_w(nr_w, nr_w, dimen), &
                    D_u(nr_u, dimen), D_v(nr_v, dimen), D_w(nr_w, dimen), Deigen(dimen, nr_total), kwargs(nvoigt+3, nc_total)
    double precision :: mean(3)
    integer :: i

    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_w.le.0)) stop 'Indices must be greater than 0'
    isPrecond = .true.
    if (methodPCG.eq.'WP') isPrecond = .false.

    call initialize_mecamat(mat, dimen, properties(1), properties(1), 1.d0, properties(2), huge(0.0d0))
    call setup_geo(mat, nc_total, invJ, detJ)
    mat%isElastic = .true.

    ! Eigen decomposition
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w)) 
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

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
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    kwargs = 0.d0; kwargs(1, :) = mat%lambda; kwargs(2, :) = mat%mu
    call setup_kwargs(mat, nc_total, kwargs)
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    if (methodPCG.eq.'JMC') call compute_mean_plasticity_3d(mat, nc_u, nc_v, nc_w)

    do i = 1, dimen
        mean = 1.d0
        if (methodPCG.eq.'JMC') mean = mat%mean(i, :)
        call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, i), D_v(:, i), D_w(:, i), mean, Deigen(i, :))
    end do
    deallocate(I_u, I_v, I_w)

    call mf_wq_elasticity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, &
                            nbIterPCG, threshold, isPrecond, x, resPCG)

end subroutine mf_wq_elasticity_3d_py
