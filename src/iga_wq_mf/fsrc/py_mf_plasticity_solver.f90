! ==============================================
! This module solves elasto-plasticity problems
! Further information in "Introduction to Nonlinear Finite Element Analysis" by Nam-Ho kim
! Moreover, it uses Newton-Raphson algorithm to solve non-linear equation
! Further information in "Computational Inelasticity" by Simo and Hughes.
!
! Remarks :: tensor notation is used (it is NOT Voigt notation)
! Further information in "B free" by J. Planas, I. Romero and J.M. Sancho
! We save some memory storing only the upper triangular of a symmetric matrix
! For example:
! Strain e = [e11, e22, e33, e12, e13, e23]
! Stress s = [s11, s22, s33, s12, s13, s23]
! In this way, we try to avoid some misconceptions when using Voigt notation
! ==============================================

subroutine interpolate_strain_2d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, invJ, u, isvoigt, strain)
    !! Computes strain in 3D (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------  
    integer, parameter :: dimen = 2, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nr_v, nc_u, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: invJ, u
    dimension :: invJ(dimen, dimen, nc_total), u(dimen, nr_u*nr_v)
    logical, intent(in) :: isvoigt

    double precision, intent(out) :: strain
    dimension :: strain(nvoigt, nc_total)

    ! Local data
    !-----------
    integer :: indi_T_u, indi_T_v
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
    integer :: indj_T_u, indj_T_v
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    integer :: i, j, k, beta(dimen)
    double precision :: dersu, dersx
    dimension :: dersu(dimen, dimen, nc_total), dersx(dimen, dimen)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    ! Compute jacobian of displacement field with respect to u (in parametric space)
    do i = 1, dimen
        do k = 1, dimen
            beta = 1; beta(k) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
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

end subroutine interpolate_strain_2d

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

subroutine wq_get_forcevol_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    integer :: i

    do i = 1, dimen
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            coefs(i, :), array_out(i, :))
    end do

end subroutine wq_get_forcevol_2d

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

subroutine wq_get_forcesurf_2d(vforce, JJ, nc_total, nr_u, nc_u, nnz_u, indi_u, indj_u, data_W_u, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nnz_u
    double precision, intent(in) :: vforce, JJ
    dimension :: vforce(dimen, nc_total), JJ(dimen, dimen-1, nc_total)
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_W_u
    dimension :: data_W_u(nnz_u, 4)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, W00
    dimension :: coefs(dimen, nc_total), v1(dimen), W00(nr_u, nc_total)
    integer :: i

    if (nc_total.ne.nc_u) stop 'Not possible finde force surf'

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        dsurf = sqrt(dot_product(v1, v1))
        coefs(:, i) = vforce(:, i) * dsurf
    end do

    ! Compute force
    W00 = 0.d0
    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_total, W00)

    do i = 1, dimen
        array_out(i, :) = matmul(W00, coefs(i, :))
    end do
    
end subroutine wq_get_forcesurf_2d

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

subroutine wq_get_intforce_2d(stress, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v,  &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, &
                            invJ, detJ, array_out)
    !! Computes internal force vector in 3D 
    !! Probably correct (?)
    !! IN CSR FORMAT

    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: stress
    dimension :: stress(nvoigt, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total) 

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    type(mecamat), pointer :: mat
    integer :: nr_total

    allocate(mat)
    mat%dimen = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    nr_total = nr_u*nr_v
    call wq_intforce_2d(mat, stress, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_out)

end subroutine wq_get_intforce_2d

subroutine wq_get_intforce_3d(stress, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            invJ, detJ, array_out)
    !! Computes internal force vector in 3D 
    !! Probably correct (?)
    !! IN CSR FORMAT

    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: stress
    dimension :: stress(nvoigt, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total) 

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    type(mecamat), pointer :: mat
    integer :: nr_total

    allocate(mat)
    mat%dimen = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    nr_total = nr_u*nr_v*nr_w
    call wq_intforce_3d(mat, stress, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_out)

end subroutine wq_get_intforce_3d

subroutine fd_elasticity_2d(nr_total, nr_u, nr_v, U_u, U_v, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    use solverplasticity2
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nr_u, nr_v
    double precision, intent(in) :: U_u, U_v, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), &
                    eigen_diag(dimen, nr_total), array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call setup_preconditionerdiag(solv, dimen, nr_total, eigen_diag)
    solv%isdiagblocks = .true.
    call applyfastdiag(solv, nr_total, nr_u, nr_v, U_u, U_v, array_in, array_out)
    
end subroutine fd_elasticity_2d

subroutine fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    use solverplasticity3
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), U_w(nr_w, nr_w, dimen), &
                    eigen_diag(dimen, nr_total), array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call setup_preconditionerdiag(solv, dimen, nr_total, eigen_diag)
    solv%isdiagblocks = .true.
    call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
    
end subroutine fd_elasticity_3d

subroutine mf_wq_get_su_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            invJ, detJ, CepArgs, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix. 
    !! This function is adapted to python and ONLY for elastric materials
    !! IN CSR FORMAT

    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    double precision :: invJ, detJ, CepArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), CepArgs(nvoigt+3, nc_total)

    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data 
    ! ----------
    type(mecamat), pointer :: mat
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    allocate(mat)
    mat%dimen  = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, CepArgs)
    call mf_wq_stiffness_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

end subroutine mf_wq_get_su_2d

subroutine mf_wq_get_su_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            invJ, detJ, CepArgs, array_in, array_out)
    !! Computes S.u in 3D where S is stiffness matrix. 
    !! This function is adapted to python and ONLY for elastric materials
    !! IN CSR FORMAT

    use matrixfreeplasticity
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

    double precision :: invJ, detJ, CepArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), CepArgs(nvoigt+3, nc_total)

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
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    allocate(mat)
    mat%dimen  = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, CepArgs)
    call mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_wq_get_su_3d

subroutine mf_wq_elasticity_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, b, &
                            ndu, ndv, dod_u, dod_v, table, invJ, detJ, properties, CepArgs, &
                            nbIterPCG, threshold, methodPCG, x, resPCG)

    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! This function is adapted to python 
    !! IN CSR FORMAT 

    use matrixfreeplasticity
    use solverplasticity2
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    integer, intent(in) :: ndu, ndv
    integer, intent(in) :: table, dod_u, dod_v
    dimension :: table(dimen, 2, dimen), dod_u(ndu), dod_v(ndv)

    double precision, intent(in) :: invJ, detJ, properties(3), CepArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), CepArgs(nvoigt+3, nc_total) 
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold 
    logical :: isdiagblocks = .true.

    double precision, intent(in) :: b
    dimension :: b(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(mecamat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision, dimension(:, :, :), allocatable :: U_u, U_v
    double precision, dimension(:, :), allocatable :: Deigen, D_u, D_v

    integer :: i, j, c, ddl
    double precision :: I_u, I_v
    dimension :: I_u(nr_u), I_v(nr_v)

    ! Csr format
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)
    
    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'

    call initialize_mecamat(mat, dimen, properties(1), properties(2), properties(3))
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, CepArgs)
    allocate(solv)

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                data_W_u, data_W_v, ndu, ndv, dod_u, dod_v, nbIterPCG, threshold, b, x, resPCG)

    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C')) then

        if ((isdiagblocks).and.(methodPCG.eq.'JMC')) then
            call compute_mean_diagblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v/))
        else if (.not.(isdiagblocks).and.(methodPCG.eq.'JMC')) then
            call compute_mean_allblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v/))
        end if

        if (isdiagblocks.or.(methodPCG.eq.'C')) then

            ddl = dimen
            allocate(U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl), D_u(nr_u, ddl), D_v(nr_v, ddl))
            call eigendecomp_plasticity_2d(nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                    data_B_u, data_B_v, data_W_u, data_W_v, table, &
                                    ddl, U_u, U_v, D_u, D_v)

            allocate(Deigen(ddl, nr_total))
            I_u = 1.d0; I_v = 1.d0
            do i = 1, dimen
                call find_parametric_diag_2d(nr_u, nr_v, I_u, I_v, D_u(:, i), D_v(:, i), mat%mean(i, i, :), Deigen(i, :))
            end do
        
        else 

            ddl = dimen*dimen
            allocate(U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl), D_u(nr_u, ddl), D_v(nr_v, ddl))
            call eigendecomp_plasticity_2d(nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                    data_B_u, data_B_v, data_W_u, data_W_v, table, &
                                    ddl, U_u, U_v, D_u, D_v)

            allocate(Deigen(ddl, nr_total))
            I_u = 1.d0; I_v = 1.d0; c = 0
            do i = 1, dimen
                do j = 1, dimen
                    c = c + 1
                    call find_parametric_diag_2d(nr_u, nr_v, I_u, I_v, D_u(:, c), D_v(:, c), mat%mean(i, j, :), Deigen(c, :))
                end do
            end do

        end if

        call setup_preconditionerdiag(solv, ddl, nr_total, Deigen)
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, U_u, U_v, ndu, ndv, dod_u, dod_v, &
                        nbIterPCG, threshold, b, x, resPCG)
    else 
        stop 'Unknown method'                    
    end if

end subroutine mf_wq_elasticity_2d

subroutine mf_wq_elasticity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, b, &
                            ndu, ndv, ndw, dod_u, dod_v, dod_w, table, invJ, detJ, properties, CepArgs, &
                            nbIterPCG, threshold, methodPCG, x, resPCG)

    !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
    !! This algorithm solve S x = F, where S is the stiffness matrix
    !! Moreover, it considers Dirichlet boundaries are zero
    !! This function is adapted to python 
    !! IN CSR FORMAT 

    use matrixfreeplasticity
    use solverplasticity3
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

    double precision, intent(in) :: invJ, detJ, properties(3), CepArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), CepArgs(nvoigt+3, nc_total) 
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold 
    logical :: isdiagblocks = .true.

    double precision, intent(in) :: b
    dimension :: b(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(mecamat), pointer :: mat
    type(cgsolver), pointer :: solv
    double precision, dimension(:, :, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:, :), allocatable :: Deigen, D_u, D_v, D_w

    integer :: i, j, c, ddl
    double precision :: I_u, I_v, I_w
    dimension :: I_u(nr_u), I_v(nr_v), I_w(nr_w)

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

    if (any(dod_u.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_v.le.0)) stop 'Indices must be greater than 0'
    if (any(dod_w.le.0)) stop 'Indices must be greater than 0'

    call initialize_mecamat(mat, dimen, properties(1), properties(2), properties(3))
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, CepArgs)
    allocate(solv)

    if (methodPCG.eq.'WP') then 
        ! Set solver
        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, ndu, ndv, ndw, dod_u, dod_v, dod_w, nbIterPCG, threshold, b, x, resPCG)

    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C')) then

        if ((isdiagblocks).and.(methodPCG.eq.'JMC')) then
            call compute_mean_diagblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v, nc_w/))
        else if (.not.(isdiagblocks).and.(methodPCG.eq.'JMC')) then
            call compute_mean_allblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v, nc_w/))
        end if

        if (isdiagblocks.or.(methodPCG.eq.'C')) then

            ddl = dimen
            allocate(U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl), U_w(nr_w, nr_w, ddl), &
                    D_u(nr_u, ddl), D_v(nr_v, ddl), D_w(nr_w, ddl))
            call eigendecomp_plasticity_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                                        data_W_u, data_W_v, data_W_w, table, ddl, U_u, U_v, U_w, D_u, D_v, D_w)

            allocate(Deigen(ddl, nr_total))
            I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
            do i = 1, dimen
                call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, i), D_v(:, i), D_w(:, i), &
                                            mat%mean(i, i, :), Deigen(i, :))
            end do
        
        else 

            ddl = dimen*dimen
            allocate(U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl), U_w(nr_w, nr_w, ddl), &
                    D_u(nr_u, ddl), D_v(nr_v, ddl), D_w(nr_w, ddl))
            call eigendecomp_plasticity_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                                        data_W_u, data_W_v, data_W_w, table, ddl, U_u, U_v, U_w, D_u, D_v, D_w)

            allocate(Deigen(ddl, nr_total))
            I_u = 1.d0; I_v = 1.d0; I_w = 1.d0; c = 0
            do i = 1, dimen
                do j = 1, dimen
                    c = c + 1
                    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u(:, c), D_v(:, c), D_w(:, c), &
                                            mat%mean(i, j, :), Deigen(c, :))
                end do
            end do

        end if

        call setup_preconditionerdiag(solv, ddl, nr_total, Deigen)
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                        U_u, U_v, U_w, ndu, ndv, ndw, dod_u, dod_v, dod_w, nbIterPCG, threshold, b, x, resPCG)
    else 
        stop 'Unknown method'                   
    end if

end subroutine mf_wq_elasticity_3d
