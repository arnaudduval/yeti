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
                                data_B_u, data_B_v, invJ, u, strain)
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

end subroutine interpolate_strain_2d

subroutine interpolate_strain_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, invJ, u, strain)
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

end subroutine interpolate_strain_3d

subroutine get_forcevol_2d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, detJ, prop, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(dimen, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp(nc_total)

    do i = 1, dimen
        tmp = prop(i, :)*detJ
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            tmp, array_out(i, :))
    end do

end subroutine get_forcevol_2d

subroutine get_forcevol_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, detJ, prop, array_out)
    !! Computes volumetric force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    double precision, intent(in) :: detJ, prop
    dimension :: detJ(nc_total), prop(dimen, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp(nc_total)

    do i = 1, dimen
        tmp = prop(i, :)*detJ
        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            tmp, array_out(i, :))
    end do

end subroutine get_forcevol_3d

subroutine get_forcesurf_2d(nc_total, nr_u, nc_u, nnz_u, indi_u, indj_u, data_W_u, JJ, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nnz_u
    integer, intent(in) :: indi_u, indj_u
    dimension :: indi_u(nr_u+1), indj_u(nnz_u)
    double precision, intent(in) :: data_W_u
    dimension :: data_W_u(nnz_u, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(dimen, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u)

    ! Local data
    ! ----------
    double precision :: coefs, dsurf, v1, W00
    dimension :: coefs(dimen, nc_total), v1(dimen), W00(nr_u, nc_total)
    integer :: i

    if (nc_total.ne.nc_u) stop 'Not possible find force surf'

    ! Compute coefficients
    do i = 1, nc_total
        v1 = JJ(:, 1, i)
        dsurf = sqrt(dot_product(v1, v1))
        coefs(:, i) = prop(:, i)*dsurf
    end do

    ! Compute force
    W00 = 0.d0
    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W00)

    do i = 1, dimen
        array_out(i, :) = matmul(W00, coefs(i, :))
    end do
    
end subroutine get_forcesurf_2d

subroutine get_forcesurf_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, JJ, prop, array_out)
    !! Computes boundary force vector in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) ::  nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: JJ, prop
    dimension :: JJ(dimen, dimen-1, nc_total), prop(dimen, nc_total)

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
        coefs(:, i) = prop(:, i) * dsurf
    end do

    ! Compute force
    do i = 1, dimen
        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, &
                        data_W_u(:, 1), nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs(i, :), array_out(i, :))
    end do
    
end subroutine get_forcesurf_3d

subroutine get_intforce_2d(nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v,  &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, &
                            invJ, detJ, stress, array_out)
    !! Computes internal force vector in 3D 
    !! Probably correct (?)
    !! IN CSR FORMAT

    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
    double precision, intent(in) :: invJ, detJ, stress
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), stress(nvoigt, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    type(mecamat) :: mat
    integer :: nr_total

    mat%dimen = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    nr_total = nr_u*nr_v
    call intforce_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, stress, array_out)

end subroutine get_intforce_2d

subroutine get_intforce_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                            invJ, detJ, stress, array_out)
    !! Computes internal force vector in 3D 
    !! Probably correct (?)
    !! IN CSR FORMAT

    use matrixfreeplasticity
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
    double precision, intent(in) :: invJ, detJ, stress
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), stress(nvoigt, nc_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    type(mecamat) :: mat
    integer :: nr_total

    mat%dimen = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    nr_total = nr_u*nr_v*nr_w
    call intforce_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, stress, array_out)

end subroutine get_intforce_3d

subroutine mf_get_su_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            invJ, detJ, MechArgs, array_in, array_out)
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

    double precision :: invJ, detJ, MechArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), MechArgs(nvoigt+3, nc_total)

    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data 
    ! ----------
    type(mecamat) :: mat
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    mat%dimen  = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, MechArgs)
    call mf_stiffness_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

end subroutine mf_get_su_2d

subroutine mf_get_su_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            invJ, detJ, MechArgs, array_in, array_out)
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

    double precision :: invJ, detJ, MechArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), MechArgs(nvoigt+3, nc_total)

    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data 
    ! ----------
    type(mecamat) :: mat
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    mat%dimen  = dimen
    mat%nvoigt = nvoigt
    call setup_geometry(mat, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_jacobiennormal(mat, MechArgs)
    call mf_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

end subroutine mf_get_su_3d

subroutine solver_elasticity_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            ndu, ndv, dod_u, dod_v, table, invJ, detJ, properties, mechArgs, &
                            Fext, nbIterPCG, threshold, methodPCG, x, resPCG)

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
    integer, intent(in) :: dod_u, dod_v
    dimension :: dod_u(ndu), dod_v(ndv)
    logical, intent(in) :: table
    dimension :: table(dimen, 2, dimen) 

    double precision, intent(in) :: invJ, detJ, properties(3), mechArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), mechArgs(nvoigt+3, nc_total) 
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold 

    double precision, intent(in) :: Fext
    dimension :: Fext(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(mecamat) :: mat
    type(cgsolver) :: solv

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
    call setup_jacobiennormal(mat, mechArgs)

    if (methodPCG.eq.'WP') then 

        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                data_W_u, data_W_v, ndu, ndv, dod_u, dod_v, nbIterPCG, threshold, Fext, x, resPCG)

    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C')) then

        if (methodPCG.eq.'JMC') then
            call compute_mean_diagblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v/))
        end if

        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, table, mat%mean)

        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndu, ndv, dod_u, dod_v, &
                        nbIterPCG, threshold, Fext, x, resPCG)
    else 
        stop 'Unknown method'                    
    end if

end subroutine solver_elasticity_2d

subroutine solver_elasticity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            ndu, ndv, ndw, dod_u, dod_v, dod_w, table, invJ, detJ, properties, mechArgs, &
                            Fext, nbIterPCG, threshold, methodPCG, x, resPCG)

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
    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)
    logical, intent(in) :: table
    dimension :: table(dimen, 2, dimen) 

    double precision, intent(in) :: invJ, detJ, properties(3), mechArgs
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), mechArgs(nvoigt+3, nc_total) 
    character(len=10), intent(in) :: methodPCG
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold 

    double precision, intent(in) :: Fext
    dimension :: Fext(dimen, nr_total)
    
    double precision, intent(out) :: x, resPCG
    dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)

    ! Local data
    ! -----------
    type(mecamat) :: mat
    type(cgsolver) :: solv

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
    call setup_jacobiennormal(mat, mechArgs)

    if (methodPCG.eq.'WP') then 

        call BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, ndu, ndv, ndw, dod_u, dod_v, dod_w, nbIterPCG, threshold, Fext, x, resPCG)

    else if ((methodPCG.eq.'JMC').or.(methodPCG.eq.'C')) then

        if (methodPCG.eq.'JMC') then
            call compute_mean_diagblocks(mat, mat%dimen, mat%nvoigt, (/nc_u, nc_v, nc_w/))
        end if

        call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w, table, mat%mean)
        
        call PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                        ndu, ndv, ndw, dod_u, dod_v, dod_w, nbIterPCG, threshold, Fext, x, resPCG)
    else 
        stop 'Unknown method'                   
    end if

end subroutine solver_elasticity_3d
