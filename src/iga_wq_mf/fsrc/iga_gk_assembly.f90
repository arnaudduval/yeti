! ==========================
! module :: assembly for IGA-Galerkin 
! author :: Joaquin Cornejo
! 
! This module computes matrices and vectors using sum-factorization algorithms. 
! These algorithms exploits tensor-product structure of shape functions.
! ATTENTION :   this file is not supposed to be used since we want to optimize 
!               the calculation time by not forming the matrix.
! ==========================

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine iga_get_capacity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, &
                                nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)
    !! Computes capacity matrix in 3D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, W_u, data_B_v, W_v, data_B_w, W_w
    dimension ::    data_B_u(nnz_u, 2), W_u(nc_u), &
                    data_B_v(nnz_v, 2), W_v(nc_v), &
                    data_B_w(nnz_w, 2), W_w(nc_w)
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1), indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ----------
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    call iga2wq3d(nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, indj_u, indj_v, indj_w, &
                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, data_W_u, data_W_v, data_W_w)

    call wq_get_capacity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)

end subroutine iga_get_capacity_3d

subroutine iga_get_conductivity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                    data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, &
                                    nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, W_u, data_B_v, W_v, data_B_w, W_w
    dimension ::    data_B_u(nnz_u, 2), W_u(nc_u),&
                    data_B_v(nnz_v, 2), W_v(nc_v), &
                    data_B_w(nnz_w, 2), W_w(nc_w)
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nr_u*nr_v*nr_w+1), &
                    indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ----------
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    call iga2wq3d(nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w, indj_u, indj_v, indj_w, &
                data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, data_W_u, data_W_v, data_W_w)

    call wq_get_conductivity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)

end subroutine iga_get_conductivity_3d

subroutine iga_get_source_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, W_u, W_v, W_w, array_out)
    !! Computes source vector in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: W_u, W_v, W_w
    dimension :: W_u(nc_u), W_v(nc_v), W_w(nc_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    integer :: i
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v), data_W_w(nnz_w)
    
    ! Calculate equivalent weight
    do i = 1, nnz_u
        data_W_u(i) = data_B_u(i, 1) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i) = data_B_v(i, 1) * W_v(indj_v(i))
    end do

    do i = 1, nnz_w
        data_W_w(i) = data_B_w(i, 1) * W_w(indj_w(i))
    end do

    ! Compute vector 
    call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, indi_u, indj_u, data_W_u, &
                        nnz_v, indi_v, indj_v, data_W_v, &
                        nnz_w, indi_w, indj_w, data_W_w, &
                        coefs, array_out)

end subroutine iga_get_source_3d

! ----------------------------------------
! Assembly in 2D
! ----------------------------------------

subroutine iga_get_capacity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, W_u, W_v, &
                                nnz_I_u, nnz_I_v, data_result, indi_result, indj_result)
    !! Computes capacity matrix in 2D 
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, W_u, data_B_v, W_v
    dimension ::    data_B_u(nnz_u, 2), W_u(nc_u), &
                    data_B_v(nnz_v, 2), W_v(nc_v)
    integer, intent(in) :: nnz_I_u, nnz_I_v

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v+1), indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ----------
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    call iga2wq2d(nc_u, nc_v, nnz_u, nnz_v, indj_u, indj_v, &
                data_B_u, data_B_v, W_u, W_v, data_W_u, data_W_v)

    call wq_get_capacity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            nnz_I_u, nnz_I_v, data_result, indi_result, indj_result)

end subroutine iga_get_capacity_2d

subroutine iga_get_conductivity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                    indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, W_u, W_v, &
                                    nnz_I_u, nnz_I_v, data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 2D 
    !! IN CSR FORMAT
                
    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 2 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, W_u, data_B_v, W_v
    dimension ::    data_B_u(nnz_u, 2), W_u(nc_u), &
                    data_B_v(nnz_v, 2), W_v(nc_v)
    integer, intent(in) :: nnz_I_u, nnz_I_v

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v+1), indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ----------
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    call iga2wq2d(nc_u, nc_v, nnz_u, nnz_v, indj_u, indj_v, &
                data_B_u, data_B_v, W_u, W_v, data_W_u, data_W_v)

    call wq_get_conductivity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, &
                                nnz_I_u, nnz_I_v, data_result, indi_result, indj_result)

end subroutine iga_get_conductivity_2d

subroutine iga_get_source_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, W_u, W_v, array_out)
    !! Computes source vector in 2D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2)
    double precision, intent(in) :: W_u, W_v
    dimension :: W_u(nc_u), W_v(nc_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data
    ! ----------
    integer :: i
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v)

    do i = 1, nnz_u
        data_W_u(i) = data_B_u(i, 1) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i) = data_B_v(i, 1) * W_v(indj_v(i))
    end do

    ! Compute vector 
    call sumproduct2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, data_W_u, &
                        nnz_v, indi_v, indj_v, data_W_v, coefs, array_out)

end subroutine iga_get_source_2d
