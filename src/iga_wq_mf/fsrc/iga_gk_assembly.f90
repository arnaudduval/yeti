! ==========================
! module :: assembly for IGA-Galerkin 
! author :: Joaquin Cornejo
! 
! This module computes matrices and vectors using sum-factorization algorithms. 
! These algorithms exploits tensor-product structure of shape functions.
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
    integer :: nnz_I_u_t, nnz_I_v_t, nnz_I_w_t
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v), data_W_w(nnz_w)
    integer :: i

    ! Verify if number of non-zero values are well-defined
    nnz_I_u_t = -1; nnz_I_v_t = -1; nnz_I_w_t = -1
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w_t, indi_I_w, indj_I_w)
    if (nnz_I_u_t.ne.nnz_I_u) stop 'Something went wrong computing indices'
    if (nnz_I_v_t.ne.nnz_I_v) stop 'Something went wrong computing indices'
    if (nnz_I_w_t.ne.nnz_I_w) stop 'Something went wrong computing indices'

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u_t), indj_I_v(nnz_I_v_t), indj_I_w(nnz_I_w_t))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w_t, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                            nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                            nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                            size(data_result), indi_result, indj_result)

    do i = 1, nnz_u
        data_W_u(i) = data_B_u(i, 1) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i) = data_B_v(i, 1) * W_v(indj_v(i))
    end do

    do i = 1, nnz_w
        data_W_w(i) = data_B_w(i, 1) * W_w(indj_w(i))
    end do

    ! Get values
    call csr_get_matrix_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u(:, 1), data_B_v(:, 1), data_B_w(:, 1), &
                        data_W_u, data_W_v, data_W_w, nnz_I_u, nnz_I_v, nnz_I_w, & 
                        indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size(data_result), data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

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
    integer :: nnz_I_u_t, nnz_I_v_t, nnz_I_w_t
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 2), data_W_v(nnz_v, 2), data_W_w(nnz_w, 2)
    double precision :: data_temp
    dimension :: data_temp(nnz_I_u*nnz_I_v*nnz_I_w)
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)

    ! Verify if number of non-zero values are well-defined
    nnz_I_u_t = -1; nnz_I_v_t = -1; nnz_I_w_t = -1
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w_t, indi_I_w, indj_I_w)
    if (nnz_I_u_t.ne.nnz_I_u) stop 'Something went wrong computing indices'
    if (nnz_I_v_t.ne.nnz_I_v) stop 'Something went wrong computing indices'
    if (nnz_I_w_t.ne.nnz_I_w) stop 'Something went wrong computing indices'

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u_t), indj_I_v(nnz_I_v_t), indj_I_w(nnz_I_w_t))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w_t, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                                nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                size(data_result), indi_result, indj_result)

    do i = 1, nnz_u
        data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
    end do

    do i = 1, nnz_w
        data_W_w(i, :) = data_B_w(i, :) * W_w(indj_w(i))
    end do

    ! Compute non zero data
    data_result = 0.d0
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            call csr_get_matrix_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                        data_W_u(:, alpha(1)), data_W_v(:, alpha(2)), data_W_w(:, alpha(3)), &
                        nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                        indi_result, size(data_result), data_temp)
            data_result = data_result + data_temp
        end do
    end do
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

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
    double precision :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v), data_W_w(nnz_w)
    integer :: i

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
    integer :: nnz_I_u_t, nnz_I_v_t
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v)
    integer :: i

    ! Verify if number of non-zero values are well-defined
    nnz_I_u_t = -1; nnz_I_v_t = -1
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    if (nnz_I_u_t.ne.nnz_I_u) stop 'Something went wrong computing indices'
    if (nnz_I_v_t.ne.nnz_I_v) stop 'Something went wrong computing indices'

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u_t), indj_I_v(nnz_I_v_t))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                    nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                    size(data_result), indi_result, indj_result)

    do i = 1, nnz_u
        data_W_u(i) = data_B_u(i, 1) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i) = data_B_v(i, 1) * W_v(indj_v(i))
    end do

    ! Get values
    call csr_get_matrix_2d(coefs, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B_u(:, 1), data_B_v(:, 1), data_W_u, data_W_v, &
                        nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size(data_result), data_result)
    deallocate(indj_I_u, indj_I_v)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

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
    integer :: nnz_I_u_t, nnz_I_v_t
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 2), data_W_v(nnz_v, 2)
    double precision :: data_temp
    dimension :: data_temp(nnz_I_u*nnz_I_v)
    integer :: i, j, alpha, beta
    dimension :: alpha(d), beta(d)

    ! Verify if number of non-zero values are well-defined
    nnz_I_u_t = -1; nnz_I_v_t = -1
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)
    if (nnz_I_u_t.ne.nnz_I_u) stop 'Something went wrong computing indices'
    if (nnz_I_v_t.ne.nnz_I_v) stop 'Something went wrong computing indices'

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u_t), indj_I_v(nnz_I_v_t))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u_t, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v_t, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                    nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                    size(data_result), indi_result, indj_result)

    do i = 1, nnz_u
        data_W_u(i, :) = data_B_u(i, :) * W_u(indj_u(i))
    end do

    do i = 1, nnz_v
        data_W_v(i, :) = data_B_v(i, :) * W_v(indj_v(i))
    end do

    ! Compute non zero data
    data_result = 0.d0
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            call csr_get_matrix_2d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B_u(:, beta(1)), data_B_v(:, beta(2)), &
                        data_W_u(:, alpha(1)), data_W_v(:, alpha(2)), &
                        nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, size(data_result), data_temp)
            data_result = data_result + data_temp
        end do
    end do
    deallocate(indj_I_u, indj_I_v)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine iga_get_conductivity_2d

subroutine iga_get_source_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, W_u, W_v, array_out)
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
    double precision :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u), data_W_v(nnz_v)
    integer :: i

    ! Calculate equivalent weight
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
