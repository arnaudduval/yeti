! ====================================================
! This modules aims to assemble matrices using optimized algorithms 
! exploiting the tensor-product structure of shape functions.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! ====================================================

subroutine test(n, x)

    implicit none 
    integer, intent(in) :: n
    integer, dimension(n), intent(inout) :: x
    integer :: i
    integer :: xt
    do i = 1, n 
        xt = x(i)
        if (xt.gt.0) then
            x(i) = xt*2
        end if
    end do

end subroutine

subroutine wq_get_capacity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            nnz_I_list, nnz, data_result, indi_result, indj_result)
    !! Computes capacity matrix in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, nnz
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    integer, intent(inout) :: nnz_I_list
    dimension :: nnz_I_list(3)

    double precision, intent(out) :: data_result(nnz)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1), indj_result(nnz)

    ! Local data
    ! ----------
    integer :: nnz_I_u, nnz_I_v, nnz_I_w
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w

    data_result = 0.d0; indi_result = 0; indj_result = 0
    if (product(nnz_I_list).ne.nnz) then
        nnz_I_list = -1
        call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_list(1), indi_I_u, indj_I_u)
        call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_list(2), indi_I_v, indj_I_v)
        call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_list(3), indi_I_w, indj_I_w)
        return
    else
        ! Get indices of I in each dimension
        nnz_I_u = nnz_I_list(1); nnz_I_v = nnz_I_list(2); nnz_I_w = nnz_I_list(3) 
        allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
        call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
        call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
        call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

        call get_indices_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                                    nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                    nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                    nnz, indi_result, indj_result)

        ! Compute nonzero data 
        call csr_get_matrix_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u(:, 1), data_B_v(:, 1), data_B_w(:, 1), &
                                data_W_u(:, 1), data_W_v(:, 1), data_W_w(:, 1), &
                                nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, &
                                indj_I_u, indj_I_v, indj_I_w, indi_result, nnz, data_result)
        deallocate(indj_I_u, indj_I_v, indj_I_w)

    end if

end subroutine wq_get_capacity_3d

subroutine wq_get_conductivity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                nnz_I_list, nnz, data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 3D
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, nnz
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
    integer, intent(inout) :: nnz_I_list
    dimension :: nnz_I_list(3)

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1), indj_result(nnz)

    ! Local data
    ! ----------
    integer :: nnz_I_u, nnz_I_v, nnz_I_w
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    double precision, allocatable, dimension(:) :: tmp
    integer :: i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)

    data_result = 0.d0; indi_result = 0; indj_result = 0
    if (product(nnz_I_list).ne.nnz) then
        nnz_I_list = -1
        call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_list(1), indi_I_u, indj_I_u)
        call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_list(2), indi_I_v, indj_I_v)
        call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_list(3), indi_I_w, indj_I_w)
        return
    else
        ! Get indices of I in each dimension
        nnz_I_u = nnz_I_list(1); nnz_I_v = nnz_I_list(2); nnz_I_w = nnz_I_list(3) 
        allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
        call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
        call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
        call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

        call get_indices_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                                    nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                    nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                    nnz, indi_result, indj_result)

        ! Compute nonzero data
        allocate(tmp(product(nnz_I_list)))
        data_result = 0.d0
        do j = 1, d
            do i = 1, d
                alpha = 1; alpha(i) = 2
                beta = 1; beta(j) = 2
                zeta = beta + (alpha - 1)*2
                call csr_get_matrix_3d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_B_w(:, beta(3)), &
                            data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), data_W_w(:, zeta(3)), &
                            nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, nnz, tmp)
                data_result = data_result + tmp
            end do
        end do
        deallocate(indj_I_u, indj_I_v, indj_I_w)

    end if

end subroutine wq_get_conductivity_3d