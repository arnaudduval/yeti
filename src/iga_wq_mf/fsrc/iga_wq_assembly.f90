! ==========================
! module :: assembly for IGA-WQ (using sum-factorization)
! author :: Joaquin Cornejo
! ==========================

! ----------------------------------------
! Assembly in 3D
! ----------------------------------------

subroutine wq_get_capacity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)
    !! Computes a matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1), indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ---------------
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: nnz_result, dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                                nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data 
    nnz_result = nnz_I_u*nnz_I_v*nnz_I_w
    call csr_get_matrix_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u(:, 1), data_B_v(:, 1), data_B_w(:, 1), &
                            data_W_u(:, 1), data_W_v(:, 1), data_W_w(:, 1), &
                            nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, nnz_result, data_result)
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_3d

subroutine wq_get_conductivity_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w 

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v*nnz_I_w)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1), indj_result(nnz_I_u*nnz_I_v*nnz_I_w)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(nnz_I_u*nnz_I_v*nnz_I_w)
    integer :: indi_I_u, indi_I_v, indi_I_w
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v, indj_I_w
    integer :: nnz_result, dummy1, dummy2, dummy3, i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)
    call get_I_csr(nr_w, nc_w, nnz_w, indi_w, indj_w, nnz_I_w, indi_I_w, indj_I_w)

    call get_indexes_kron3_product(nr_w, nr_w, nnz_I_w, indi_I_w, indj_I_w, &
                                nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v*nnz_I_w
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
                        indi_result, nnz_result, data_result_temp)
            data_result = data_result + data_result_temp
        end do
    end do
    deallocate(indj_I_u, indj_I_v, indj_I_w)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_conductivity_3d

subroutine wq_get_source_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, result)
    !! Computes source vector in 3D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_W_u, data_W_v, data_W_w
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

    double precision, intent(out) :: result
    dimension :: result(nr_u*nr_v*nr_w)

    ! Find vector
    call tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                        coefs, result)

end subroutine wq_get_source_3d

subroutine wq_get_stiffness_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            nnz_I_u, nnz_I_v, nnz_I_w, data_result, indi_result, indj_result)
    !! Computes stiffness matrix in 3D case

    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision :: coefs
    dimension :: coefs(d*d, d*d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v, nnz_I_w 

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v*nnz_I_w*d*d)
    integer, intent(out) :: indi_result, indj_result
    dimension ::    indi_result(nnz_I_u*nnz_I_v*nnz_I_w*d*d), &
                    indj_result(nnz_I_u*nnz_I_v*nnz_I_w*d*d)

    ! Local data 
    !-------------
    integer :: size_data_result, nr_total
    double precision, allocatable, dimension(:) :: data_result_temp
    integer, allocatable, dimension(:) :: indi_result_temp_csr, indj_result_temp, indi_result_temp_coo
    integer :: i, j, k, l, nnz, count
    integer :: init, fin

    size_data_result = nnz_I_u*nnz_I_v*nnz_I_w
    nr_total = nr_u*nr_v*nr_w
    allocate(data_result_temp(size_data_result), indi_result_temp_csr(nr_total+1), &
            indj_result_temp(size_data_result), indi_result_temp_coo(size_data_result))
    
    do i = 1, d
        do j = 1, d
            call wq_get_conductivity_3D(coefs((i-1)*d+1:i*d, (j-1)*d+1:j*d, :), nc_total, &
                                    nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                                    data_W_u, data_W_v, data_W_w, nnz_I_u, nnz_I_v, nnz_I_w, &
                                    data_result_temp, indi_result_temp_csr, indj_result_temp)

            if ((i.eq.1).and.(j.eq.1)) then
                count = 1
                indi_result_temp_coo = 0
                do k = 1, nr_total
                    nnz = indi_result_temp_csr(k+1) - indi_result_temp_csr(k)
                    do l = 1, nnz
                        indi_result_temp_coo(count) = k-1
                        count = count + 1
                    end do
                end do
            end if

            init = (j + (i - 1)*d - 1)*size_data_result + 1
            fin = (j + (i - 1)*d)*size_data_result

            data_result(init : fin) = data_result_temp    
            indi_result(init : fin) = indi_result_temp_coo + (i - 1)*nr_total
            indj_result(init : fin) = indj_result_temp + (j - 1)*nr_total

        end do 
    end do

    deallocate(data_result_temp, indi_result_temp_csr, indi_result_temp_coo, indj_result_temp)

end subroutine wq_get_stiffness_3d

! ----------------------------------------
! Assembly in 2D
! ----------------------------------------

subroutine wq_get_capacity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            nnz_I_u, nnz_I_v, & 
                            data_result, indi_result, indj_result)
    !! Computes a matrix in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output 
    ! ------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v

    double precision, intent(out) :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v+1), indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ---------------
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: nnz_result, dummy1, dummy2, dummy3

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v
    call csr_get_matrix_2d(coefs, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u(:, 1), data_B_v(:, 1), data_W_u(:, 1), data_W_v(:, 1), &
                            nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                            indi_result, nnz_result, data_result)
    deallocate(indj_I_u, indj_I_v)

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_capacity_2d

subroutine wq_get_conductivity_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, &
                                nnz_I_u, nnz_I_v, & 
                                data_result, indi_result, indj_result)
    !! Computes conductivity matrix in 2D case
    !! IN CSR FORMAT

    use omp_lib
    use tensor_methods
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 2 
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    integer, intent(in) :: nnz_I_u, nnz_I_v 

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_I_u*nnz_I_v)
    integer, intent(out) :: indi_result, indj_result
    dimension :: indi_result(nr_u*nr_v+1), indj_result(nnz_I_u*nnz_I_v)

    ! Local data
    ! ---------------
    double precision :: data_result_temp
    dimension :: data_result_temp(nnz_I_u*nnz_I_v)
    integer :: indi_I_u, indi_I_v
    dimension :: indi_I_u(nr_u+1), indi_I_v(nr_v+1)
    integer, allocatable, dimension(:) :: indj_I_u, indj_I_v
    integer :: nnz_result, dummy1, dummy2, dummy3, i, j, alpha, beta, zeta
    dimension :: alpha(d), beta(d), zeta(d)

    ! Get indices of I in each dimension
    allocate(indj_I_u(nnz_I_u), indj_I_v(nnz_I_v))
    call get_I_csr(nr_u, nc_u, nnz_u, indi_u, indj_u, nnz_I_u, indi_I_u, indj_I_u)
    call get_I_csr(nr_v, nc_v, nnz_v, indi_v, indj_v, nnz_I_v, indi_I_v, indj_I_v)

    call get_indexes_kron2_product(nr_v, nr_v, nnz_I_v, indi_I_v, indj_I_v, &
                                    nr_u, nr_u, nnz_I_u, indi_I_u, indj_I_u, &
                                    dummy1, dummy2, dummy3, indi_result, indj_result)

    ! Compute non zero data
    nnz_result = nnz_I_u*nnz_I_v
    data_result = 0.d0
    do j = 1, d
        do i = 1, d
            alpha = 1; alpha(i) = 2
            beta = 1; beta(j) = 2
            zeta = beta + (alpha - 1)*2
            call csr_get_matrix_2d(coefs(i, j, :), nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, &
                        data_B_u(:, beta(1)), data_B_v(:, beta(2)), data_W_u(:, zeta(1)), data_W_v(:, zeta(2)), &
                        nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                        indi_result, nnz_result, data_result_temp)
            data_result = data_result + data_result_temp
        end do
    end do
    deallocate(indj_I_u, indj_I_v)   

    ! Fortran to python 
    indi_result = indi_result - 1
    indj_result = indj_result - 1

end subroutine wq_get_conductivity_2d

subroutine wq_get_source_2d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, result)
    !! Computes source vector in 2D case
    !! IN CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_total)
    integer, intent(in) :: nnz_u, nnz_v
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_W_u, data_W_v
    dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

    double precision, intent(out) :: result
    dimension :: result(nr_u*nr_v)

    ! Find vector
    call tensor2d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                        nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                        coefs, result)

end subroutine wq_get_source_2d