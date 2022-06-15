! ==========================
! module :: Algebra 
! author :: Joaquin Cornejo
! modules :: ---------
! ==========================
module constants_iga_wq_mf
    
    implicit none
    ! Definition of some constants used in all this project

    integer, parameter :: r = 2 ! For WQ analysis
    double precision, parameter :: span_tol = 1.d-8 ! For spans precision
    double precision, parameter :: tol = 1.d-15 ! For other kind of precision

end module constants_iga_wq_mf

! -------------------
! Vector and matrices
! -------------------
subroutine linspace(x0, xf, n, array) 
    !! Evaluates N equidistant points given the first and last points 

    implicit none
    ! Input / output data
    ! -------------------
    double precision, intent(in) :: x0, xf 
    integer, intent(in) :: n 

    double precision, intent(out) :: array 
    dimension :: array(n)

    ! Local data
    ! -------------
    integer :: i
    double precision :: h 

    ! Define spacing 
    h = (xf - x0)/dble(n - 1)

    ! Assign first and last values
    array(1) = x0 
    array(n) = xf
    
    ! Assign values
    do i = 2, n - 1
        array(i) = x0 + dble(i-1)*h
    end do

end subroutine linspace

subroutine product_AWB(mode, nrA, ncA, A, nrB, ncB, B, W, nrR, ncR, R)
    !! Matrix multiplication 
    !! 1: type:  R = A.diag(W).BT
    !! 2: type:  R = AT.diag(W).B
    !! Matrix A = (nb_rows_A, nb_columns_A)
    !! Array W = (*) it depends 
    !! Matrix B = (nb_rows_B, nb_columns_B)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: mode, nrA, ncA, nrB, ncB, nrR, ncR
    double precision, intent(in) :: A, B, W
    dimension :: A(nrA, ncA), B(nrB, ncB), W(*)

    double precision, intent(out) :: R
    dimension :: R(nrR, ncR)

    ! Local data
    ! -------------
    integer :: i, j, k
    double precision :: s

    ! Initiliaze 
    R = 0.d0

    if (mode.eq.1) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB

        do j = 1, nrB
            do i = 1, nrA
                s = 0.d0
                do k = 1, ncA
                    s = s + A(i, k)*W(k)*B(j, k)
                end do
                R(i, j) = s
            end do
        end do

    else if (mode.eq.2) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB
        ! Only diagonal

        do i = 1, nrA
            s = 0.d0
            do k = 1, ncA
                s = s + A(i, k)*W(k)*B(i, k)
            end do
            R(i, i) = s
        end do

    else if (mode.eq.3) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB

        do j = 1, ncB
            do i = 1, ncA
                s = 0.d0
                do k = 1, nrA
                    s = s + A(k, i)*W(k)*B(k, j)
                end do
                R(i, j) = s
            end do
        end do

    else if (mode.eq.4) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB
        ! Only diagonal

        do i = 1, ncA
            s = 0.d0
            do k = 1, nrA
                s = s + A(k, i)*W(k)*B(k, i)
            end do
            R(i, i) = s
        end do

    end if

end subroutine product_AWB

subroutine solve_system(nr, nc, A, b, x)
    !! Solves equation system A.x = b
    !! Matrix A = (nb_rows, nb_columns)
    !! Vector b = (nb_rows)
    !! Vector x = (nb_cols)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc
    double precision, intent(in) :: A, b
    dimension :: A(nr, nc), b(nr)

    double precision, intent(out) :: x
    dimension :: x(nc)

    ! Local data
    ! ------------------
    double precision :: Atemp
    dimension :: Atemp(nr, nc)
    double precision, allocatable, dimension(:) :: xtemp

    ! Lapack
    integer :: ipiv
    dimension :: ipiv(nr)
    double precision, allocatable, dimension(:) :: work
    integer :: i, info, lwork

    ! Initialize 
    Atemp = A
    x = 0.d0
    
    if (nr.le.nc) then ! Case 1

        allocate(xtemp(nc))
        xtemp = 0.d0
        
        do i = 1, nr
            xtemp(i) = b(i)
        end do

        if (nr.eq.nc) then

            ! Determmined system: 
            call dgesv(nr, 1, Atemp, nr, ipiv, xtemp, nr, info)

        elseif (nr.lt.nc) then

            ! Under-determined system: 
            lwork = 2*nr
            allocate(work(lwork))
            call dgels('N', nr, nc, 1, Atemp, nr, xtemp, nc, work, lwork, info)
                
        end if

        do i = 1, nc
            x(i) = xtemp(i)
        end do

    else if (nr.gt.nc) then ! Case 2

        allocate(xtemp(nr))
        xtemp = 0.d0
            
        do i = 1, nr
            xtemp(i) = b(i)
        end do
        
        ! Over-determined system: 
        lwork = 2*nc
        allocate(work(lwork))

        call dgels('N', nr, nc, 1, Atemp, nr, xtemp, nr, work, lwork, info)

        do i = 1, nc
            x(i) = xtemp(i)
        end do

    end if

end subroutine solve_system

subroutine crossproduct(v1, v2, v3)
    !! Computes cross product in a 3D Euclidean space

    ! Input / output data
    ! -------------------
    double precision, intent(in) :: v1, v2
    dimension :: v1(3), v2(3)

    double precision, intent(out) ::  v3
    dimension :: v3(3)

    ! Definition of cross product
    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

end subroutine crossproduct

subroutine kron_product_2vec(size_A, A, size_B, B, C, alpha)
    !! Evaluates kron product A x B = C (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output
    ! ---------------- 
    integer, intent(in) :: size_A, size_B
    double precision, intent(in) :: A(size_A), B(size_B)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: C(size_A*size_B)

    ! Local data
    ! ------------
    integer :: i1, i2, j, nb_tasks

    !$OMP PARALLEL PRIVATE(j)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, size_A*size_B/nb_tasks)
    do i1 = 1, size_A
        do i2 = 1, size_B
            j = i2 + (i1-1)*size_B
            C(j) = C(j) + alpha * A(i1) * B(i2)
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_2vec 

subroutine kron_product_3vec(size_A, A, size_B, B, size_C, C, D, alpha)
    !! Returns the result of A x B x C = D, where x is kronecker product

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: size_A,size_B, size_C
    double precision, intent(in) :: A, B, C
    dimension :: A(size_A), B(size_B), C(size_C)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: D
    dimension :: D(size_A*size_B*size_C)

    ! Local data
    ! -------------
    integer :: i1, i2, i3, j, nb_tasks

    !$OMP PARALLEL PRIVATE(j)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size_A*size_B*size_C/nb_tasks)
    do i1 = 1, size_A
        do i2 = 1, size_B
            do i3 = 1, size_C
                j = i3 + (i2-1)*size_C + (i1-1)*size_C*size_B
                D(j) = D(j) + alpha * A(i1) * B(i2) * C(i3)
            end do
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_3vec

subroutine spM_dot_dM(nrA, nnz_A, data_A, indi_A, indj_A, nrB, ncB, B, C)
    !! Returns the dot produt between a sparse matrix and a dense matrix
    !! It seems that it is not optimized

    implicit none
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nrA, nnz_A, nrB, ncB
    integer, intent(in) :: indi_A, indj_A
    dimension :: indi_A(nrA+1), indj_A(nnz_A)
    double precision, intent(in) :: data_A, B
    dimension :: data_A(nnz_A), B(nrB, ncB)

    double precision, intent(out) :: C
    dimension :: C(nrA, ncB)

    ! Local data
    ! --------------
    integer :: i, nnz_row, offset, table_nnz
    dimension :: table_nnz(nrA)
    integer, allocatable, dimension(:) :: indj_row
    double precision, allocatable, dimension(:) :: data_row

    ! Get number of non zeros in each row
    do i = 1, nrA
        table_nnz(i) = indi_A(i+1) - indi_A(i)
    end do

    ! Allocate indj and data with maximum value of non zeros
    nnz_row = maxval(table_nnz)
    allocate(indj_row(nnz_row), data_row(nnz_row))

    ! Compute result
    do i = 1, nrA
        ! Get information for each row
        nnz_row = table_nnz(i)
        offset = indi_A(i)

        ! Set the data for each row
        indj_row(1:nnz_row) = indj_A(offset:offset+nnz_row-1)
        data_row(1:nnz_row) = data_A(offset:offset+nnz_row-1)

        ! Do the dot product
        C(i, :) = matmul(data_row, B(indj_row, :))
    end do

end subroutine spM_dot_dM

! -------------
! Indices
! -------------
subroutine coo2csr(nr, nnz, a_in, indi_coo, indj_coo, a_out, indj_csr, indi_csr)
    !! Change COO format to CSR format

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nr, nnz
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)

    double precision, intent(out) :: a_out
    dimension :: a_out(nnz)
    integer, intent(out) :: indj_csr, indi_csr
    dimension :: indj_csr(nnz), indi_csr(nr+1)

    ! Local data
    ! -------------
    double precision :: x
    integer :: i, j, k, k0, iad

    ! Initialize
    indi_csr = 0

    ! Get number of non-zero terms for each row
    do  k = 1, nnz
        indi_csr(indi_coo(k)) = indi_csr(indi_coo(k)) + 1
    end do

    ! Get CSR format for i-indices 
    k = 1
    do j = 1, nr+1
        k0 = indi_csr(j)
        indi_csr(j) = k
        k = k + k0
    end do

    ! Sort COO non-zero values and update values of CSR j-indices
    do k = 1, nnz
        i = indi_coo(k)
        j = indj_coo(k)
        x = a_in(k)
        iad = indi_csr(i)
        a_out(iad) =  x
        indj_csr(iad) = j
        indi_csr(i) = iad + 1
    end do

    ! Update i-indices
    do  j = nr, 1, -1
        indi_csr(j+1) = indi_csr(j)
    end do
    indi_csr(1) = 1

end subroutine coo2csr

subroutine coo2matrix(nnz, indi_coo, indj_coo, a_in, nr, nc, A_out)
    !! Gives a dense matrix from a COO format
    !! Repeated positions are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: A_out
    dimension :: A_out(nr, nc)

    ! Local data 
    ! -------------
    integer :: i, j, k

    ! Initialize matrix values
    A_out = 0.d0

    ! Update values
    do k = 1, nnz 
        i = indi_coo(k)
        j = indj_coo(k)
        A_out(i, j) = A_out(i, j) + a_in(k)
    end do

end subroutine coo2matrix

subroutine csr2matrix(nnz, indi_csr, indj_csr, a_in, nr, nc, A_out)
    !! Gives a dense matrix from a CSR format
    !! Repeated positions are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: A_out
    dimension :: A_out(nr, nc)

    ! Local data
    ! -------------
    integer :: i, j, k
    integer :: nnz_col, offset

    ! Initialize matrix values
    A_out = 0.d0
    
    ! Update values
    do i = 1, nr
        nnz_col = indi_csr(i+1) - indi_csr(i)
        offset = indi_csr(i)
        do k = 1, nnz_col
            j = indj_csr(k+offset-1)
            A_out(i, j) = A_out(i, j) + a_in(k+offset-1)
        end do
    end do
    
end subroutine csr2matrix

subroutine matrix2csr(nr, nc, A_in, nnz, indi_csr, indj_csr)
    !! Returns CSR format from matrix but not the values
    !! Only for integers

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nr, nc, nnz
    integer, intent(in) :: A_in 
    dimension :: A_in(nr, nc)

    integer, intent(out) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)

    ! Local data
    ! -----------
    integer :: i, j, k, l

    ! Initialize
    k = 1
    indi_csr(1) = 1

    ! Update CSR format
    do i = 1, nr
        l = 0
        do j = 1, nc
            ! Save only values greater than zero
            if (abs(A_in(i, j)).gt.0) then
                indj_csr(k) = j
                k = k + 1
                l = l + 1
            end if
        end do
        indi_csr(i+1) = indi_csr(i) + l
    end do

end subroutine matrix2csr

subroutine csr2csc(nr, nc, nnz, a_in, indj_csr, indi_csr, a_out, indj_csc, indi_csc)
    !! Gets CSC format from CSR format. 
    !! (CSC format can be interpreted as the transpose)

    implicit none
    ! Input / output data 
    ! ----------------------
    integer, intent(in) :: nnz, nr, nc
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    integer, intent(out) :: indi_csc, indj_csc
    dimension :: indi_csc(nc+1), indj_csc(nnz)
    double precision, intent(out) :: a_out
    dimension :: a_out(nnz)

    ! Local data
    ! --------------
    integer :: indi_coo 
    dimension :: indi_coo(nnz)
    integer :: i, j, c

    ! We assume that csr is close to coo format. The only thing that change is indi
    c = 0
    do i = 1, nr
        do j = indi_csr(i), indi_csr(i+1) - 1
            c = c + 1
            indi_coo(c) = i
        end do
    end do

    ! Do COO to CSR format (inverting order to CSC)
    call coo2csr(nc, nnz, a_in, indj_csr, indi_coo, a_out, indj_csc, indi_csc)

end subroutine csr2csc

subroutine get_indexes_kron2_product(nr_A, nc_A, nnz_A, & 
                                indi_A, indj_A, &
                                nr_B, nc_B, nnz_B, &
                                indi_B, indj_B, &  
                                nr_C, nc_C, nnz_C, &
                                indi_C, indj_C)
    !! Returns indices of A x B = C (x : kronecker product)
    !! Where A and B are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! ----------------------
    integer, intent(in) ::  nr_A, nc_A, nnz_A, &
                            nr_B, nc_B, nnz_B, &
                            nnz_C
    integer, intent(in) :: indi_A, indj_A, indi_B, indj_B
    dimension ::    indi_A(nr_A+1), indj_A(nnz_A), &
                    indi_B(nr_B+1), indj_B(nnz_B)

    integer, intent(out) :: nr_C, nc_C
    integer, intent(out) :: indi_C, indj_C
    dimension :: indi_C(nr_A*nr_B+1), indj_C(nnz_C)

    ! Loca data
    ! -----------
    integer :: i1, i2, k, j1, j2, nb_tasks
    integer :: nnz_row_A, nnz_row_B, nnz_row_C
    integer :: count
    integer, allocatable, dimension(:) :: indj_C_temp

    ! Set new number of rows
    nr_C = nr_A * nr_B
    nc_C = nc_A * nc_B

    ! Set indices i in CSR format
    indi_C(1) = 1
    do i1 = 1, nr_A
        do i2 = 1, nr_B
            ! Find C's row position
            k = i2 + (i1 - 1)*nr_B

            ! Set number of non-zero elements of A's i-row  
            nnz_row_A = indi_A(i1+1) - indi_A(i1) 

            ! Set number of non-zero elements of B's j-row  
            nnz_row_B = indi_B(i2+1) - indi_B(i2)

            ! Set number of non-zero elements of C's k-row 
            nnz_row_C = nnz_row_A * nnz_row_B

            ! Update value 
            indi_C(k+1) = indi_C(k) + nnz_row_C 
        end do
    end do

    !$OMP PARALLEL PRIVATE(count,k,j1,j2,indj_C_temp) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nr_A*nr_B/nb_tasks)
    ! Set indices j in csr format
    do i1 = 1, nr_A
        do i2 = 1, nr_B
            ! Select row
            k = (i1 - 1)*nr_B + i2
            allocate(indj_C_temp(indi_C(k+1) - indi_C(k)))
            
            ! Get values of C's k-row
            count = 0
            do j1 = indi_A(i1), indi_A(i1+1) - 1        
                do j2 = indi_B(i2), indi_B(i2+1) - 1
                    count = count + 1
                    indj_C_temp(count) = (indj_A(j1) - 1)*nc_B + indj_B(j2)
                end do
            end do

            ! Update values
            indj_C(indi_C(k): indi_C(k+1)-1) = indj_C_temp
            deallocate(indj_C_temp)
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
end subroutine get_indexes_kron2_product

subroutine get_indexes_kron3_product(nr_A, nc_A, nnz_A, & 
                            indi_A, indj_A, &
                            nr_B, nc_B, nnz_B, &
                            indi_B, indj_B, &  
                            nr_C, nc_C, nnz_C, &
                            indi_C, indj_C, &
                            nr_D, nc_D, nnz_D, &
                            indi_D, indj_D)
    !! Returns indices of A x B x C = D (x : kronecker product)
    !! Where A, B and C are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! ----------------------
    integer, intent(in) ::  nr_A, nc_A, nnz_A, &
                            nr_B, nc_B, nnz_B, &
                            nr_C, nc_C, nnz_C, nnz_D
    integer, intent(in) :: indi_A, indj_A, indi_B, indj_B, indi_C, indj_C
    dimension ::    indi_A(nr_A+1), indj_A(nnz_A), &
                    indi_B(nr_B+1), indj_B(nnz_B), &
                    indi_C(nr_C+1), indj_C(nnz_C)

    integer, intent(out) :: nr_D, nc_D
    integer, intent(out) :: indi_D, indj_D
    dimension :: indi_D(nr_A*nr_B*nr_C+1), indj_D(nnz_D)

    ! Loca data
    ! -----------
    integer :: i1, i2, i3, k, j1, j2, j3, nb_tasks
    integer :: nnz_row_A, nnz_row_B, nnz_row_C, nnz_row_D
    integer :: count
    integer, allocatable, dimension(:) :: indj_D_temp

    ! Set new number of rows
    nr_D = nr_A * nr_B * nr_C
    nc_D = nc_A * nc_B * nc_D

    ! Set indices i in CSR format
    indi_D(1) = 1
    do i1 = 1, nr_A
        do i2 = 1, nr_B
            do i3 = 1, nr_C
                ! Find D's row position
                k = i3 + (i2-1)*nr_C + (i1 - 1)*nr_C*nr_B

                ! Set number of non-zero elements of A's i1-row  
                nnz_row_A = indi_A(i1+1) - indi_A(i1) 

                ! Set number of non-zero elements of B's i2-row  
                nnz_row_B = indi_B(i2+1) - indi_B(i2)

                ! Set number of non-zero elements of C's i3-row  
                nnz_row_C = indi_C(i3+1) - indi_C(i3)

                ! Set number of non-zero elements of D's k-row 
                nnz_row_D = nnz_row_A * nnz_row_B * nnz_row_C

                ! Update value 
                indi_D(k+1) = indi_D(k) + nnz_row_D
            end do
        end do
    end do

    !$OMP PARALLEL PRIVATE(count,k,j1,j2,j3,indj_D_temp) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_A*nr_B*nr_C/nb_tasks)
    ! Set indices j in csr format
    do i1 = 1, nr_A
        do i2 = 1, nr_B
            do i3 = 1, nr_C
                ! Select row
                k = i3 + (i2-1)*nr_C + (i1 - 1)*nr_C*nr_B
                allocate(indj_D_temp(indi_D(k+1) - indi_D(k)))
                
                ! Get values of D's k-row
                count = 0
                do j1 = indi_A(i1), indi_A(i1+1) - 1        
                    do j2 = indi_B(i2), indi_B(i2+1) - 1
                        do j3 = indi_C(i3), indi_C(i3+1) - 1
                            count = count + 1
                            indj_D_temp(count) = (indj_A(j1)-1)*nc_B*nc_C + (indj_B(j2)-1)*nc_C + indj_C(j3)
                        end do
                    end do
                end do

                ! Update values
                indj_D(indi_D(k): indi_D(k+1)-1) = indj_D_temp
                deallocate(indj_D_temp)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
end subroutine get_indexes_kron3_product