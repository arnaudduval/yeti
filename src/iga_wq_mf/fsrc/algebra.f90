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

subroutine product_AWB(nb_rows_A, nb_rows_B, nb_cols, A, W, B, AWB)
    !! Matrix multiplication type: A.diag(W).transpose(B)
    !! Matrix A = (nb_rows_A, nb_columns)
    !! Array W = (nb_columns)
    !! Matrix B = (nb_rows_B, nb_columns)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_rows_A, nb_rows_B, nb_cols
    double precision, intent(in) :: A, B, W
    dimension ::    A(nb_rows_A, nb_cols), &
                    B(nb_rows_B, nb_cols), &
                    W(nb_cols)

    double precision, intent(out) :: AWB
    dimension :: AWB(nb_rows_A, nb_rows_B)

    ! Local data
    ! -------------
    double precision :: AW
    dimension :: AW(nb_rows_A, nb_cols)
    integer :: i, j

    ! Evaluate AW = A * diag(W)
    do j = 1, nb_cols
        do i = 1, nb_rows_A
            AW(i, j) = A(i, j) * W(j)
        end do
    end do

    ! Evaluate AB = A * B.T
    AWB = matmul(AW, transpose(B))

end subroutine product_AWB

subroutine solve_system(nb_rows, nb_cols, A, b, x)
    !! Solves equation system A.x = b
    !! Matrix A = (nb_rows, nb_columns)
    !! Vector b = (nb_rows)
    !! Vector x = (nb_cols)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_rows, nb_cols
    double precision, intent(in) :: A, b
    dimension :: A(nb_rows, nb_cols), b(nb_rows)

    double precision, intent(out) :: x
    dimension :: x(nb_cols)

    ! Local data
    ! ------------------
    integer :: i
    double precision :: A_copy
    dimension :: A_copy(nb_rows, nb_cols)

    double precision :: xcase1, xcase2
    dimension :: xcase1(nb_cols), xcase2(nb_rows)

    ! Lapack
    integer :: IPIV, INFO, LWORK
    dimension :: IPIV(nb_rows)
    double precision, allocatable :: WORK(:)

    ! Set a copy of matrix A 
    A_copy = A

    ! Set true solution x and intermediate solutions
    x = 0.d0
    xcase1 = 0.d0
    xcase2 = 0.d0

    if (nb_rows.le.nb_cols) then ! Case 1
        
        do i = 1, nb_rows
            xcase1(i) = b(i)
        end do

        if (nb_rows.eq.nb_cols) then

            ! Determmined system: 
            call dgesv(nb_rows, 1, A_copy, nb_rows, IPIV, xcase1, nb_rows, INFO)

        elseif (nb_rows.lt.nb_cols) then

            ! Under-determined system: 
            LWORK = 2 * nb_rows
            allocate(WORK(LWORK))
            call dgels('N', nb_rows, nb_cols, 1, A_copy, nb_rows, xcase1, nb_cols, WORK, LWORK, INFO)
                
        end if

        do i = 1, nb_cols
            x(i) = xcase1(i)
        end do

    else if (nb_rows.gt.nb_cols) then ! Case 2
            
        do i = 1, nb_rows
            xcase2(i) = b(i)
        end do
        
        ! Over-determined system: 
        LWORK = 2 * nb_cols
        allocate(WORK(LWORK))

        call dgels('N', nb_rows, nb_cols, 1, A_copy, nb_rows, xcase2, nb_rows, WORK, LWORK, INFO)

        do i = 1, nb_cols
            x(i) = xcase2(i)
        end do

    end if

end subroutine solve_system

subroutine vector_kron_vector(size_A, A, size_B, B, C)
    !! Evaluates kron product A x B (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output
    ! ---------------- 
    integer, intent(in) :: size_A, size_B
    double precision, intent(in) :: A(size_A), B(size_B)

    double precision, intent(out) :: C(size_A*size_B)

    ! Local data
    ! ------------
    integer :: i1, i2, indj, nb_tasks

    ! Initialize vector
    C = 0.d0

    ! Complete result
    !$OMP PARALLEL PRIVATE(i1, i2, indj)
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size_A*size_B/nb_tasks)
    do i1 = 1, size_A
        do i2 = 1, size_B
            indj = (i1-1)*size_B + i2
            C(indj) = A(i1) * B(i2)
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine vector_kron_vector 

subroutine kron_product_3vec(size_A, A, size_B, B, size_C, C, D, alpha)
    !! Returns the result of A x B x C, where x is kronecker product

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
    integer :: i, nb_tasks
    double precision, allocatable, dimension(:) :: AB, Dtemp

    ! Compute A x B
    allocate(AB(size_A*size_B))
    call vector_kron_vector(size_A, A, size_B, B, AB)

    ! Compute (A x B) x C
    allocate(Dtemp(size_A*size_B*size_C))
    call vector_kron_vector(size(AB), AB, size_C, C, Dtemp)
    deallocate(AB)

    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, size(Dtemp)/nb_tasks)
    do i = 1, size(Dtemp)
        D(i) = D(i) + alpha*Dtemp(i)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_3vec

! -------------
! Indexes
! -------------
subroutine coo2csr(nb_rows, nnz, a_in, indi_coo, indj_coo, a_out, indj_csr, indi_csr)
    !! Change COO format to CSR format
    !! Algorithm adapted from f70 to f90 (sparskit library)

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_rows, nnz
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)

    double precision, intent(out) :: a_out
    dimension :: a_out(nnz)
    integer, intent(out) :: indj_csr, indi_csr
    dimension :: indj_csr(nnz), indi_csr(nb_rows+1)

    ! Local data
    ! -------------
    double precision :: x
    integer :: i, j, k
    integer :: k0, iad

    do k = 1, nb_rows+1
        indi_csr(k) = 0
    end do

    do  k = 1, nnz
        indi_csr(indi_coo(k)) = indi_csr(indi_coo(k)) + 1
    end do

    k = 1
    do j = 1, nb_rows+1
        k0 = indi_csr(j)
        indi_csr(j) = k
        k = k + k0
    end do

    do k = 1, nnz
        i = indi_coo(k)
        j = indj_coo(k)
        x = a_in(k)
        iad = indi_csr(i)
        a_out(iad) =  x
        indj_csr(iad) = j
        indi_csr(i) = iad + 1
    end do

    do  j = nb_rows, 1, -1
        indi_csr(j+1) = indi_csr(j)
    end do

    indi_csr(1) = 1

end subroutine coo2csr

subroutine coo2matrix(nnz, indi_coo, indj_coo, a_in, nb_rows, nb_cols, A_out)
    !! Gives a dense matrix from a COO format
    !! Repeated positions are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nb_rows, nb_cols 
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: A_out
    dimension :: A_out(nb_rows, nb_cols)

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

subroutine csr2matrix(nnz, indi_csr, indj_csr, a_in, nb_rows, nb_cols, A_out)
    !! Gives a dense matrix from a CSR format
    !! Repeated positions are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nb_rows, nb_cols 
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nb_rows+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: A_out
    dimension :: A_out(nb_rows, nb_cols)

    ! Local data
    ! -------------
    integer :: i, j, k
    integer :: nnz_col, offset

    ! Initialize matrix values
    A_out = 0.d0
    
    ! Update values
    do i = 1, nb_rows
        nnz_col = indi_csr(i+1) - indi_csr(i)
        offset = indi_csr(i)
        do k = 1, nnz_col
            j = indj_csr(k+offset-1)
            A_out(i,j) = A_out(i,j) + a_in(k+offset-1)
        end do
    end do
    
end subroutine csr2matrix

subroutine matrix2csr(nb_rows, nb_cols, A_in, nnz, indi_csr, indj_csr)

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_rows, nb_cols, nnz
    integer, intent(in) :: A_in 
    dimension :: A_in(nb_rows, nb_cols)

    integer, intent(out) :: indi_csr, indj_csr
    dimension :: indi_csr(nb_rows+1), indj_csr(nnz)

    ! Local data
    ! -----------
    integer :: i, j, k, l

    ! Initialize
    k = 1
    indi_csr(1) = 1

    ! Update CSR format
    do i = 1, nb_rows
        l = 0
        do j = 1, nb_cols
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

subroutine csr2csc(nb_rows, nb_cols, nnz, a_in, indj_csr, indi_csr, a_out, indj_csc, indi_csc)
    !! Gets CSC format from CSR format. 
    !! (CSC format can be interpreted as the transpose)

    implicit none
    ! Input / output data 
    ! ----------------------
    integer, intent(in) :: nnz, nb_rows, nb_cols
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nb_rows+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    integer, intent(out) :: indi_csc, indj_csc
    dimension :: indi_csc(nb_cols+1), indj_csc(nnz)
    double precision, intent(out) :: a_out
    dimension :: a_out(nnz)

    ! Local data
    ! --------------
    integer :: indi_coo 
    dimension :: indi_coo(nnz)
    integer :: i, j, c

    ! We assume that csr is close to coo format. The only thing that change is indi
    c = 0
    do i = 1, nb_rows
        do j = indi_csr(i), indi_csr(i+1) - 1
            c = c + 1
            indi_coo(c) = i
        end do
    end do

    ! Do COO to CSR format (inverting order to CSC)
    call coo2csr(nb_cols, nnz, a_in, indj_csr, indi_coo, a_out, indj_csc, indi_csc)

end subroutine csr2csc

subroutine get_indexes_kron_product(nb_rows_A, nb_cols_A, nnz_A, & 
                                indi_A, indj_A, &
                                nb_rows_B, nb_cols_B, nnz_B, &
                                indi_B, indj_B, &  
                                nb_rows_C, nb_cols_C, nnz_C, &
                                indi_C, indj_C)
    !! Returns indexes of A x B = C (x : kronecker product)
    !! Where A and B are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! ----------------------
    integer, intent(in) ::  nb_rows_A, nb_cols_A, nnz_A, &
                            nb_rows_B, nb_cols_B, nnz_B, &
                            nnz_C
    integer, intent(in) :: indi_A, indj_A, indi_B, indj_B
    dimension ::    indi_A(nb_rows_A+1), indj_A(nnz_A), &
                    indi_B(nb_rows_B+1), indj_B(nnz_B)

    integer, intent(out) :: nb_rows_C, nb_cols_C
    integer, intent(out) :: indi_C, indj_C
    dimension :: indi_C(nb_rows_A*nb_rows_B+1), indj_C(nnz_C)

    ! Loca data
    ! -----------
    integer :: i, j, k, m, n, nb_tasks
    integer :: nnz_row_A, nnz_row_B, nnz_row_C
    integer :: count
    integer, allocatable, dimension(:) :: indj_C_temp

    ! Set new number of rows
    nb_rows_C = nb_rows_A * nb_rows_B
    nb_cols_C = nb_cols_A * nb_cols_B

    ! Set indexes i in CSR format
    indi_C(1) = 1
    do i = 1, nb_rows_A
        do j = 1, nb_rows_B
            ! Find C's row position
            k = (i - 1)*nb_rows_B + j

            ! Set number of non-zero elements of A's i-row  
            nnz_row_A = indi_A(i+1) - indi_A(i) 

            ! Set number of non-zero elements of B's j-row  
            nnz_row_B = indi_B(j+1) - indi_B(j)

            ! Set number of non-zero elements of C's k-row 
            nnz_row_C = nnz_row_A * nnz_row_B

            ! Update value 
            indi_C(k+1) = indi_C(k) + nnz_row_C 
        end do
    end do

    !$OMP PARALLEL PRIVATE(count,k,m,n,indj_C_temp) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nb_rows_A*nb_rows_B/nb_tasks)
    ! Set indexes j in csr format
    do i = 1, nb_rows_A
        do j = 1, nb_rows_B
            ! Select row
            k = (i - 1)*nb_rows_B + j
            allocate(indj_C_temp(indi_C(k+1) - indi_C(k)))
            
            ! Get values of C's k-row
            count = 0
            do m = indi_A(i), indi_A(i+1) - 1        
                do n = indi_B(j), indi_B(j+1) - 1
                    count = count + 1
                    indj_C_temp(count) = (indj_A(m) - 1)*nb_cols_B + indj_B(n)
                end do
            end do

            ! Update values
            indj_C(indi_C(k): indi_C(k+1)-1) = indj_C_temp
            deallocate(indj_C_temp)
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
end subroutine get_indexes_kron_product