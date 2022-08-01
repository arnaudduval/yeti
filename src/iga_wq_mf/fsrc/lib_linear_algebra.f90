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

subroutine scale_vector(nnz, factor_up, factor_down, vector)
    !! Scaling in fast diagonalization

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz
    double precision, intent(in) :: factor_up, factor_down
    dimension :: factor_up(nnz), factor_down(nnz)

    double precision, intent(inout) :: vector
    dimension :: vector(nnz)

    ! Local data
    ! -------------
    integer :: i, nb_tasks

    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks)
    do i = 1, nnz
        vector(i) = sqrt(factor_up(i)/factor_down(i)) * vector(i) 
    end do  
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine scale_vector

subroutine diff_vector(times, nnz, vector_in, vector_out)
    !! Returns the difference between elements of array

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nnz, times
    double precision, intent(in) :: vector_in
    dimension :: vector_in(nnz)

    double precision, intent(out) :: vector_out
    dimension :: vector_out(nnz)

    ! Local data
    ! -----------
    double precision :: vector_temp
    dimension :: vector_temp(nnz)
    integer :: i, j

    ! Initialize
    vector_out = vector_in
    vector_temp = vector_out

    do j = 1, times

        do i = 1, nnz-j
            vector_out(i) = vector_temp(i+1) - vector_temp(i)
        end do

        do i = nnz-j+1, nnz
            vector_out(i) = 0.d0
        end do

        vector_temp = vector_out
    end do
    
end subroutine diff_vector

subroutine find_unique_vector(nnz, vec, vec_unique)
    !! Returns the non-repreated values of a vector

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz 
    double precision, intent(in) :: vec
    dimension :: vec(nnz)

    double precision, intent(out) :: vec_unique
    dimension :: vec_unique(nnz+1)

    ! Local data
    ! --------------
    integer :: i, num, nnz_nr
    logical, dimension(nnz) :: mask

    ! Initialize
    vec_unique = 0.d0

    ! Define mask 
    mask = .false.

    do i = 1, nnz

        ! Count the number of occurrences of this element:  
        num = count(vec(i).eq.vec)
    
        if (num.eq.1) then  
            ! There is only one, flag it:  
            mask(i) = .true.  
        else  
            !  Flag this value only if it hasn't already been flagged:  
            if (.not. any(vec(i).eq.vec .and. mask)) then
                mask(i) = .true.  
            end if
        end if
    
    end do

    ! Return only flagged elements
    nnz_nr = count(mask)
    vec_unique(1:nnz_nr) = pack(vec, mask)
    vec_unique(nnz+1) = dble(nnz_nr)

end subroutine find_unique_vector

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

subroutine gemm_AWB(mode, nrA, ncA, A, nrB, ncB, B, W, nrR, ncR, R)
    !! General matrix multiplication M1 Diag M2
    !! 1: R = A.diag(W).BT
    !! 2: R = A.diag(W).BT only diagonal
    !! 3: R = AT.diag(W).B
    !! 4: R = AT.diag(W).B only diagonal
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
    double precision, allocatable, dimension(:, :) :: TT
    integer :: i, j, k
    double precision :: s

    ! Initiliaze 
    R = 0.d0

    if (mode.eq.1) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB

        allocate(TT(nrA, ncA))
        forall (i = 1 : nrA, j = 1 : ncA) 
            TT(i, j) = A(i, j) * W(j)
        end forall
        R = matmul(TT, transpose(B))

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

        allocate(TT(nrB, ncB))
        forall (i = 1 : nrB, j = 1 : ncB) 
            TT(i, j) = W(i) * B(i, j) 
        end forall
        R = matmul(transpose(A), TT)

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

end subroutine gemm_AWB

subroutine solve_linear_system(nr, nc, A, b, x)
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

end subroutine solve_linear_system

subroutine inverse_matrix(nnz, A)
    ! Computes the inverse of a matrix using LU decomposition
    
    implicit none
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nnz
    double precision, intent(inout) :: A
    dimension :: A(nnz, nnz)

    ! Local data
    ! ----------
    double precision :: work
    integer :: ipiv, info
    dimension :: work(nnz), ipiv(nnz)

    ! Compute LU decomposition
    call dgetrf(nnz, nnz, A, nnz, ipiv, info)

    ! Compute inverse of A with LU decomposition
    call dgetri(nnz, A, nnz, ipiv, work, nnz, info)

end subroutine inverse_matrix

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

subroutine kron_product_2vec(nnz_A1, A1, nnz_A2, A2, R, alpha)
    !! Evaluates kron product R = alpha*R + A1 x A2 (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output
    ! ---------------- 
    integer, intent(in) :: nnz_A1, nnz_A2
    double precision, intent(in) :: A1(nnz_A1), A2(nnz_A2)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: R(nnz_A1*nnz_A2)

    ! Local data
    ! ------------
    integer :: i1, i2, j, nb_tasks

    !$OMP PARALLEL PRIVATE(j)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nnz_A1*nnz_A2/nb_tasks)
    do i1 = 1, nnz_A1
        do i2 = 1, nnz_A2
            j = i2 + (i1-1)*nnz_A2
            R(j) = R(j) + alpha * A1(i1) * A2(i2)
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_2vec 

subroutine kron_product_3vec(nnz_A1, A1, nnz_A2, A2, nnz_A3, A3, R, alpha)
    !! Returns the result of R = alpha*R + A1 x A2 x A3, where x is kronecker product

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz_A1,nnz_A2, nnz_A3
    double precision, intent(in) :: A1, A2, A3
    dimension :: A1(nnz_A1), A2(nnz_A2), A3(nnz_A3)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: R
    dimension :: R(nnz_A1*nnz_A2*nnz_A3)

    ! Local data
    ! -------------
    integer :: i1, i2, i3, j, nb_tasks

    !$OMP PARALLEL PRIVATE(j)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nnz_A1*nnz_A2*nnz_A3/nb_tasks)
    do i1 = 1, nnz_A1
        do i2 = 1, nnz_A2
            do i3 = 1, nnz_A3
                j = i3 + (i2-1)*nnz_A3 + (i1-1)*nnz_A3*nnz_A2
                R(j) = R(j) + alpha * A1(i1) * A2(i2) * A3(i3)
            end do
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_3vec

subroutine polar_decomposition(A, Q, H, Sigma, onlyH, onlyDiag)
    !! Polar decompoition of 3-by-3 matrix A = Q*H, where Q is orthogonal
    !! and H is symmetric positive semidefinite.
    !! It uses SDV: A = U S VT = U VT . V S VT then Q = U VT and H = V S VT

    implicit none
    ! Input / output data
    ! --------------------
    double precision, intent(in) :: A
    integer, intent(in) :: onlyH, onlyDiag 
    dimension :: A(3, 3)

    double precision, intent(out) :: Q, H, Sigma
    dimension :: Q(3, 3), H(3, 3), Sigma(3)

    ! Local data
    ! --------------
    double precision :: U, VT, work
    dimension :: work(15), U(3, 3), VT(3, 3)
    integer :: info

    ! Compute Singular value decomposition
    if (onlyH.eq.1) then
        ! We do not compute Q = U VT, then U = zeros
        call dgesvd('N', 'A', 3, 3, A, 3, Sigma, U, 3, VT, 3, work, 15, info)
        Q = 0.d0
    else 
        ! We do compute Q = U VT
        call dgesvd('A', 'A', 3, 3, A, 3, Sigma, U, 3, VT, 3, work, 15, info)
        Q = matmul(U, VT)
    end if

    if (onlyDiag.eq.1) then
        call gemm_AWB(4, 3, 3, VT, 3, 3, VT, Sigma, 3, 3, H)
    else
        call gemm_AWB(3, 3, 3, VT, 3, 3, VT, Sigma, 3, 3, H)
    end if

end subroutine polar_decomposition

! -------------
! Indices
! -------------
subroutine coo2csr(nba, nr, nnz, a_in, indi_coo, indj_coo, a_out, indj_csr, indi_csr)
    !! Change COO format to CSR format

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nba, nr, nnz
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz, nba)
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)

    double precision, intent(out) :: a_out
    dimension :: a_out(nnz, nba)
    integer, intent(out) :: indj_csr, indi_csr
    dimension :: indj_csr(nnz), indi_csr(nr+1)

    ! Local data
    ! -------------
    double precision :: x
    dimension :: x(nba)
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
        x = a_in(k, :)
        iad = indi_csr(i)
        a_out(iad, :) =  x
        indj_csr(iad) = j
        indi_csr(i) = iad + 1
    end do

    ! Update i-indices
    do  j = nr, 1, -1
        indi_csr(j+1) = indi_csr(j)
    end do
    indi_csr(1) = 1

end subroutine coo2csr

subroutine csr2csc(nba, nr, nc, nnz, a_in, indj_csr, indi_csr, a_out, indj_csc, indi_csc)
    !! Gets CSC format from CSR format. 
    !! (CSC format can be interpreted as the transpose)

    implicit none
    ! Input / output data 
    ! ----------------------
    integer, intent(in) :: nba, nnz, nr, nc
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz, nba)

    integer, intent(out) :: indi_csc, indj_csc
    dimension :: indi_csc(nc+1), indj_csc(nnz)
    double precision, intent(out) :: a_out
    dimension :: a_out(nnz, nba)

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
    call coo2csr(nba, nc, nnz, a_in, indj_csr, indi_coo, a_out, indj_csc, indi_csc)

end subroutine csr2csc

subroutine coo2dense(nnz, indi_coo, indj_coo, a_in, nr, nc, A_out)
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

end subroutine coo2dense

subroutine csr2dense(nnz, indi_csr, indj_csr, a_in, nr, nc, A_out)
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
    
end subroutine csr2dense

subroutine dense2csr(nr, nc, A_in, nnz, indi_csr, indj_csr)
    !! Returns CSR format from matrix but not the values
    !! Only for integers

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nr, nc, nnz
    double precision, intent(in) :: A_in 
    dimension :: A_in(nr, nc)

    integer, intent(out) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)

    ! Local data
    ! -----------
    integer :: i, j, k, l

    ! Initialize
    k = 1
    indi_csr(1) = 1
    indj_csr = 0

    ! Update CSR format
    do i = 1, nr
        l = 0
        do j = 1, nc
            ! Save only values greater than zero
            if (abs(A_in(i, j)).gt.tol) then
                indj_csr(k) = j
                k = k + 1
                l = l + 1
            end if
        end do
        indi_csr(i+1) = indi_csr(i) + l
    end do

end subroutine dense2csr

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
                            indj_D_temp(count) = (indj_A(j1)-1)*nc_B*nc_C &
                                                + (indj_B(j2)-1)*nc_C + indj_C(j3)
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
