! ==========================
! module :: Linear algebra 
! author :: Joaquin Cornejo
! 
! This modules is intended to have all the linear algebra functions necessary for the good
! performance of the rest of functions.
! Some functions in this file calls functions defined in YETI or LAPACK libraries.
! ==========================

! --------------------
! Vector and matrices
! --------------------

subroutine find_interpolation_span(size_array, array, x, span, tol)
    !! Finds the interpolation span of the given knot. 
    !! Ex: Given the nodes {0, 0.5, 1} and x = 0.25, the interpolation span is 1

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_array
    double precision, intent(in) :: array, x, tol
    dimension :: array(size_array)

    integer, intent(out) :: span 

    ! Set first value of span
    span = 2
    
    ! Find span
    do while ((span.lt.size_array) &
            .and.((array(span)-x).le.tol))

        ! Update value 
        span = span + 1
    end do
    
    ! Set result
    span = span - 1 

end subroutine find_interpolation_span

subroutine linear_interpolation(nr, table, nnz, x, y, tol)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nnz
    double precision, intent(in) :: table, x, tol
    dimension :: table(nr, 2), x(nnz)

    double precision, intent(out) :: y
    dimension :: y(nnz)

    ! Local data
    ! ----------
    integer :: i, span
    double precision :: x1, x2, y1, y2, xi, yi

    ! Verify if table is well defined (values in ascending order)
    do i = 1, nr-1
        if (table(i+1,1)-table(i,1).le.tol) stop 'Table is not well - defined'
    end do

    do i = 1, nnz
        xi = x(i)
        if (xi.le.table(1,1)) then
            ! First case: x is less than the first point
            yi = table(1, 2)
        else if (xi.ge.table(nr,1)) then
            ! Second case: x is greatr than the last point
            yi = table(nr, 2)
        else
            ! Third case: x is within a knot defined in the table
            ! Find the span of x 
            call find_interpolation_span(nr, table(:, 1), xi, span, tol)

            ! Linear interpolation
            x1 = table(span, 1); x2 = table(span+1, 1)
            y1 = table(span, 2); y2 = table(span+1, 2)
            yi = y1 + (xi-x1)*(y2-y1)/(x2-x1)
        end if
        y(i) = yi
    end do

end subroutine linear_interpolation

subroutine diff_vector(times, nnz, vector_in, vector_out)
    !! Returns the difference between elements of array. 
    !! Ex. given a vector [0, 1, 3, 10] and times = 1, the output is [1, 2, 7, 0]. 
    !! If times = 2, the output is [1, 5, 0, 0]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz, times
    double precision, intent(in) :: vector_in
    dimension :: vector_in(nnz)

    double precision, intent(out) :: vector_out
    dimension :: vector_out(nnz)

    ! Local data
    ! ----------
    double precision :: vector_temp
    dimension :: vector_temp(nnz)
    integer :: i, j

    if (times.ge.nnz) stop 'Times can not be greater than the vector size'
    if (times.lt.0) stop 'Times can not be less than zero'

    ! Initialize
    vector_out = vector_in

    do j = 1, times

        ! Save temporal vector
        vector_temp = vector_out

        ! Get the difference
        do i = 1, nnz-j
            vector_out(i) = vector_temp(i+1) - vector_temp(i)
        end do

        ! Extra data is set to 0
        do i = nnz-j+1, nnz
            vector_out(i) = 0.d0
        end do

    end do
    
end subroutine diff_vector

subroutine find_unique_vector(nnz, vec, vec_unique)
    !! Gets the non-repreated values of a vector
    !! Ex: Given the vector [0, 1, 1, 2, 3, 3, 3], the unique vector is [0, 1, 2, 3, 0, 0, 0, 4]
    !! The last value is the number of unique elements

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz 
    double precision, intent(in) :: vec
    dimension :: vec(nnz)

    double precision, intent(out) :: vec_unique
    dimension :: vec_unique(nnz+1)

    ! Local data
    ! ----------
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
            !  Flag this value only if it has not already been flagged
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
    !! Evaluates n equidistant points given the first and last points 

    implicit none
    ! Input / output data
    ! -------------------
    double precision, intent(in) :: x0, xf 
    integer, intent(in) :: n 

    double precision, intent(out) :: array 
    dimension :: array(n)

    ! Local data
    ! ----------
    integer :: i
    double precision :: h 

    if (n.le.1) stop 'Size of array can not be 0 or less'

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

subroutine gemm_AWB(gemm_case, nrA, ncA, A, nrB, ncB, B, W, nrR, ncR, R)
    !! General matrix multiplication of the kind A diag(W) B
    !! 1: R = A.diag(W).BT
    !! 2: R = A.diag(W).BT, returns only the diagonal
    !! 3: R = AT.diag(W).B
    !! 4: R = AT.diag(W).B, returns only the diagonal
    !! Matrix A = (nb_rows_A, nb_columns_A)
    !! Array W = (*) it depends on the mode-multiplication
    !! Matrix B = (nb_rows_B, nb_columns_B)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: gemm_case, nrA, ncA, nrB, ncB, nrR, ncR
    double precision, intent(in) :: A, B, W
    dimension :: A(nrA, ncA), B(nrB, ncB), W(*)

    double precision, intent(out) :: R
    dimension :: R(nrR, ncR)

    ! Local data
    ! ----------
    integer :: i, j, k
    double precision :: s
    double precision, allocatable, dimension(:, :) :: MM

    ! Initiliaze 
    R = 0.d0

    if (gemm_case.eq.1) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB

        if (ncA.ne.ncB) stop 'Please verify size of matrices'
        allocate(MM(nrA, ncA))
        forall (i = 1 : nrA, j = 1 : ncA) 
            MM(i, j) = A(i, j) * W(j)
        end forall
        R = matmul(MM, transpose(B))

    else if (gemm_case.eq.2) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB
        ! Only diagonal

        if (ncA.ne.ncB) stop 'Please verify size of matrices'
        do i = 1, nrA
            s = 0.d0
            do k = 1, ncA
                s = s + A(i, k)*W(k)*B(i, k)
            end do
            R(i, i) = s
        end do

    else if (gemm_case.eq.3) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB

        if (nrA.ne.nrB) stop 'Please verify size of matrices'
        allocate(MM(nrB, ncB))
        forall (i = 1 : nrB, j = 1 : ncB) 
            MM(i, j) = W(i) * B(i, j) 
        end forall
        R = matmul(transpose(A), MM)

    else if (gemm_case.eq.4) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB
        ! Only diagonal

        if (nrA.ne.nrB) stop 'Please verify size of matrices'
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
    !! Solves linear equation system A.x = b
    !! It was adapted to solve any kind of linear systems (under, well and over determined systems).
    !! Determined system: LU decomposition with partial pivoting is used.
    !! Under-determined system: minimum norm solution method is used.
    !! Over-determined system: least squares solution method is used.
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
    ! ----------
    double precision :: Atemp
    dimension :: Atemp(nr, nc)
    double precision, allocatable, dimension(:) :: xtemp

    integer :: i, info, lwork, ipiv
    dimension :: ipiv(nr)
    double precision, allocatable, dimension(:) :: work

    ! Initialize 
    Atemp = A; x = 0.d0
    
    if (nr.le.nc) then 

        ! Initialize solution
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

        ! Save solution
        x = xtemp

    else if (nr.gt.nc) then 
        
        ! Initialize solution
        allocate(xtemp(nr))
        xtemp = b
        
        ! Over-determined system: 
        lwork = 2*nc
        allocate(work(lwork))
        call dgels('N', nr, nc, 1, Atemp, nr, xtemp, nr, work, lwork, info)

        ! Save solution
        do i = 1, nc
            x(i) = xtemp(i)
        end do

    end if

end subroutine solve_linear_system

subroutine inverse_matrix(nr, A, invA)
    !! Computes the inverse of a square matrix using LU decomposition
    
    implicit none
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nr
    double precision, intent(in) :: A
    dimension :: A(nr, nr)

    double precision, intent(out) :: invA
    dimension :: invA(nr, nr)

    ! Local data
    ! ----------
    double precision :: work
    integer :: ipiv, info
    dimension :: work(nr), ipiv(nr)

    ! Initialize
    invA = A

    ! Compute LU decomposition
    call dgetrf(nr, nr, invA, nr, ipiv, info)

    ! Compute inverse of A with LU decomposition
    call dgetri(nr, invA, nr, ipiv, work, nr, info)

end subroutine inverse_matrix

subroutine crossproduct(v1, v2, v3)
    !! Computes cross product in a 3D Euclidean space, v3 =  v1 x v2 (x is cross product)

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

subroutine kron_product_2vec(nnz_A, A, nnz_B, B, R, alpha)
    !! Computes kron product of 2 vectors and returns R = R + alpha*(A x B) (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz_A, nnz_B
    double precision, intent(in) :: A(nnz_A), B(nnz_B)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: R(nnz_A*nnz_B)

    ! Local data
    ! ----------
    integer :: iA, iB, genPos, nb_tasks

    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nnz_A*nnz_B/nb_tasks)
    do iA = 1, nnz_A
        do iB = 1, nnz_B
            genPos = iB + (iA-1)*nnz_B
            R(genPos) = R(genPos) + alpha * A(iA) * B(iB)
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_2vec 

subroutine kron_product_3vec(nnz_A, A, nnz_B, B, nnz_C, C, R, alpha)
    !! Computes the kron product of 3 vectors and returns R = R + alpha*(Au x Av x Aw) (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz_A, nnz_B, nnz_C
    double precision, intent(in) :: A, B, C
    dimension :: A(nnz_A), B(nnz_B), C(nnz_C)
    double precision, intent(in) :: alpha

    double precision, intent(inout) :: R
    dimension :: R(nnz_A*nnz_B*nnz_C)

    ! Local data
    ! ----------
    integer :: iA, iB, iC, genPos, nb_tasks

    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nnz_A*nnz_B*nnz_C/nb_tasks)
    do iA = 1, nnz_A
        do iB = 1, nnz_B
            do iC = 1, nnz_C
                genPos = iC + (iB-1)*nnz_C + (iA-1)*nnz_C*nnz_B
                R(genPos) = R(genPos) + alpha * A(iA) * B(iB) * C(iC)
            end do
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kron_product_3vec

subroutine polar_decomposition(A, Q, H, Sigma, onlyH, onlyDiag)
    !! Computes the polar decomposition of 3-by-3 matrix A = Q*H, where Q is orthogonal
    !! and H is symmetric positive semidefinite.
    !! It uses SDV: A = U S VT = (U VT) . (V S VT) then Q = U VT and H = V S VT

    implicit none
    ! Input / output data
    ! -------------------
    integer, parameter :: d = 3
    double precision, intent(in) :: A
    dimension :: A(d, d)
    logical, intent(in) :: onlyH, onlyDiag 

    double precision, intent(out) :: Q, H, Sigma
    dimension :: Q(d, d), H(d, d), Sigma(d)

    ! Local data
    ! ----------
    double precision :: U, VT, work
    dimension :: U(d, d), VT(d, d), work(5*d)
    integer :: info

    ! Compute singular value decomposition (SVD)
    if (onlyH) then
        ! We do not compute Q = U VT, then U = zeros
        call dgesvd('N', 'A', d, d, A, d, Sigma, U, d, VT, d, work, size(work), info)
        Q = 0.d0
    else 
        ! We do compute Q = U VT
        call dgesvd('A', 'A', d, d, A, d, Sigma, U, d, VT, d, work, size(work), info)
        Q = matmul(U, VT)
    end if

    ! Computes H
    if (onlyDiag) then
        call gemm_AWB(4, d, d, VT, d, d, VT, Sigma, d, d, H)
    else
        call gemm_AWB(3, d, d, VT, d, d, VT, Sigma, d, d, H)
    end if

end subroutine polar_decomposition

subroutine spMdotdV(nr, nc, nnz, indi, indj, A, vin, vout)
    !! Computes the dot product of sparse matrix with dense vector. It returns a dense vector
    !! Sparse matrix in CSR format

    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: A, vin
    dimension :: A(nnz), vin(nc)

    double precision, intent(out) :: vout
    dimension :: vout(nr)

    ! Local data
    ! ----------
    integer :: i , j, k
    double precision :: sum

    ! Initialize
    vout = 0.d0

    ! Compute the result at each row
    do i = 1, nr
        sum = 0.d0
        do j = indi(i), indi(i+1)-1
            k = indj(j)
            sum = sum + A(j)*vin(k)
        end do
        vout(i) = sum
    end do

end subroutine spMdotdV

! --------
! Indices
! --------
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
    ! ----------
    integer :: i, j, k, c, iad
    double precision :: x
    dimension :: x(nba)
    
    ! Initialize
    indi_csr = 0

    ! Get number of non-zero terms for each row
    do  k = 1, nnz
        indi_csr(indi_coo(k)) = indi_csr(indi_coo(k)) + 1
    end do

    ! Get CSR format for i-indices 
    k = 1
    do j = 1, nr+1
        c = indi_csr(j)
        indi_csr(j) = k
        k = k + c
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

    implicit none
    ! Input / output data 
    ! -------------------
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
    ! ----------
    integer :: indi_coo 
    dimension :: indi_coo(nnz)
    integer :: i, j, c

    ! We assume that csr is close to coo format. The only thing that change is the row indices (indi)
    c = 0
    do i = 1, nr
        do j = indi_csr(i), indi_csr(i+1) - 1
            c = c + 1
            indi_coo(c) = i
        end do
    end do

    ! Change COO format to CSR format (inverting order to CSC format)
    call coo2csr(nba, nc, nnz, a_in, indj_csr, indi_coo, a_out, indj_csc, indi_csc)

end subroutine csr2csc

subroutine coo2dense(nnz, indi_coo, indj_coo, a_in, nr, nc, a_out)
    !! Gets a dense matrix from COO format
    !! Repeated indices (i, j) are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: a_out
    dimension :: a_out(nr, nc)

    ! Local data 
    ! ----------
    integer :: i, j, k

    ! Initialize 
    a_out = 0.d0

    ! Update values
    do k = 1, nnz
        i = indi_coo(k)
        j = indj_coo(k)
        a_out(i, j) = a_out(i, j) + a_in(k)
    end do

end subroutine coo2dense

subroutine csr2dense(nnz, indi_csr, indj_csr, a_in, nr, nc, a_out)
    !! Gets a dense matrix from a CSR format
    !! Repeated positions are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz)

    double precision, intent(out) :: a_out
    dimension :: a_out(nr, nc)

    ! Local data
    ! -------------
    integer :: i, j, k
    integer :: nnz_col, offset

    ! Initialize matrix values
    a_out = 0.d0
    
    ! Update values
    do i = 1, nr
        nnz_col = indi_csr(i+1) - indi_csr(i)
        offset = indi_csr(i)
        do k = 1, nnz_col
            j = indj_csr(k+offset-1)
            a_out(i, j) = a_out(i, j) + a_in(k+offset-1)
        end do
    end do
    
end subroutine csr2dense

subroutine dense2csr(nr, nc, MM, nnz, indi_csr, indj_csr, data_csr)
    !! Gets CSR format from matrix 
    !! If nnz < 0, the it only computes the number of non zero values

    implicit none 
    ! Input / output data
    ! --------------------
    double precision, parameter :: tol = 1.d-15
    integer, intent(in) :: nr, nc
    double precision, intent(in) :: MM 
    dimension :: MM(nr, nc)

    integer, intent(inout) :: nnz
    integer, intent(out) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(*)
    double precision, intent(out) :: data_csr
    dimension :: data_csr(*)

    ! Local data
    ! ----------
    integer :: i, j, k, l

    if (nnz.le.0) then 
        ! Gets only the number of non zeros values of the matrix 
        nnz = 0
        do i = 1, nr
            do j = 1, nc
                if (abs(MM(i, j)).gt.tol) then
                    nnz = nnz + 1
                end if
            end do
        end do

    else
        ! Initialize
        k = 1
        indi_csr(1) = 1

        ! Update CSR format
        do i = 1, nr
            l = 0
            do j = 1, nc
                ! Save only values greater than the threshold
                if (abs(MM(i, j)).gt.tol) then
                    data_csr(k) = MM(i, j)
                    indj_csr(k) = j
                    k = k + 1
                    l = l + 1
                end if
            end do
            indi_csr(i+1) = indi_csr(i) + l
        end do

    end if

end subroutine dense2csr

subroutine erase_rows_csr(nr2er, rows2er, nr_in, nba, &
                        nnz_in, indi_in, indj_in, a_in, &
                        nnz_out, indi_out, indj_out, a_out)
    !! Erase rows from CSR format

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr2er, nba, nnz_in, nr_in
    integer, intent(inout) :: nnz_out
    integer, intent(in) :: rows2er, indi_in, indj_in
    dimension :: rows2er(nr2er), indi_in(nr_in+1), indj_in(nnz_in)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz_in, nba)

    integer, intent(out) :: indi_out, indj_out
    dimension :: indi_out(nr_in+1-nr2er), indj_out(nnz_out)
    double precision, intent(out) :: a_out
    dimension :: a_out(nnz_out, nba)

    ! Local data
    ! ----------
    integer :: i, j, k, c, nr_out, nnzrow
    integer, allocatable, dimension(:) :: rows2save

    ! Verify that rows to erase are consistant
    do i = 1, nr2er
        k = rows2er(i)
        if ((k.lt.1).or.(k.gt.nr_in)) then 
            stop 'Rows to erase are not well defined'
        end if
    end do

    ! Verify that rows to erase are different 
    do i = 1, nr2er
        c = 0
        k = rows2er(i)
        do j = 1, nr2er
            if (k.eq.rows2er(j)) then 
                c = c + 1
            end if
        end do
        if (c.ge.2) then 
            stop 'Rows to erase are repeated'
        end if
    end do

    if (nnz_out.le.0) then 
        ! Gets the number of non zeros values of output
        nnz_out = nnz_in

        do i = 1, nr2er
            j = rows2er(i)
            nnzrow = indi_in(j+1) - indi_in(j)
            nnz_out = nnz_out - nnzrow
        end do

    else
        ! Compute rows to save
        nr_out = nr_in - nr2er
        allocate(rows2save(nr_out))
        i = 1; j = 1
        do while ((j.le.nr_in).and.(i.le.nr_out))
            if (any(rows2er.eq.j)) then
                continue
            else
                rows2save(i) = j
                i = i + 1 
            end if
            j = j + 1
        end do

        ! Compute indices
        c = 1
        indi_out(1) = 1
        do i = 1, nr_out
            j = rows2save(i)
            nnzrow = indi_in(j+1) - indi_in(j)
            do k = indi_in(j), indi_in(j+1) - 1
                indj_out(c) = indj_in(k)
                a_out(c, :) = a_in(k, :)
                c = c + 1
            end do
            indi_out(i+1) = indi_out(i) + nnzrow
        end do

    end if

end subroutine erase_rows_csr

subroutine get_indexes_kron2_product(nr_A, nc_A, nnz_A, indi_A, indj_A, &
                                    nr_B, nc_B, nnz_B, indi_B, indj_B, &  
                                    nr_C, nc_C, nnz_C, indi_C, indj_C)
    !! Gets indices of A x B = C (x : kronecker product)
    !! Where A and B are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_A, nc_A, nnz_A, &
                            nr_B, nc_B, nnz_B, nnz_C
    integer, intent(in) :: indi_A, indj_A, indi_B, indj_B
    dimension ::    indi_A(nr_A+1), indj_A(nnz_A), &
                    indi_B(nr_B+1), indj_B(nnz_B)

    integer, intent(out) :: nr_C, nc_C
    integer, intent(out) :: indi_C, indj_C
    dimension :: indi_C(nr_A*nr_B+1), indj_C(nnz_C)

    ! Loca data
    ! ---------
    integer :: iA, iB, genPos, jA, jB, c, nb_tasks
    integer :: nnz_row_A, nnz_row_B, nnz_row_C
    integer, allocatable, dimension(:) :: indj_temp

    ! Set new number of rows
    nr_C = nr_A * nr_B
    nc_C = nc_A * nc_B

    ! Set indices i in CSR format
    indi_C(1) = 1
    do iA = 1, nr_A
        do iB = 1, nr_B
            ! Find C's row position
            genPos = iB + (iA - 1)*nr_B

            ! Set number of non-zero elements of A's i-row  
            nnz_row_A = indi_A(iA+1) - indi_A(iA) 

            ! Set number of non-zero elements of B's j-row  
            nnz_row_B = indi_B(iB+1) - indi_B(iB)

            ! Set number of non-zero elements of C's k-row 
            nnz_row_C = nnz_row_A * nnz_row_B

            ! Update value 
            indi_C(genPos+1) = indi_C(genPos) + nnz_row_C 
        end do
    end do

    !$OMP PARALLEL PRIVATE(c,genPos,jA,jB,indj_temp) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nr_A*nr_B/nb_tasks)
    ! Set indices j in csr format
    do iA = 1, nr_A
        do iB = 1, nr_B
            ! Select row
            genPos = (iA - 1)*nr_B + iB
            nnz_row_C = indi_C(genPos+1) - indi_C(genPos)
            allocate(indj_temp(nnz_row_C))
            
            ! Get values of C's row
            c = 0
            do jA = indi_A(iA), indi_A(iA+1) - 1        
                do jB = indi_B(iB), indi_B(iB+1) - 1
                    c = c + 1
                    indj_temp(c) = (indj_A(jA) - 1)*nc_B + indj_B(jB)
                end do
            end do

            ! Update values
            indj_C(indi_C(genPos):indi_C(genPos+1)-1) = indj_temp
            deallocate(indj_temp)
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
end subroutine get_indexes_kron2_product

subroutine get_indexes_kron3_product(nr_A, nc_A, nnz_A, indi_A, indj_A, &
                                    nr_B, nc_B, nnz_B, indi_B, indj_B, &  
                                    nr_C, nc_C, nnz_C, indi_C, indj_C, &
                                    nr_D, nc_D, nnz_D, indi_D, indj_D)
    !! Returns indices of A x B x C = D (x : kronecker product)
    !! Where A, B and C are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
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
    ! ---------
    integer :: iA, iB, iC, genPos, jA, jB, jC, c, nb_tasks
    integer :: nnz_row_A, nnz_row_B, nnz_row_C, nnz_row_D
    integer, allocatable, dimension(:) :: indj_temp

    ! Set new number of rows
    nr_D = nr_A * nr_B * nr_C
    nc_D = nc_A * nc_B * nc_D

    ! Set indices i in CSR format
    indi_D(1) = 1
    do iA = 1, nr_A
        do iB = 1, nr_B
            do iC = 1, nr_C
                ! Find D's row position
                genPos = iC + (iB-1)*nr_C + (iA - 1)*nr_C*nr_B

                ! Set number of non-zero elements of A's i1-row  
                nnz_row_A = indi_A(iA+1) - indi_A(iA) 

                ! Set number of non-zero elements of B's i2-row  
                nnz_row_B = indi_B(iB+1) - indi_B(iB)

                ! Set number of non-zero elements of C's i3-row  
                nnz_row_C = indi_C(iC+1) - indi_C(iC)

                ! Set number of non-zero elements of D's k-row 
                nnz_row_D = nnz_row_A * nnz_row_B * nnz_row_C

                ! Update value 
                indi_D(genPos+1) = indi_D(genPos) + nnz_row_D
            end do
        end do
    end do

    !$OMP PARALLEL PRIVATE(c,genPos,jA,jB,jC,indj_temp) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_A*nr_B*nr_C/nb_tasks)
    ! Set indices j in csr format
    do iA = 1, nr_A
        do iB = 1, nr_B
            do iC = 1, nr_C
                ! Select row
                genPos = iC + (iB-1)*nr_C + (iA - 1)*nr_C*nr_B
                nnz_row_D = indi_D(genPos+1) - indi_D(genPos)
                allocate(indj_temp(nnz_row_D))
                
                ! Get values of D's row
                c = 0
                do jA = indi_A(iA), indi_A(iA+1) - 1        
                    do jB = indi_B(iB), indi_B(iB+1) - 1
                        do jC = indi_C(iC), indi_C(iC+1) - 1
                            c = c + 1
                            indj_temp(c) = (indj_A(jA)-1)*nc_B*nc_C &
                                            + (indj_B(jB)-1)*nc_C + indj_C(jC)
                        end do
                    end do
                end do

                ! Update values
                indj_D(indi_D(genPos):indi_D(genPos+1)-1) = indj_temp
                deallocate(indj_temp)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 
    
end subroutine get_indexes_kron3_product
