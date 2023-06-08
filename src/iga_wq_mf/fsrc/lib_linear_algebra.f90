! ==============================================
! This modules is intended to have all the linear algebra functions necessary for the good
! performance of the rest of the files.
! Some functions in this file calls functions defined in YETI (see documentation) or LAPACK libraries.
! ==============================================

! --------------------
! Vector and matrices
! --------------------

subroutine diff_array(repeat, nnz, array_in, array_out)
    !! Returns the difference between the elements of an array 
    !! Ex. given an array [0, 1, 3, 10] and repeat = 1, the output is [1, 2, 7, 0]. 
    !! If repeat = 2, the output is [1, 5, 0, 0]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz, repeat
    double precision, intent(in) :: array_in
    dimension :: array_in(nnz)

    double precision, intent(out) :: array_out
    dimension :: array_out(nnz)

    ! Local data
    ! ----------
    double precision :: array_temp
    dimension :: array_temp(nnz)
    integer :: i, j

    if (repeat.ge.nnz) stop 'It can not be greater than the array size'
    if (repeat.lt.0) stop 'It can not be less than zero'

    array_out = array_in

    do j = 1, repeat

        array_temp = array_out

        ! Get the difference
        do i = 1, nnz-j
            array_out(i) = array_temp(i+1) - array_temp(i)
        end do

        ! Extra data is set to 0
        do i = nnz-j+1, nnz
            array_out(i) = 0.d0
        end do

    end do
    
end subroutine diff_array

subroutine find_unique_array(nnz, array, array_unique)
    !! Gets the non-repreated values of an array
    !! Ex: Given the array [0, 1, 1, 2, 3, 3, 3], the unique array is [0, 1, 2, 3, 0, 0, 0, 4]
    !! The last value is the number of unique elements

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz 
    double precision, intent(in) :: array
    dimension :: array(nnz)

    double precision, intent(out) :: array_unique
    dimension :: array_unique(nnz+1)

    ! Local data
    ! ----------
    integer :: i, c, nnz_unique
    logical :: mask
    dimension :: mask(nnz)

    array_unique = 0.d0
    mask = .false.

    do i = 1, nnz

        ! Count the number of occurrences of this element
        c = count(array(i).eq.array)
    
        if (c.eq.1) then  
            ! There is only one, flag it:  
            mask(i) = .true.  
        else  
            !  Flag this value only if it has not already been flagged
            if (.not. any(array(i).eq.array .and. mask)) then
                mask(i) = .true.  
            end if
        end if
    
    end do

    ! Return only flagged elements
    nnz_unique = count(mask)
    array_unique(1:nnz_unique) = pack(array, mask)
    array_unique(nnz+1) = dble(nnz_unique)

end subroutine find_unique_array

subroutine linspace(x0, xf, nnz, array) 
    !! Evaluates n equidistant points given the first and last points 

    implicit none
    ! Input / output data
    ! -------------------
    double precision, intent(in) :: x0, xf 
    integer, intent(in) :: nnz 

    double precision, intent(out) :: array 
    dimension :: array(nnz)

    ! Local data
    ! ----------
    integer :: i
    double precision :: h 

    if (nnz.le.1) stop 'Size of array can not be 1 or less'

    array(1) = x0; array(nnz) = xf

    h = (xf - x0)/dble(nnz - 1)
    
    do i = 2, nnz - 1
        array(i) = x0 + dble(i-1)*h
    end do

end subroutine linspace

subroutine gemm_AWB(type, nrA, ncA, A, nrB, ncB, B, W, nrR, ncR, R)
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
    integer, intent(in) :: type, nrA, ncA, nrB, ncB, nrR, ncR
    double precision, intent(in) :: A, B, W
    dimension :: A(nrA, ncA), B(nrB, ncB), W(*)

    double precision, intent(out) :: R
    dimension :: R(nrR, ncR)

    ! Local data
    ! ----------
    integer :: i, j, k
    double precision :: s
    double precision, allocatable, dimension(:, :) :: MM

    R = 0.d0

    if (type.eq.1) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB

        if (ncA.ne.ncB) stop 'Please verify size of matrices'
        allocate(MM(nrA, ncA))
        forall (i = 1 : nrA, j = 1 : ncA) 
            MM(i, j) = A(i, j) * W(j)
        end forall
        R = matmul(MM, transpose(B))

    else if (type.eq.2) then 
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

    else if (type.eq.3) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB

        if (nrA.ne.nrB) stop 'Please verify size of matrices'
        allocate(MM(nrB, ncB))
        forall (i = 1 : nrB, j = 1 : ncB) 
            MM(i, j) = W(i) * B(i, j) 
        end forall
        R = matmul(transpose(A), MM)

    else if (type.eq.4) then 
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
    !! Solves linear equation system A x = b
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
    double precision :: Acopy
    dimension :: Acopy(nr, nc)
    double precision, allocatable, dimension(:) :: xsol

    integer :: i, info, ipiv
    dimension :: ipiv(nr)
    double precision, allocatable, dimension(:) :: work

    x = 0.d0; Acopy = A 
    
    if (nr.le.nc) then 

        allocate(xsol(nc))
        xsol = 0.d0
        do i = 1, nr
            xsol(i) = b(i)
        end do

        if (nr.eq.nc) then

            ! Determmined system: 
            call dgesv(nr, 1, Acopy, nr, ipiv, xsol, nr, info)

        elseif (nr.lt.nc) then

            ! Under-determined system: 
            allocate(work(2*nr))
            call dgels('N', nr, nc, 1, Acopy, nr, xsol, nc, work, size(work), info)
                
        end if

        x = xsol

    else if (nr.gt.nc) then 
        
        allocate(xsol(nr))
        xsol = b
        
        ! Over-determined system: 
        allocate(work(2*nc))
        call dgels('N', nr, nc, 1, Acopy, nr, xsol, nr, work, size(work), info)

        do i = 1, nc
            x(i) = xsol(i)
        end do

    end if

end subroutine solve_linear_system

subroutine inverse_matrix(nr, A)
    !! Computes the inverse of a square matrix using LU decomposition
    
    implicit none
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nr
    double precision, intent(inout) :: A
    dimension :: A(nr, nr)

    ! Local data
    ! ----------
    integer :: ipiv, info
    double precision :: work, Acopy
    dimension :: work(nr), ipiv(nr), Acopy(nr, nr)

    Acopy = A

    ! Compute LU decomposition
    call dgetrf(nr, nr, Acopy, nr, ipiv, info)

    ! Compute inverse from LU decomposition
    call dgetri(nr, Acopy, nr, ipiv, work, nr, info)

    A = Acopy

end subroutine inverse_matrix

subroutine crossproduct(v1, v2, v3)
    !! Computes cross product in a 3D Euclidean space, v3 =  v1 x v2 (x is cross product)

    ! Input / output data
    ! -------------------
    double precision, intent(in) :: v1, v2
    dimension :: v1(3), v2(3)

    double precision, intent(out) :: v3
    dimension :: v3(3)

    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

end subroutine crossproduct

subroutine kronvec2d(nnz_A, A, nnz_B, B, R, alpha)
    !! Computes kron product of 2 vectors and returns R = R + alpha*(A x B) (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz_A, nnz_B
    double precision, intent(in) :: A, B, alpha
    dimension :: A(nnz_A), B(nnz_B)

    double precision, intent(inout) :: R
    dimension :: R(nnz_A*nnz_B)

    ! Local data
    ! ----------
    integer :: iA, iB, genPos, nb_tasks

    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, size(R)/nb_tasks)
    do iA = 1, nnz_A
        do iB = 1, nnz_B
            genPos = iB + (iA-1)*nnz_B
            R(genPos) = R(genPos) + alpha*A(iA)*B(iB)
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kronvec2d 

subroutine kronvec3d(nnz_A, A, nnz_B, B, nnz_C, C, R, alpha)
    !! Computes the kron product of 3 vectors and returns R = R + alpha*(Au x Av x Aw) (x : tensor product)

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz_A, nnz_B, nnz_C
    double precision, intent(in) :: A, B, C, alpha
    dimension :: A(nnz_A), B(nnz_B), C(nnz_C)

    double precision, intent(inout) :: R
    dimension :: R(nnz_A*nnz_B*nnz_C)

    ! Local data
    ! ----------
    integer :: iA, iB, iC, genPos, nb_tasks

    !$OMP PARALLEL PRIVATE(genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(R)/nb_tasks)
    do iA = 1, nnz_A
        do iB = 1, nnz_B
            do iC = 1, nnz_C
                genPos = iC + (iB-1)*nnz_C + (iA-1)*nnz_C*nnz_B
                R(genPos) = R(genPos) + alpha*A(iA)*B(iB)*C(iC)
            end do
        end do 
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

end subroutine kronvec3d

subroutine spmat_dot_dvec(nr, nc, nnz, indi, indj, A, array_in, array_out)
    !! Computes the dot product of sparse matrix with dense vector. It returns a dense vector
    !! Sparse matrix in CSR format

    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr, nc, nnz 
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: A, array_in
    dimension :: A(nnz), array_in(nc)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr)

    ! Local data
    ! ----------
    integer :: i , j
    double precision :: sum

    array_out = 0.d0

    do i = 1, nr
        sum = 0.d0
        do j = indi(i), indi(i+1)-1
            sum = sum + A(j)*array_in(indj(j))
        end do
        array_out(i) = sum
    end do

end subroutine spmat_dot_dvec

subroutine trapezoidal_rule_2d(nru, nrv, tensor, result)
    !! Computes an integral inside a unitary square using trapezoidal rule.
    !! It supposes that points are equidistant

    implicit none
    ! Input / output data
    ! --------------------
    integer :: nru, nrv
    double precision, intent(in) :: tensor
    dimension :: tensor(nru*nrv)

    double precision, intent(out) :: result
    
    ! Local data
    ! ----------
    integer :: iu, iv
    integer :: indu, indv
    dimension :: indu(2), indv(2)
    double precision :: ttensor
    dimension :: ttensor(nru, nrv)

    result = 0.d0
    indu = (/1, nru/); indv = (/1, nrv/)

    ttensor = reshape(tensor, (/nru, nrv/))

    ! Get internal points
    do iv = 2, nrv-1
        do iu = 2, nru-1
            result = result + ttensor(iu, iv)
        end do
    end do

    ! Get boundary edges
    do iv = 1, 2
        do iu = 2, nru-1
            result = result + ttensor(iu, indv(iv))/2.d0
        end do
    end do

    do iv = 2, nrv-1
        do iu = 1, 2
            result = result + ttensor(indu(iu), iv)/2.d0
        end do
    end do

    ! Get corner points
    do iv = 1, 2
        do iu = 1, 2
            result = result + ttensor(indu(iu), indv(iv))/4.d0
        end do
    end do

    ! Update integral
    result = result/((nru - 1)*(nrv - 1))

end subroutine trapezoidal_rule_2d

subroutine trapezoidal_rule_3d(nru, nrv, nrw, tensor, result)
    !! Computes an integral inside a unitary cube using trapezoidal rule.
    !! It supposes that points are equidistant

    implicit none
    ! Input / output data
    ! --------------------
    integer :: nru, nrv, nrw
    double precision, intent(in) :: tensor
    dimension :: tensor(nru*nrv*nrw)

    double precision, intent(out) :: result
    
    ! Local data
    ! ----------
    integer :: iu, iv, iw
    integer :: indu, indv, indw
    dimension :: indu(2), indv(2), indw(2)
    double precision :: ttensor
    dimension :: ttensor(nru, nrv, nrw)

    result = 0.d0
    indu = (/1, nru/); indv = (/1, nrv/); indw = (/1, nrw/)

    ttensor = reshape(tensor, (/nru, nrv, nrw/))

    ! Get internal points
    do iw = 2, nrw-1
        do iv = 2, nrv-1
            do iu = 2, nru-1
                result = result + ttensor(iu, iv, iw)
            end do
        end do
    end do

    ! Get bounding surface internal points
    do iw = 2, nrw-1
        do iv = 2, nrv-1
            do iu = 1, 2
                result = result + ttensor(indu(iu), iv, iw)/2.d0
            end do
        end do
    end do

    do iw = 2, nrw-1
        do iv = 1, 2
            do iu = 2, nru-1
                result = result + ttensor(iu, indv(iv), iw)/2.d0
            end do
        end do
    end do

    do iw = 1, 2
        do iv = 2, nrv-1
            do iu = 2, nru-1
                result = result + ttensor(iu, iv, indw(iw))/2.d0
            end do
        end do
    end do

    ! Get boundary edges
    do iw = 1, 2
        do iv = 1, 2
            do iu = 2, nru-1
                result = result + ttensor(iu, indv(iv), indw(iw))/4.d0
            end do
        end do
    end do

    do iw = 1, 2
        do iv = 2, nrv-1
            do iu = 1, 2
                result = result + ttensor(indu(iu), iv, indw(iw))/4.d0
            end do
        end do
    end do

    do iw = 2, nrw-1
        do iv = 1, 2
            do iu = 1, 2
                result = result + ttensor(indu(iu), indv(iv), iw)/4.d0
            end do
        end do
    end do

    ! Get corner points
    do iw = 1, 2
        do iv = 1, 2
            do iu = 1, 2
                result = result + ttensor(indu(iu), indv(iv), indw(iw))/8.d0
            end do
        end do
    end do

    ! Update integral
    result = result/((nru - 1)*(nrv - 1)*(nrw - 1))

end subroutine trapezoidal_rule_3d

subroutine symtensor2array(dimen, nvoigt, matrix, array)
    !! Returns the upper triangular part of a matrix

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nvoigt
    double precision, intent(in) :: matrix
    dimension :: matrix(dimen, dimen)
    double precision, intent(out) :: array
    dimension :: array(nvoigt)

    ! Local data
    ! ----------
    integer :: i, j, k

    array = 0.d0; k = 0
    do i = 1, dimen
        k = k + 1
        array(k) = matrix(i, i)
    end do
    
    do i = 1, dimen-1
        do j = i+1, dimen
            k = k + 1
            array(k) = matrix(i, j)
        end do
    end do

end subroutine symtensor2array

subroutine array2symtensor(dimen, nvoigt, array, matrix)
    !! Returns the matrix built from the upper triangular part

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nvoigt
    double precision, intent(in) :: array
    dimension :: array(nvoigt)
    double precision, intent(out) :: matrix
    dimension :: matrix(dimen, dimen)
    
    ! Local data
    ! ----------
    integer :: i, j, k

    matrix = 0.d0; k = 0
    do i = 1, dimen
        k = k + 1
        matrix(i, i) = array(k)
    end do
    
    do i = 1, dimen-1
        do j = i+1, dimen
            k = k + 1
            matrix(i, j) = array(k)
            matrix(j, i) = array(k)
        end do
    end do

end subroutine array2symtensor

! --------
! Indices
! --------

subroutine coo2csr(nm, nr, nnz, a_coo, indi_coo, indj_coo, a_csr, indj_csr, indi_csr)
    !! Change COO format to CSR format

    implicit none 
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nm, nr, nnz
    double precision, intent(in) :: a_coo
    dimension :: a_coo(nnz, nm)
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)

    double precision, intent(out) :: a_csr
    dimension :: a_csr(nnz, nm)
    integer, intent(out) :: indj_csr, indi_csr
    dimension :: indj_csr(nnz), indi_csr(nr+1)

    ! Local data
    ! ----------
    integer :: i, j, k, c, iad
    double precision :: x
    dimension :: x(nm)
    
    ! Get number of nonzero values for each row
    indi_csr = 0
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

    ! Sort COO nonzero values and update values of CSR j-indices
    do k = 1, nnz
        i = indi_coo(k)
        j = indj_coo(k)
        x = a_coo(k, :)
        iad = indi_csr(i)
        a_csr(iad, :) =  x
        indj_csr(iad) = j
        indi_csr(i) = iad + 1
    end do

    ! Update i-indices
    do  j = nr, 1, -1
        indi_csr(j+1) = indi_csr(j)
    end do
    indi_csr(1) = 1

end subroutine coo2csr

subroutine csr2csc(nm, nr, nc, nnz, a_csr, indj_csr, indi_csr, a_csc, indj_csc, indi_csc)
    !! Gets CSC format from CSR format. 
    !! It only works when csr is close to coo format, i.e., nonzero values are sorted in order of occurrence

    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nm, nnz, nr, nc
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_csr
    dimension :: a_csr(nnz, nm)

    integer, intent(out) :: indi_csc, indj_csc
    dimension :: indi_csc(nc+1), indj_csc(nnz)
    double precision, intent(out) :: a_csc
    dimension :: a_csc(nnz, nm)

    ! Local data
    ! ----------
    integer :: indi_coo 
    dimension :: indi_coo(nnz)
    integer :: i, j, c

    c = 0
    do i = 1, nr
        do j = indi_csr(i), indi_csr(i+1) - 1
            c = c + 1
            indi_coo(c) = i
        end do
    end do

    call coo2csr(nm, nc, nnz, a_csr, indj_csr, indi_coo, a_csc, indj_csc, indi_csc)

end subroutine csr2csc

subroutine coo2dense(nnz, indi_coo, indj_coo, a_coo, nr, nc, AA)
    !! Gets a dense matrix from COO format
    !! Repeated indices (i, j) are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_coo, indj_coo
    dimension :: indi_coo(nnz), indj_coo(nnz)
    double precision, intent(in) :: a_coo
    dimension :: a_coo(nnz)

    double precision, intent(out) :: AA
    dimension :: AA(nr, nc)

    ! Local data 
    ! ----------
    integer :: i, j, k

    AA = 0.d0

    do k = 1, nnz
        i = indi_coo(k)
        j = indj_coo(k)
        AA(i, j) = AA(i, j) + a_coo(k)
    end do

end subroutine coo2dense

subroutine csr2dense(nnz, indi_csr, indj_csr, a_csr, nr, nc, AA)
    !! Gets a dense matrix from a CSR format
    !! Repeated indices (i, j) are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nnz, nr, nc 
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_csr
    dimension :: a_csr(nnz)

    double precision, intent(out) :: AA
    dimension :: AA(nr, nc)

    ! Local data
    ! -------------
    integer :: i, j, k

    AA = 0.d0

    do i = 1, nr
        do k = indi_csr(i), indi_csr(i+1)-1
            j = indj_csr(k)
            AA(i, j) = AA(i, j) + a_csr(k)
        end do
    end do
    
end subroutine csr2dense

subroutine dense2csr(nr, nc, AA, nnz, indi_csr, indj_csr, a_csr)
    !! Gets CSR format from matrix 
    !! If nnz < 0, the it only computes the number of nonzero values

    implicit none 
    ! Input / output data
    ! --------------------
    double precision, parameter :: threshold = 1.d-15
    integer, intent(in) :: nr, nc
    double precision, intent(in) :: AA 
    dimension :: AA(nr, nc)

    integer, intent(inout) :: nnz
    integer, intent(out) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(*)
    double precision, intent(out) :: a_csr
    dimension :: a_csr(*)

    ! Local data
    ! ----------
    integer :: i, j, k, l

    if (nnz.le.0) then 
        ! Gets only the number of nonzero values of the matrix 
        nnz = 0
        do i = 1, nr
            do j = 1, nc
                if (abs(AA(i, j)).gt.threshold) then
                    nnz = nnz + 1
                end if
            end do
        end do

    else
        ! Get CSR format from dense matrix 
        k = 1; indi_csr(1) = 1
        do i = 1, nr
            l = 0
            do j = 1, nc
                if (abs(AA(i, j)).gt.threshold) then
                    a_csr(k) = AA(i, j)
                    indj_csr(k) = j
                    k = k + 1; l = l + 1
                end if
            end do
            indi_csr(i+1) = indi_csr(i) + l
        end do

    end if

end subroutine dense2csr

subroutine get_indices_kron2_product(nr_A, nc_A, nnz_A, indi_A, indj_A, &
                                    nr_B, nc_B, nnz_B, indi_B, indj_B, &  
                                    nnz_C, indi_C, indj_C)
    !! Gets indices of the kronecker product A x B = C 
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

    integer, intent(out) :: indi_C, indj_C
    dimension :: indi_C(nr_A*nr_B+1), indj_C(nnz_C)

    ! Loca data
    ! ---------
    integer :: iA, iB, genPos, jA, jB, c, nb_tasks
    integer :: dummy, nnz_row_A, nnz_row_B, nnz_row_C
    integer, allocatable, dimension(:) :: indj_temp

    dummy = nc_A

    ! Set indices i in csr format
    indi_C(1) = 1
    do iA = 1, nr_A
        do iB = 1, nr_B          
            ! Set number of nonzero elements 
            nnz_row_A = indi_A(iA+1) - indi_A(iA) 
            nnz_row_B = indi_B(iB+1) - indi_B(iB)
            nnz_row_C = nnz_row_A * nnz_row_B

            ! Update value 
            genPos = iB + (iA - 1)*nr_B
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
    
end subroutine get_indices_kron2_product

subroutine get_indices_kron3_product(nr_A, nc_A, nnz_A, indi_A, indj_A, &
                                    nr_B, nc_B, nnz_B, indi_B, indj_B, &  
                                    nr_C, nc_C, nnz_C, indi_C, indj_C, &
                                    nnz_D, indi_D, indj_D)
    !! Returns indices of the kronecker product A x B x C = D
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

    integer, intent(out) :: indi_D, indj_D
    dimension :: indi_D(nr_A*nr_B*nr_C+1), indj_D(nnz_D)

    ! Loca data
    ! ---------
    integer :: iA, iB, iC, genPos, jA, jB, jC, c, nb_tasks
    integer :: dummy, nnz_row_A, nnz_row_B, nnz_row_C, nnz_row_D
    integer, allocatable, dimension(:) :: indj_temp
    
    dummy = nc_A

    ! Set indices i in CSR format
    indi_D(1) = 1
    do iA = 1, nr_A
        do iB = 1, nr_B
            do iC = 1, nr_C
                ! Set number of nonzero elements   
                nnz_row_A = indi_A(iA+1) - indi_A(iA) 
                nnz_row_B = indi_B(iB+1) - indi_B(iB)
                nnz_row_C = indi_C(iC+1) - indi_C(iC)
                nnz_row_D = nnz_row_A * nnz_row_B * nnz_row_C

                ! Update value 
                genPos = iC + (iB-1)*nr_C + (iA - 1)*nr_C*nr_B
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
    
end subroutine get_indices_kron3_product
