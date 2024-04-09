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
    double precision :: tmp(nnz)
    integer :: i, j

    if ((repeat.ge.nnz).or.(repeat.lt.0)) stop 'Something went wrong in diff array'

    array_out = array_in

    do j = 1, repeat

        tmp = array_out

        ! Get the difference
        do i = 1, nnz-j
            array_out(i) = tmp(i+1) - tmp(i)
        end do

        ! Extra data is set to 0
        do i = nnz-j+1, nnz
            array_out(i) = 0.d0
        end do

    end do
    
end subroutine diff_array

subroutine find_unique_array(nnz_in, array_in, array_out)
    !! Gets the non-repreated values of an array
    !! Ex: Given the array [0, 1, 1, 2, 3, 3, 3], the unique array is [0, 1, 2, 3, 0, 0, 0, 4]
    !! The last value is the number of unique elements

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz_in 
    double precision, intent(in) :: array_in
    dimension :: array_in(nnz_in)

    double precision, intent(out) :: array_out
    dimension :: array_out(nnz_in+1)

    ! Local data
    ! ----------
    integer :: i, c, nnz_out
    logical :: mask(nnz_in)

    array_out = 0.d0
    mask = .false.

    do i = 1, nnz_in

        ! Count the number of occurrences of this element
        c = count(array_in(i).eq.array_in)
    
        if (c.eq.1) then  
            ! There is only one, flag it:  
            mask(i) = .true.  
        else  
            !  Flag this value only if it has not already been flagged
            if (.not. any(array_in(i).eq.array_in .and. mask)) then
                mask(i) = .true.  
            end if
        end if
    
    end do

    ! Return only flagged elements
    nnz_out = count(mask)
    array_out(1:nnz_out) = pack(array_in, mask)
    array_out(nnz_in+1) = dble(nnz_out)

end subroutine find_unique_array

subroutine indices2list(size_nclist, nbPts, ptsPos, nclist, sample)
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: size_nclist, nbPts
    integer, intent(in) :: ptsPos(size_nclist, nbPts), nclist(size_nclist)
    integer, intent(out) :: sample(nbPts**size_nclist)

    ! Local data
    ! ----------
    integer :: i, j, k, l, c, gp
    
    c = 1
    if (size_nclist.eq.2) then
        do j = 1, nbPts
            do i = 1, nbPts
                gp = ptsPos(1, i) + (ptsPos(2, j) - 1)*nclist(1)
                sample(c) = gp; c = c + 1
            end do
        end do
        
    else if (size_nclist.eq.3) then
        do k = 1, nbPts
            do j = 1, nbPts
                do i = 1, nbPts
                    gp = ptsPos(1, i) + (ptsPos(2, j) - 1)*nclist(1) + (ptsPos(3, k) - 1)*nclist(1)*nclist(2)
                    sample(c) = gp; c = c + 1
                end do
            end do
        end do

    else if (size_nclist.eq.4) then
        do l = 1, nbPts
            do k = 1, nbPts
                do j = 1, nbPts
                    do i = 1, nbPts
                        gp = ptsPos(1, i) + (ptsPos(2, j) - 1)*nclist(1) &
                            + (ptsPos(3, k) - 1)*nclist(1)*nclist(2) &
                            + (ptsPos(4, l) - 1)*nclist(1)*nclist(2)*nclist(3)
                        sample(c) = gp; c = c + 1
                    end do
                end do
            end do
        end do
    end if
end subroutine

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

    if (nnz.le.1) stop 'Try greater than one'

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

        if (ncA.ne.ncB) stop 'Size problem in GEMM'
        allocate(MM(nrA, ncA))
        forall (i = 1:nrA, j = 1:ncA) 
            MM(i, j) = A(i, j) * W(j)
        end forall
        R = matmul(MM, transpose(B))

    else if (type.eq.2) then 
        ! A diag(W) B.T
        ! nrR = nrA, ncR = nrB, size(W) = ncA = ncB
        ! Only diagonal

        if (ncA.ne.ncB) stop 'Size problem in GEMM'
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

        if (nrA.ne.nrB) stop 'Size problem in GEMM'
        allocate(MM(nrB, ncB))
        forall (i = 1:nrB, j = 1:ncB) 
            MM(i, j) = W(i) * B(i, j) 
        end forall
        R = matmul(transpose(A), MM)

    else if (type.eq.4) then 
        ! A.T diag(W) B
        ! nrR = ncA, ncR = ncB, size(W) = nrA = nrB
        ! Only diagonal

        if (nrA.ne.nrB) stop 'Size problem in GEMM'
        do i = 1, ncA
            s = 0.d0
            do k = 1, nrA
                s = s + A(k, i)*W(k)*B(k, i)
            end do
            R(i, i) = s
        end do

    end if

end subroutine gemm_AWB

subroutine block_dot_product(nm, nr, A, B, result)
    !! Computes dot product of A and B. Both are actually vectors arranged following each dimension
    !! Vector A is composed of [Au, Av, Aw] and B of [Bu, Bv, Bw]. 
    !! Dot product A.B = Au.Bu + Av.Bv + Aw.Bw 

    implicit none
    ! Input/ output data
    ! ------------------
    integer, intent(in) :: nm, nr
    double precision, intent(in) :: A, B
    dimension :: A(nm, nr), B(nm, nr)

    double precision :: result

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp

    result = 0.d0
    do i = 1, nm 
        tmp = dot_product(A(i, :), B(i, :))
        result = result + tmp
    end do

end subroutine block_dot_product

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

subroutine solve_uppertriangular_system(nr, A, b, x)
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double precision, intent(in) :: A, b
    dimension :: A(nr, nr), b(nr)

    double precision, intent(out) :: x
    dimension :: x(nr)

    ! Local data
    ! ----------
    integer :: nhrs, info

    x = b; nhrs = 1
    call dtrtrs('U', 'N', 'N', nr, nhrs, A, nr, x, nr, info)

end subroutine solve_uppertriangular_system

subroutine solve_complex_uppertriangular_system(nr, A, b, x)
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double complex, intent(in) :: A, b
    dimension :: A(nr, nr), b(nr)

    double complex, intent(out) :: x
    dimension :: x(nr)

    ! Local data
    ! ----------
    integer :: nhrs, info

    x = b; nhrs = 1
    call ztrtrs('U', 'N', 'N', nr, nhrs, A, nr, x, nr, info)

end subroutine solve_complex_uppertriangular_system

subroutine solve_banded_system(nr, ML, MU, ABanded, b, x)
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ML, MU
    double precision, intent(in) :: ABanded, b
    dimension :: ABanded(2*ML+MU+1, nr), b(nr)
    double precision, intent(out) :: x(nr)

    ! Local data
    ! ----------
    integer :: nhrs, info
    integer :: ipiv(nr)

    x = b; nhrs = 1
    call dgbsv(nr, ML, MU, nhrs, ABanded, size(ABanded, dim=1), ipiv, x, nr, info)

end subroutine solve_banded_system

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

subroutine convert2bandedstorage(nr, ML, MU, A, B)
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ML, MU
    double precision, intent(in) :: A
    dimension :: A(nr, nr)
    double precision, intent(out) :: B
    dimension :: B(2*ML+MU+1, nr)
    ! Local data
    ! ----------
    integer :: i, j, k

    k = ML + MU + 1
    B = 0.d0
    do j = 1, nr
        do i = max(1, j - MU), min(nr, j + ML)
            B(k+i-j, j) = A(i, j)
        end do
    end do

end subroutine convert2bandedstorage

subroutine crossproduct(vA, vB, vC)
    !! Computes cross product in a 3D Euclidean space, v3 =  v1 x v2 (x is cross product)

    ! Input / output data
    ! -------------------
    double precision, intent(in) :: vA, vB
    dimension :: vA(3), vB(3)

    double precision, intent(out) :: vC
    dimension :: vC(3)

    vC(1) = vA(2)*vB(3) - vA(3)*vB(2)
    vC(2) = vA(3)*vB(1) - vA(1)*vB(3)
    vC(3) = vA(1)*vB(2) - vA(2)*vB(1)

end subroutine crossproduct

subroutine kronvec2d(nnz_A, A, nnz_B, B, R, alpha)
    !! Computes kron product of 2 vectors and returns R = R + alpha*(A x B) (x : tensor product)

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
    integer :: iA, iB, gp

    do iA = 1, nnz_A
        do iB = 1, nnz_B
            gp = iB + (iA-1)*nnz_B
            R(gp) = R(gp) + alpha*A(iA)*B(iB)
        end do 
    end do

end subroutine kronvec2d 

subroutine kronvec3d(nnz_A, A, nnz_B, B, nnz_C, C, R, alpha)
    !! Computes the kron product of 3 vectors and returns R = R + alpha*(Au x Av x Aw) (x : tensor product)

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
    integer :: iA, iB, iC, gp

    do iA = 1, nnz_A
        do iB = 1, nnz_B
            do iC = 1, nnz_C
                gp = iC + (iB-1)*nnz_C + (iA-1)*nnz_C*nnz_B
                R(gp) = R(gp) + alpha*A(iA)*B(iB)*C(iC)
            end do
        end do 
    end do

end subroutine kronvec3d

subroutine kronvec4d(nnz_A, A, nnz_B, B, nnz_C, C, nnz_D, D, R, alpha)
    !! Computes the kron product of 3 vectors and returns R = R + alpha*(Au x Av x Aw) (x : tensor product)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz_A, nnz_B, nnz_C, nnz_D
    double precision, intent(in) :: A, B, C, D, alpha
    dimension :: A(nnz_A), B(nnz_B), C(nnz_C), D(nnz_D)

    double precision, intent(inout) :: R
    dimension :: R(nnz_A, nnz_B, nnz_C, nnz_D)

    ! Local data
    ! ----------
    integer :: iB, iC, iD

    do iD = 1, nnz_D
        do iC = 1, nnz_C
            do iB = 1, nnz_B
                R(:, iB, iC, iD) = R(:, iB, iC, iD) + alpha*A(:)*B(iB)*C(iC)*D(iD)
            end do
        end do 
    end do

end subroutine kronvec4d

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
    integer :: i , k

    array_out = 0.d0

    do i = 1, nr
        do k = indi(i), indi(i+1)-1
            array_out(i) = array_out(i) + A(k)*array_in(indj(k))
        end do
    end do

end subroutine spmat_dot_dvec

subroutine spmat_dot_dmat(nr, nnz, indi, indj, A, nr_in, nc_in, matrix_in, matrix_out)
    !! Computes the dot product of sparse matrix with dense vector. It returns a dense vector
    !! Sparse matrix in CSR format

    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr, nr_in, nc_in, nnz 
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: A, matrix_in
    dimension :: A(nnz), matrix_in(nr_in, nc_in)

    double precision, intent(out) :: matrix_out
    dimension :: matrix_out(nr, nc_in)

    ! Local data
    ! ----------
    integer :: i, indl, indr

    do i = 1, nr
        indl = indi(i)
        indr = indi(i+1) - 1
        matrix_out(i, :) = matmul(A(indl:indr), matrix_in(indj(indl:indr), :))
    end do

end subroutine spmat_dot_dmat

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
    double precision :: w_u(nru), w_v(nrv)
    double precision :: ttensor
    dimension :: ttensor(nru, nrv, 1, 1)

    ttensor = reshape(tensor, (/nru, nrv, 1, 1/))
    w_u = 1.d0; w_u(1) = 0.5d0; w_u(nru) = 0.5d0 
    w_v = 1.d0; w_v(1) = 0.5d0; w_v(nrv) = 0.5d0 

    result = 0.d0
    do iv = 1, nrv
        do iu = 1, nru
            result = result + ttensor(iu, iv, 1, 1)*w_u(iu)*w_v(iv) 
        end do
    end do
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
    double precision :: w_u(nru), w_v(nrv), w_w(nrw)
    double precision :: ttensor
    dimension :: ttensor(nru, nrv, nrw, 1)

    ttensor = reshape(tensor, (/nru, nrv, nrw, 1/))
    w_u = 1.d0; w_u(1) = 0.5d0; w_u(nru) = 0.5d0 
    w_v = 1.d0; w_v(1) = 0.5d0; w_v(nrv) = 0.5d0 
    w_w = 1.d0; w_w(1) = 0.5d0; w_w(nrw) = 0.5d0 

    result = 0.d0
    do iw = 1, nrw
        do iv = 1, nrv
            do iu = 1, nru
                result = result + ttensor(iu, iv, iw, 1)*w_u(iu)*w_v(iv)*w_w(iw)
            end do
        end do
    end do
    result = result/((nru - 1)*(nrv - 1)*(nrw - 1))

end subroutine trapezoidal_rule_3d

subroutine trapezoidal_rule_4d(nru, nrv, nrw, nrt, tensor, result)
    !! Computes an integral inside a unitary cube using trapezoidal rule.
    !! It supposes that points are equidistant

    implicit none
    ! Input / output data
    ! --------------------
    integer :: nru, nrv, nrw, nrt
    double precision, intent(in) :: tensor
    dimension :: tensor(nru*nrv*nrw*nrt)

    double precision, intent(out) :: result
    
    ! Local data
    ! ----------
    integer :: iu, iv, iw, it
    double precision :: w_u(nru), w_v(nrv), w_w(nrw), w_t(nrt)
    double precision :: ttensor
    dimension :: ttensor(nru, nrv, nrw, 1)

    ttensor = reshape(tensor, (/nru, nrv, nrw, nrt/))
    w_u = 1.d0; w_u(1) = 0.5d0; w_u(nru) = 0.5d0 
    w_v = 1.d0; w_v(1) = 0.5d0; w_v(nrv) = 0.5d0 
    w_w = 1.d0; w_w(1) = 0.5d0; w_w(nrw) = 0.5d0 
    w_t = 1.d0; w_t(1) = 0.5d0; w_t(nrt) = 0.5d0 

    result = 0.d0
    do it = 1, nrt
        do iw = 1, nrw
            do iv = 1, nrv
                do iu = 1, nru
                    result = result + ttensor(iu, iv, iw, it)*w_u(iu)*w_v(iv)*w_w(iw)*w_t(it)
                end do
            end do
        end do
    end do
    result = result/((nru - 1)*(nrv - 1)*(nrw - 1)*(nrt - 1))

end subroutine trapezoidal_rule_4d

subroutine get_indexes_kron_product(nc_A, nnz_A, indj_A, &
                                nnz_B, indj_B, indj_C)
    !! Returns indexes of A x B = C (x : kronecker product)
    !! Where A and B are sparse matrices in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! ----------------------
    integer, intent(in) ::  nc_A, nnz_A, nnz_B
    integer, intent(in) :: indj_A, indj_B
    dimension :: indj_A(nnz_A), indj_B(nnz_B)

    integer, intent(out) :: indj_C
    dimension :: indj_C(nnz_A*nnz_B)

    ! Loca data
    ! -----------
    integer :: i, j

    do i = 1, nnz_A
        do j = 1, nnz_B
            indj_C(i+(j-1)*nnz_A) = indj_A(i) + (indj_B(j)-1)*nc_A
        end do
    end do

end subroutine get_indexes_kron_product

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

subroutine erase_rows_csr(nr2er, rows2er, nm, nr_in, nnz_in, indi_in, indj_in, a_in, &
                        nnz_out, indi_out, indj_out, a_out)
    !! Erase rows from CSR format
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr2er, nr_in, nnz_in, nm
    integer, intent(in) :: rows2er, indi_in, indj_in
    dimension :: rows2er(nr2er), indi_in(nr_in+1), indj_in(nnz_in)
    double precision, intent(in) :: a_in
    dimension :: a_in(nnz_in, nm)

    integer, intent(inout) :: nnz_out
    integer, intent(out) :: indi_out, indj_out
    dimension :: indi_out(nr_in-nr2er+1), indj_out(nnz_out)
    double precision, intent(out) :: a_out
    dimension :: a_out(nnz_out, nm) 

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

! --------
! Indices
! --------

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

subroutine multicoo2csr(nm, nr, nnz, a_coo, indi_coo, indj_coo, a_csr, indj_csr, indi_csr)
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

end subroutine multicoo2csr

subroutine multicsr2csc(nm, nr, nc, nnz, a_csr, indj_csr, indi_csr, a_csc, indj_csc, indi_csc)
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

    call multicoo2csr(nm, nc, nnz, a_csr, indj_csr, indi_coo, a_csc, indj_csc, indi_csc)

end subroutine multicsr2csc

subroutine multicsr2dense(nm, nnz, indi_csr, indj_csr, a_csr, nr, nc, AA)
    !! Gets a dense matrix from a CSR format
    !! Repeated indices (i, j) are added

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nm, nnz, nr, nc 
    integer, intent(in) :: indi_csr, indj_csr
    dimension :: indi_csr(nr+1), indj_csr(nnz)
    double precision, intent(in) :: a_csr
    dimension :: a_csr(nnz, nm)

    double precision, intent(out) :: AA
    dimension :: AA(nr, nc, nm)

    ! Local data
    ! -------------
    integer :: i, j, k, l

    AA = 0.d0

    do i = 1, nr
        do k = indi_csr(i), indi_csr(i+1)-1
            j = indj_csr(k)
            do l = 1, nm
                AA(i, j, l) = AA(i, j, l) + a_csr(k, l)
            end do
        end do
    end do
    
end subroutine multicsr2dense
