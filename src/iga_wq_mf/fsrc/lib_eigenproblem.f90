! ==============================================
!! This file contains basic algorithms to compute the maximum eigenvalue of 
!! the generalized eigenproblem A x = lambda B x
!! We consider dense matrices but it is possible to adapt them with a matrix-free approach
!! More information on "Lecture Notes on Solving Large Scale Eigenvalue Problems" by Arbenz, P. 
!! Author: Joaquin Cornejo
! ==============================================

subroutine compute_eigdecomp_pdr(nr, A, B, rho, x)
    !! Computes all the eigenvalues of the generalized eigenproblem A x = rho B x
    !! using dsygvd from LAPACK libraries. 
    !! Matrix A is must be positive definite and B could be semi-positive definite real matrices
    !! Thus A = U^T U and B = U^T D U, where U and D are the eigenvector and eigenvalues respectively

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double precision, intent(in) :: A, B
    dimension :: A(nr, nr), B(nr, nr)

    double precision, intent(out) :: rho, x
    dimension :: rho(nr), x(nr, nr)

    ! Local data
    ! ----------
    integer :: info, lwork, liwork, idum(1)
    double precision, allocatable, dimension(:) :: work, iwork
    double precision :: dummy(1)

    ! Get optimal size
    call dsygvd(1, 'V', 'L', nr, A, nr, B, nr, rho, dummy, -1, idum, -1, info)
    lwork  = max(1+(6+2*nr)*nr, nint(dummy(1)))
    liwork = max(3+5*nr, idum(1))
    allocate (work(lwork), iwork(liwork))

    ! Get eigen decomposition
    x = A
    call dsygvd(1, 'V', 'L', nr, x, nr, B, nr, rho, work, lwork, iwork, liwork, info)
    deallocate(work, iwork)

end subroutine compute_eigdecomp_pdr

subroutine compute_schurdecomp_gc(nr, A, B, U, V, S, T)
    !! Computes the Schur decomposition of the pencil (A, B) such that
    !! A = U S V^H and B = U T V^H. H means the hermitian transpose (with conjugates)
    !! S and T are triangular matrices
    !! using dsygvd from LAPACK libraries. 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double precision, intent(in) :: A, B
    dimension :: A(nr, nr), B(nr, nr)

    double complex, intent(out) :: U, V, S, T
    dimension :: U(nr, nr), V(nr, nr), S(nr, nr), T(nr, nr)

    ! Local data
    ! ----------
    double complex :: alpha(nr), beta(nr)
    logical :: bwork(nr), selctg
    double precision, allocatable, dimension(:) :: rwork
    double complex, allocatable, dimension(:) :: work
    integer :: lwork, sdim, info

    S = dcmplx(A); T = dcmplx(B)
    allocate(work(2*nr), rwork(8*nr))
    lwork = size(work)
    call zgges('V', 'V', 'N', selctg, nr, S, nr, T, nr, sdim, alpha, beta, &
                U, nr, V, nr, work, lwork, rwork, bwork, info)
    
end subroutine compute_schurdecomp_gc

subroutine power_iteration(nr, A, x, rho, nbIter)
    !! Computes the maximum eigenvalue of A and its eigenvector. 
    !! The matrix A could be any matrix

    implicit none
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nr, nbIter
    double precision, intent(in) :: A
    dimension :: A(nr, nr)

    double precision, intent(out) :: x, rho
    dimension :: x(nr)

    ! Local data
    ! ----------
    integer :: i
    double precision :: norm, x1
    dimension :: x1(nr)

    call random_number(x)
    norm = norm2(x)
    x    = x/norm

    do i = 1, nbIter
        x1   = matmul(A, x)
        norm = norm2(x1)
        x    = x1/norm
    end do

    rho  = dot_product(x, matmul(A, x))
    
end subroutine power_iteration

subroutine rayleigh_submatrix(d, nr, VV1, VV2, MM)
    !! Computes the Rayleigh submatrix from vectors VV1 and VV2
    !! MM_ij = dot_product(VV1_i, VV2_j)

    implicit none
    ! Input / output data
    ! ------------------- 
    integer :: d, nr
    double precision, intent(in) :: VV1, VV2
    double precision, intent(out) :: MM
    dimension :: VV1(d, nr), VV2(d, nr), MM(d, d)

    ! Local data
    ! ----------
    integer :: i, j

    do i = 1, d
        MM(i, i) = dot_product(VV1(i, :), VV2(i, :)) 
    end do

    do i = 1, d
        do j = i+1, d
            MM(i, j) = dot_product(VV1(i, :), VV2(j, :))
            MM(j, i) = MM(i, j)
        end do 
    end do

end subroutine rayleigh_submatrix

subroutine rayleighquotient(nr, A, B, C, x, rho, nbIter, threshold)
    !! Preconditioned conjugate gradient - Rayleigh quotient maximization for the computation 
    !! of the greatest eigenvalue of A x = rho B x
    !! C is a preconditioner
    !! Matrix A must be positive definite (PD) and B could be semi-positive definite
    !! Some papers shows that even if A is not symetric (and then not necessarily PD) the algorithm converges

    implicit none
    ! Input / output data
    ! ------------------- 
    integer, parameter :: d = 2
    integer, intent(in) :: nr, nbIter
    double precision, intent(in) :: A, B, C, threshold
    dimension :: A(nr, nr), B(nr, nr), C(nr, nr)

    double precision, intent(out) :: x, rho
    dimension :: x(nr)
    
    ! Local data
    ! ----------
    integer :: k, ii
    double precision, dimension(d, nr) :: RM1, RM2, RM3
    double precision, dimension(d, d) :: AA1, BB1, qq
    double precision, dimension(d) :: ll, delta
    double precision, dimension(nr) :: u, v, g, gtil, p
    double precision :: q, norm

    call random_number(x)
    norm = norm2(x)
    x = x/norm

    u = matmul(B, x)
    q = sqrt(dot_product(x, u))
    x = x/q; u = u/q
    v = matmul(A, x)
    rho = dot_product(x, v)

    g = x; norm = 1.d0

    do k = 1, nbIter
        if (norm.le.threshold) return

        gtil = g
        call solve_linear_system(nr, nr, C, 2.d0*(v - rho*u), g)
        if (k.eq.1) p = -g
        if (k.gt.1) p = -g + dot_product(g, matmul(B, g))/dot_product(gtil, matmul(B, gtil))*p

        RM1(1, :) = x; RM1(2, :) = p; RM2(1, :) = v; RM3(1, :) = u
        RM2(2, :) = matmul(A, p); RM3(2, :) = matmul(B, p)
        call rayleigh_submatrix(d, nr, RM1, RM2, AA1)
        call rayleigh_submatrix(d, nr, RM1, RM3, BB1)

        call compute_eigdecomp_pdr(d, AA1, BB1, ll, qq)
        rho = maxval(ll); ii = maxloc(ll, dim=1)
        delta = qq(:, ii)

        x = x*delta(1) + p*delta(2)
        u = matmul(B, x)
        q = sqrt(dot_product(x, u))
        x = x/q; u = u/q
        v = matmul(A, x)
        norm = norm2(g)
    end do

end subroutine rayleighquotient

subroutine locally_optimal_block_pcg(nr, A, B, C, x, rho, ishigher, nbIter, threshold)
    !! Locally optimal block preconditioned conjugate gradient algorithm (LOBCPG) for 
    !! computing the greatest eigenvalue of A x = rho B x
    !! C is a preconditioner
    !! Matrix A must be positive definite (PD) and B could be semi-positive definite
    !! Some papers shows that even if A is not symetric (and then not necessarily PD) 
    !! the algorithm still could converges

    implicit none
    ! Input / output data
    ! ------------------- 
    integer, parameter :: d = 3
    integer, intent(in) :: nr, nbIter
    logical, intent(in) :: ishigher
    double precision, intent(in) :: A, B, C, threshold
    dimension :: A(nr, nr), B(nr, nr), C(nr, nr)

    double precision, intent(out) :: x, rho
    dimension :: x(nr)

    ! Local data
    ! ----------
    integer :: k, ii
    double precision, dimension(d, nr) :: RM1, RM2, RM3
    double precision, dimension(d, d) :: AA1, BB1
    double precision, dimension(d) :: delta
    double precision, dimension(nr) :: u, v, g, gtil, p
    double precision :: q, norm
    double precision, allocatable, dimension(:) ::  ll
    double precision, allocatable, dimension(:, :) ::  qq

    call random_number(x)
    norm = norm2(x)
    x = x/norm

    u = matmul(B, x)
    q = sqrt(dot_product(x, u))
    x = x/q; u = u/q
    v = matmul(A, x)
    rho = dot_product(x, v)
    p = 0.d0
    norm = 1.d0

    do k = 1, nbIter
        if (norm.le.threshold) return

        g = v - rho*u
        norm = norm2(g)
        call solve_linear_system(nr, nr, C, g, gtil)
        g = gtil

        RM1(1, :) = x; RM1(2, :) = -g; RM1(3, :) = p
        RM2(1, :) = v; RM2(2, :) = -matmul(A, g); RM2(3, :) = matmul(A, p)
        RM3(1, :) = u; RM3(2, :) = -matmul(B, g); RM3(3, :) = matmul(B, p)

        call rayleigh_submatrix(d, nr, RM1, RM2, AA1); AA1 = 0.5d0*(AA1 + transpose(AA1))
        call rayleigh_submatrix(d, nr, RM1, RM3, BB1); BB1 = 0.5d0*(BB1 + transpose(BB1))

        if (k.eq.1) then
            allocate(ll(d-1), qq(d-1, d-1))
            call compute_eigdecomp_pdr(size(ll), AA1(:d-1, :d-1), BB1(:d-1, :d-1), ll, qq)
        else
            allocate(ll(d), qq(d, d))
            call compute_eigdecomp_pdr(size(ll), AA1, BB1, ll, qq)
        end if

        if (ishigher) then
            rho = maxval(ll); ii = maxloc(ll, dim=1)
        else
            rho = minval(ll); ii = minloc(ll, dim=1)
        end if

        delta = 0.d0
        if (k.eq.1) then
            delta(:2) = qq(:, ii)
        else
            delta = qq(:, ii)
        end if

        p = -g*delta(2) + p*delta(3)
        x = x*delta(1) + p
        u = matmul(B, x)
        q = sqrt(dot_product(x, u))
        x = x/q; u = u/q
        v = matmul(A, x)
        norm = norm2(g)
        deallocate(ll, qq)
    end do

end subroutine locally_optimal_block_pcg

subroutine GMRES(nr, A, b, P, x, nbIter, nbRestarts, threshold)
    implicit none
    ! Input / output data
    ! ------------------- 
    double precision :: threshold
    integer, intent(in) :: nr, nbIter, nbRestarts
    double precision, intent(in) :: A, b, P
    dimension :: A(nr, nr), b(nr), P(nr, nr)
    double precision, intent(out) :: x
    dimension :: x(nr)

    ! Local data
    ! ----------
    double precision :: H, V, Z, beta, e1, y, rho
    dimension :: H(nbIter+1, nbIter), V(nbIter+1, nr), Z(nbIter+1, nr), beta(nbRestarts), e1(nbIter+1), y(nbIter)
    double precision :: r(nr), w(nr)
    integer :: i, j, k

    e1 = 0.d0; y = 0.d0; H = 0.d0; V = 0.d0; Z = 0.d0; beta = 0.d0; x = 0.d0
    do k = 1, nbRestarts
        r = b - matmul(A, x)
        beta(k) = norm2(r)
        if (beta(k).le.threshold*beta(1)) exit
        V(1, :) = r/beta(1)
        e1(1) = beta(k)

        do j = 1, nbIter
            Z(j, :) = matmul(P, V(j, :))
            w = matmul(A, Z(j, :))
            do i = 1, j
                H(i, j) = dot_product(w, V(i, :))
                w = w - H(i, j)*V(i, :)
            end do
            H(j+1, j) = norm2(w)
            if (abs(H(j+1, j)).gt.1e-10) then
                V(j+1, :) = w/H(j+1, j)
            end if
            call solve_linear_system(j+1, j, H(:j+1, :j), e1(:j+1), y(:j))
            rho = norm2(matmul(H(:j+1, :j), y(:j)) - e1(:j+1))
            if (rho.le.threshold*beta(1)) exit
        end do
        x = x + matmul(y(:j), Z(:j, :))
    end do

end subroutine GMRES