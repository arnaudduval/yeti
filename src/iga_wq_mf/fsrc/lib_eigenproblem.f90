!! This file contains basic algorithms to compute the maximum eigenvalue of 
!! the generalized eigenproblem A x = lambda B x
!! We consider dense matrices but it is possible to adapt them with a matrix-free approach
!! More information on "Lecture Notes on Solving Large Scale Eigenvalue Problems" by Arbenz, P. 
!! Author: Joaquin Cornejo

subroutine compute_eigs(nr, A, B, rho, x)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr
    double precision, intent(in)  :: A, B
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

end subroutine compute_eigs

subroutine power_iteration(nr, A, nbIter, x, rho)

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
    double precision :: prod

    do i = 1, d
        prod = dot_product(VV1(i, :), VV2(i, :))
        MM(i, i) = prod 
    end do

    do i = 1, d
        do j = i+1, d
            prod = dot_product(VV1(i, :), VV2(j, :))
            MM(i, j) = prod; MM(j, i) = prod
        end do 
    end do

end subroutine rayleigh_submatrix

subroutine rayleighquotient_max(nr, A, B, C, x, rho, nbIter, threshold)
    !! Preconditioned conjugate gradient - Rayleigh quotient maximization for the computation 
    !! of the greatest eigenvalue of A x = lambda B x
    !! C is a preconditioner

    implicit none
    ! Input / output data
    ! ------------------- 
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr, nbIter
    double precision, intent(in) :: A, B, C, threshold
    dimension :: A(nr, nr), B(nr, nr), C(nr, nr)

    double precision, intent(out) :: x, rho
    dimension :: x(nr)
    
    ! Local data
    ! ----------
    integer :: k, ii
    double precision, dimension(dimen, nr) :: RM1, RM2, RM3
    double precision, dimension(dimen, dimen) :: AA1, BB1, qq
    double precision, dimension(dimen) :: ll, delta
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
        call rayleigh_submatrix(dimen, nr, RM1, RM2, AA1)
        call rayleigh_submatrix(dimen, nr, RM1, RM3, BB1)

        call compute_eigs(dimen, AA1, BB1, ll, qq)
        rho = maxval(ll); ii = maxloc(ll, dim=1)
        delta = qq(:, ii)

        x = x*delta(1) + p*delta(2)
        u = matmul(B, x)
        q = sqrt(dot_product(x, u))
        x = x/q; u = u/q
        v = matmul(A, x)
        norm = norm2(g)
    end do

end subroutine rayleighquotient_max

subroutine locally_optimal_block_pcg(nr, A, B, C, x, rho, nbIter, threshold)
    !! Locally optimal block preconditioned conjugate gradient algorithm (LOBCPG) for 
    !! computing the greatest eigenvalue of A x = lambda B x
    !! C is a preconditioner

    implicit none
    ! Input / output data
    ! ------------------- 
    integer, parameter :: dimen = 3
    double precision, parameter :: tol_singular = 1.d-10
    integer, intent(in) :: nr, nbIter
    double precision, intent(in) :: A, B, C, threshold
    dimension :: A(nr, nr), B(nr, nr), C(nr, nr)

    double precision, intent(out) :: x, rho
    dimension :: x(nr)

    ! Local data
    ! ----------
    integer :: k, ii
    double precision, dimension(dimen, nr) :: RM1, RM2, RM3
    double precision, dimension(dimen, dimen) :: AA1, BB1, qq
    double precision, dimension(dimen) :: ll, delta
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

        call rayleigh_submatrix(dimen, nr, RM1, RM2, AA1)
        call rayleigh_submatrix(dimen, nr, RM1, RM3, BB1)

        if (norm2(p).lt.tol_singular) then
            qq = 0.d0; ll = 0.d0
            call compute_eigs(dimen-1, AA1(:2, :2), BB1(:2, :2), ll(:2), qq(:2, :2))
        else
            call compute_eigs(dimen, AA1, BB1, ll, qq)
        end if

        rho = maxval(ll); ii = maxloc(ll, dim=1)
        delta = qq(:, ii)

        p = -g*delta(2) + p*delta(3)
        x = x*delta(1) + p
        u = matmul(B, x)
        q = sqrt(dot_product(x, u))
        x = x/q; u = u/q
        v = matmul(A, x)
        norm = norm2(g)
    end do

end subroutine locally_optimal_block_pcg
