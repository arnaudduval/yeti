subroutine reset_dirichletbound2(nr, A, ndu, ndv, dod_u, dod_v)
    !! Set to 0 (Dirichlet condition) the values of an array using the dod indices in each dimension
    !! A is actually a vector arranged following each dimension [Au, Av, Aw]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ndu, ndv
    double precision, intent(inout) :: A
    dimension :: A(2, nr)

    integer, intent(in) :: dod_u, dod_v
    dimension :: dod_u(ndu), dod_v(ndv)

    A(1, dod_u) = 0.d0 
    A(2, dod_v) = 0.d0 

end subroutine reset_dirichletbound2

subroutine reset_dirichletbound3(nr, A, ndu, ndv, ndw, dod_u, dod_v, dod_w)
    !! Set to 0 (Dirichlet condition) the values of an array using the dod indices in each dimension
    !! A is actually a vector arranged following each dimension [Au, Av, Aw]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ndu, ndv, ndw
    double precision, intent(inout) :: A
    dimension :: A(3, nr)

    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    A(1, dod_u) = 0.d0 
    A(2, dod_v) = 0.d0 
    A(3, dod_w) = 0.d0 

end subroutine reset_dirichletbound3

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

module solverplasticity2

    use matrixfreeplasticity
    use datastructure
    
    type cgsolver
        integer :: dimen = 2
        type(structure) :: disp_struct(2)
        end type cgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(solv%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        call mf_stiffness_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, table, mean)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_B_u, data_B_v, data_W_u, data_W_v
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2, solv%dimen)
        double precision, intent(in) :: mean
        dimension :: mean(solv%dimen, solv%dimen)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call init_2datastructure(solv%disp_struct(i), nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                    data_B_u, data_B_v, data_W_u, data_W_v)
            call update_datastructure(solv%disp_struct(i), solv%dimen, table(:, :, i))
            call eigendecomposition(solv%disp_struct(i), mean(i, :))
        end do
        
    end subroutine initializefastdiag

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(solv%dimen, nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, i
        double precision, allocatable, dimension(:) :: tmp, tmp2

        array_out  = 0.d0
        
        do i = 1, solv%dimen
            ! Compute (Uv x Uu)'.array_in
            nr_u = solv%disp_struct(i)%nrows(1)
            nr_v = solv%disp_struct(i)%nrows(2)
            allocate(tmp(nr_u*nr_v))
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, transpose(solv%disp_struct(i)%eigvec(1, 1:nr_u, 1:nr_u)), &
                    transpose(solv%disp_struct(i)%eigvec(2, 1:nr_v, 1:nr_v)), array_in(i, solv%disp_struct(i)%dof), tmp)

            tmp = tmp/solv%disp_struct(i)%Deigen

            ! Compute (Uv x Uu).array_tmp
            allocate(tmp2(nr_u*nr_v))
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, solv%disp_struct(i)%eigvec(1, 1:nr_u, 1:nr_u), &
                    solv%disp_struct(i)%eigvec(2, 1:nr_v, 1:nr_v), tmp, tmp2)
            array_out(i, solv%disp_struct(i)%dof) = tmp2            
            deallocate(tmp, tmp2)
        end do

    end subroutine applyfastdiag

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndu, ndv, dod_u, dod_v, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        integer, intent(in) :: ndu, ndv
        integer, intent(in) :: dod_u, dod_v
        dimension :: dod_u(ndu), dod_v(ndv)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), &
                        s(solv%dimen, nr_total), Ap(solv%dimen, nr_total), As(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound2(nr_total, r, ndu, ndv, dod_u, dod_v) 
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb  = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, p, Ap)
            call reset_dirichletbound2(nr_total, Ap, ndu, ndv, dod_u, dod_v) 
            call block_dot_product(solv%dimen, nr_total, Ap, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, s, As)
            call reset_dirichletbound2(nr_total, As, ndu, ndv, dod_u, dod_v) 
            call block_dot_product(solv%dimen, nr_total, As, s, prod)
            call block_dot_product(solv%dimen, nr_total, As, As, prod2)
            omega = prod/prod2
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndu, ndv, dod_u, dod_v, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        integer, intent(in) :: ndu, ndv
        integer, intent(in) :: dod_u, dod_v
        dimension :: dod_u(ndu), dod_v(ndv)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), & 
                        s(solv%dimen, nr_total), ptilde(solv%dimen, nr_total), Aptilde(solv%dimen, nr_total), &
                        Astilde(solv%dimen, nr_total), stilde(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound2(nr_total, r, ndu, ndv, dod_u, dod_v) 
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb  = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call reset_dirichletbound2(nr_total, ptilde, ndu, ndv, dod_u, dod_v) 
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ptilde, Aptilde)
            call reset_dirichletbound2(nr_total, Aptilde, ndu, ndv, dod_u, dod_v) 
            call block_dot_product(solv%dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call reset_dirichletbound2(nr_total, stilde, ndu, ndv, dod_u, dod_v) 
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, stilde, Astilde)
            call reset_dirichletbound2(nr_total, Astilde, ndu, ndv, dod_u, dod_v) 
            
            call block_dot_product(solv%dimen, nr_total, Astilde, s, prod)
            call block_dot_product(solv%dimen, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module solverplasticity2

module solverplasticity3

    use matrixfreeplasticity
    use datastructure

    type cgsolver
        integer :: dimen = 3
        type(structure) :: disp_struct(3)
    end type cgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat)  :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(solv%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        call mf_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, mean)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2, solv%dimen)
        double precision, intent(in) :: mean
        dimension :: mean(solv%dimen, solv%dimen)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call init_3datastructure(solv%disp_struct(i), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)
            call update_datastructure(solv%disp_struct(i), solv%dimen, table(:, :, i))
            call eigendecomposition(solv%disp_struct(i), mean(i, :))
        end do
        
    end subroutine initializefastdiag

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(solv%dimen, nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w, i
        double precision, allocatable, dimension(:) :: tmp, tmp2

        array_out = 0.d0

        do i = 1, solv%dimen
            ! Compute (Uw x Uv x Uu)'.array_in
            nr_u = solv%disp_struct(i)%nrows(1)
            nr_v = solv%disp_struct(i)%nrows(2)
            nr_w = solv%disp_struct(i)%nrows(3)
            allocate(tmp(nr_u*nr_v*nr_w))

            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, transpose(solv%disp_struct(i)%eigvec(1, 1:nr_u, 1:nr_u)), &
            transpose(solv%disp_struct(i)%eigvec(2, 1:nr_v, 1:nr_v)), transpose(solv%disp_struct(i)%eigvec(3, 1:nr_w, 1:nr_w)), &
            array_in(i, solv%disp_struct(i)%dof), tmp)

            tmp = tmp/solv%disp_struct(i)%Deigen

            ! Compute (Uw x Uv x Uu).array_tmp
            allocate(tmp2(nr_u*nr_v*nr_w))
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, solv%disp_struct(i)%eigvec(1, 1:nr_u, 1:nr_u), &
            solv%disp_struct(i)%eigvec(2, 1:nr_v, 1:nr_v), solv%disp_struct(i)%eigvec(3, 1:nr_w, 1:nr_w), tmp, tmp2)
            array_out(i, solv%disp_struct(i)%dof) = tmp2
            deallocate(tmp, tmp2)
        end do

    end subroutine applyfastdiag

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndu, ndv, ndw, dod_u, dod_v, dod_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        integer, intent(in) :: ndu, ndv, ndw
        integer, intent(in) :: dod_u, dod_v, dod_w
        dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), &
                        s(solv%dimen, nr_total), Ap(solv%dimen, nr_total), As(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound3(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call reset_dirichletbound3(nr_total, Ap, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(solv%dimen, nr_total, Ap, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            call reset_dirichletbound3(nr_total, As, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(solv%dimen, nr_total, As, s, prod)
            call block_dot_product(solv%dimen, nr_total, As, As, prod2)
            omega = prod/prod2
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndu, ndv, ndw, dod_u, dod_v, dod_w, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        integer, intent(in) :: ndu, ndv, ndw
        integer, intent(in) :: dod_u, dod_v, dod_w
        dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), & 
                        s(solv%dimen, nr_total), ptilde(solv%dimen, nr_total), Aptilde(solv%dimen, nr_total), &
                        Astilde(solv%dimen, nr_total), stilde(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound3(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde) 
            call reset_dirichletbound3(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call reset_dirichletbound3(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call block_dot_product(solv%dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call reset_dirichletbound3(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call reset_dirichletbound3(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            
            call block_dot_product(solv%dimen, nr_total, Astilde, s, prod)
            call block_dot_product(solv%dimen, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module solverplasticity3
