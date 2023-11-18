subroutine reset_dirichletbound1(nr, A, ndod, dod)
    !! Set to 0 (Dirichlet condition) the values of an array using the dod indices in each dimension
    !! A is actually a vector arranged following each dimension [Au, Av, Aw]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ndod
    double precision, intent(inout) :: A
    dimension :: A(nr)

    integer, intent(in) :: dod
    dimension :: dod(ndod)

    A(dod) = 0.d0 

end subroutine reset_dirichletbound1

module solverheat2

    use matrixfreeheat
    use datastructure
    type cgsolver
        logical :: withdiag = .true.
        integer :: matrixfreetype = 1, dimen = 2
        type(structure) :: temp_struct
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
        type(thermomat) :: mat
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
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        if (solv%matrixfreetype.eq.1) then

            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        else 
            call mf_uv_gradugradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, &
                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                data_B_u, data_B_v, data_W_u, data_W_v, table_dirichlet, mean)

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
        logical, intent(in) :: table_dirichlet
        dimension :: table_dirichlet(solv%dimen, 2)
        double precision, intent(in) :: mean
        dimension :: mean(solv%dimen)

        call init_2datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v)
        call update_datastructure(solv%temp_struct, solv%dimen, table_dirichlet)
        call eigendecomposition(solv%temp_struct, mean)
    
    end subroutine initializefastdiag

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v
        double precision, allocatable, dimension(:) :: tmp, tmp2

        ! Compute (Uw x Uv x Uu)'.array_in
        nr_u = solv%temp_struct%nrows(1)
        nr_v = solv%temp_struct%nrows(2)
        allocate(tmp(nr_u*nr_v))
        call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, transpose(solv%temp_struct%eigvec(1, 1:nr_u, 1:nr_u)), &
                transpose(solv%temp_struct%eigvec(2, 1:nr_v, 1:nr_v)), array_in(solv%temp_struct%dof), tmp)
        
        if (solv%withdiag) then
            tmp = tmp/solv%temp_struct%Deigen
        end if
    
        ! Compute (Uv x Uu).array_tmp
        allocate(tmp2(nr_u*nr_v))
        call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, solv%temp_struct%eigvec(1, 1:nr_u, 1:nr_u), &
                solv%temp_struct%eigvec(2, 1:nr_v, 1:nr_v), tmp, tmp2)
        array_out(solv%temp_struct%dof) = tmp2            
        deallocate(tmp, tmp2)

    end subroutine applyfastdiag

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), &
                        s(nr_total), Ap(nr_total), As(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound1(nr_total, r, ndod, dod) 
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, p, Ap)
            call reset_dirichletbound1(nr_total, Ap, ndod, dod)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, s, As)
            call reset_dirichletbound1(nr_total, As, ndod, dod)
            omega = dot_product(As, s)/dot_product(As, As)
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndod, dod, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
                        ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound1(nr_total, r, ndod, dod)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call reset_dirichletbound1(nr_total, ptilde, ndod, dod)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ptilde, Aptilde)
            call reset_dirichletbound1(nr_total, Aptilde, ndod, dod)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call reset_dirichletbound1(nr_total, stilde, ndod, dod)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, stilde, Astilde)
            call reset_dirichletbound1(nr_total, Astilde, ndod, dod)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
    
            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

    subroutine LOBPCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ndod, dod, ishigher, nbIterPCG, threshold, eigenvec, eigenval)
        !! Using LOBPCG algorithm to compute the stability of the transient heat problem
        
        implicit none
        ! Input / output data
        ! -------------------
        integer, parameter :: d = 3
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        logical, intent(in) :: ishigher
        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(nr_total)

        ! Local data
        ! ----------
        integer :: k, ii
        double precision, dimension(d, nr_total) :: RM1, RM2, RM3
        double precision, dimension(d, d) :: AA1, BB1
        double precision, dimension(d) :: delta
        double precision, dimension(nr_total) :: u, v, g, gtil, p, tmp
        double precision :: q, norm
        double precision, allocatable, dimension(:) :: ll
        double precision, allocatable, dimension(:, :) :: qq
    
        call random_number(eigenvec)
        call reset_dirichletbound1(nr_total, eigenvec, ndod, dod)
        norm = norm2(eigenvec)
        eigenvec = eigenvec/norm

        call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, u)
        call reset_dirichletbound1(nr_total, u, ndod, dod)

        q = sqrt(dot_product(eigenvec, u))
        eigenvec = eigenvec/q; u = u/q
        call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, v)
        call reset_dirichletbound1(nr_total, v, ndod, dod)
        eigenval = dot_product(eigenvec, v)
        p = 0.d0
        norm = 1.d0

        do k = 1, nbIterPCG
            if (norm.le.threshold) return
    
            g = v - eigenval*u
            norm = norm2(g)
            call applyfastdiag(solv, nr_total, g, gtil)
            call reset_dirichletbound1(nr_total, gtil, ndod, dod)
            g = gtil

            RM1(1, :) = eigenvec; RM1(2, :) = -g; RM1(3, :) = p
            RM2(1, :) = v; RM3(1, :) = u;
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v,  -g, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM2(2, :) = tmp
            
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, p, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM2(3, :) = tmp
                        
            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, -g, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM3(2, :) = tmp

            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, p, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM3(3, :) = tmp
            
            call rayleigh_submatrix(d, nr_total, RM1, RM2, AA1); AA1 = 0.5d0*(AA1 + transpose(AA1))
            call rayleigh_submatrix(d, nr_total, RM1, RM3, BB1); BB1 = 0.5d0*(BB1 + transpose(BB1))

            if (k.eq.1) then
                allocate(ll(d-1), qq(d-1, d-1))
                call compute_eigdecomp_pdr(size(ll), AA1(:d-1, :d-1), BB1(:d-1, :d-1), ll, qq)
            else
                allocate(ll(d), qq(d, d))
                call compute_eigdecomp_pdr(size(ll), AA1, BB1, ll, qq)
            end if
            
            if (ishigher) then 
                eigenval = maxval(ll); ii = maxloc(ll, dim=1)
            else
                eigenval = minval(ll); ii = minloc(ll, dim=1)
            end if

            delta = 0.d0
            if (k.eq.1) then
                delta(:2) = qq(:, ii)
            else
                delta = qq(:, ii)
            end if
    
            p = -g*delta(2) + p*delta(3)
            eigenvec = eigenvec*delta(1) + p
            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v,&
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, eigenvec, u)
            call reset_dirichletbound1(nr_total, u, ndod, dod)
            
            q = sqrt(dot_product(eigenvec, u))
            eigenvec = eigenvec/q; u = u/q
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, v)
            call reset_dirichletbound1(nr_total, v, ndod, dod)
            
            norm = norm2(g)
            deallocate(ll, qq)
        end do

    end subroutine LOBPCGSTAB

end module solverheat2

module solverheat3

    use matrixfreeheat
    use datastructure
    type cgsolver
        logical :: withdiag = .true.
        integer :: matrixfreetype = 1, dimen = 3
        type(structure) :: temp_struct
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
        type(thermomat) :: mat
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
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        if (solv%matrixfreetype.eq.1) then

            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.3) then

            call mf_uv_gradugradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else 
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, Kmean)

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
        dimension :: table(solv%dimen, 2)
        double precision, intent(in) :: Kmean
        dimension :: Kmean(solv%dimen)

        call init_3datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)
        call update_datastructure(solv%temp_struct, solv%dimen, table)
        call eigendecomposition(solv%temp_struct, Kmean)
    
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
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w
        double precision, allocatable, dimension(:) :: tmp, tmp2

        array_out = 0.d0

        ! Compute (Uw x Uv x Uu)'.array_in
        nr_u = solv%temp_struct%nrows(1)
        nr_v = solv%temp_struct%nrows(2)
        nr_w = solv%temp_struct%nrows(3)
        allocate(tmp(nr_u*nr_v*nr_w))

        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, transpose(solv%temp_struct%eigvec(1, 1:nr_u, 1:nr_u)), &
        transpose(solv%temp_struct%eigvec(2, 1:nr_v, 1:nr_v)), transpose(solv%temp_struct%eigvec(3, 1:nr_w, 1:nr_w)), &
        array_in(solv%temp_struct%dof), tmp)

        tmp = tmp/solv%temp_struct%Deigen

        ! Compute (Uw x Uv x Uu).array_tmp
        allocate(tmp2(nr_u*nr_v*nr_w))
        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, solv%temp_struct%eigvec(1, 1:nr_u, 1:nr_u), &
        solv%temp_struct%eigvec(2, 1:nr_v, 1:nr_v), solv%temp_struct%eigvec(3, 1:nr_w, 1:nr_w), tmp, tmp2)
        array_out(solv%temp_struct%dof) = tmp2
        deallocate(tmp, tmp2)

    end subroutine applyfastdiag

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), &
                        s(nr_total), Ap(nr_total), As(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound1(nr_total, r, ndod, dod) 
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call reset_dirichletbound1(nr_total, Ap, ndod, dod)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            call reset_dirichletbound1(nr_total, As, ndod, dod)
            omega = dot_product(As, s)/dot_product(As, As)
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndod, dod, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
                        ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call reset_dirichletbound1(nr_total, r, ndod, dod)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call reset_dirichletbound1(nr_total, ptilde, ndod, dod)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call reset_dirichletbound1(nr_total, Aptilde, ndod, dod)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call reset_dirichletbound1(nr_total, stilde, ndod, dod)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call reset_dirichletbound1(nr_total, Astilde, ndod, dod)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
    
            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

    subroutine LOBPCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ndod, dod, ishigher, nbIterPCG, threshold, eigenvec, eigenval)
        !! Using LOBPCG algorithm to compute the stability of the transient heat problem
        
        implicit none
        ! Input / output data
        ! -------------------
        integer, parameter :: d = 3
        type(cgsolver) :: solv
        type(thermomat) :: mat
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

        integer, intent(in) :: ndod
        integer, intent(in) :: dod
        dimension :: dod(ndod)

        logical, intent(in) :: ishigher
        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(nr_total)

        ! Local data
        ! ----------
        integer :: k, ii
        double precision, dimension(d, nr_total) :: RM1, RM2, RM3
        double precision, dimension(d, d) :: AA1, BB1
        double precision, dimension(d) :: delta
        double precision, dimension(nr_total) :: u, v, g, gtil, p, tmp
        double precision :: q, norm
        double precision, allocatable, dimension(:) ::  ll
        double precision, allocatable, dimension(:, :) ::  qq

        call random_number(eigenvec)
        call reset_dirichletbound1(nr_total, eigenvec, ndod, dod)
        norm = norm2(eigenvec)
        eigenvec = eigenvec/norm

        call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, u)
        call reset_dirichletbound1(nr_total, u, ndod, dod)
        
        q = sqrt(dot_product(eigenvec, u))
        eigenvec = eigenvec/q; u = u/q
        call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, v)
        call reset_dirichletbound1(nr_total, v, ndod, dod)
        
        eigenval = dot_product(eigenvec, v)
        p = 0.d0
        norm = 1.d0

        do k = 1, nbIterPCG
            if (norm.le.threshold) return
    
            g = v - eigenval*u
            norm = norm2(g)
            call applyfastdiag(solv, nr_total, g, gtil)
            call reset_dirichletbound1(nr_total, gtil, ndod, dod)
            g = gtil

            RM1(1, :) = eigenvec; RM1(2, :) = -g; RM1(3, :) = p
            RM2(1, :) = v; RM3(1, :) = u;
            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, -g, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM2(2, :) = tmp

            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM2(3, :) = tmp

            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, -g, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM3(2, :) = tmp
            
            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, tmp)
            call reset_dirichletbound1(nr_total, tmp, ndod, dod)
            RM3(3, :) = tmp
            
            call rayleigh_submatrix(d, nr_total, RM1, RM2, AA1); AA1 = 0.5d0*(AA1 + transpose(AA1))
            call rayleigh_submatrix(d, nr_total, RM1, RM3, BB1); BB1 = 0.5d0*(BB1 + transpose(BB1))

            if (k.eq.1) then
                allocate(ll(d-1), qq(d-1, d-1))
                call compute_eigdecomp_pdr(size(ll), AA1(:d-1, :d-1), BB1(:d-1, :d-1), ll, qq)
            else
                allocate(ll(d), qq(d, d))
                call compute_eigdecomp_pdr(size(ll), AA1, BB1, ll, qq)
            end if
            
            if (ishigher) then 
                eigenval = maxval(ll); ii = maxloc(ll, dim=1)
            else
                eigenval = minval(ll); ii = minloc(ll, dim=1)
            end if

            delta = 0.d0
            if (k.eq.1) then
                delta(:2) = qq(:, ii)
            else
                delta = qq(:, ii)
            end if
    
            p = -g*delta(2) + p*delta(3)
            eigenvec = eigenvec*delta(1) + p
            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, eigenvec, u)
            call reset_dirichletbound1(nr_total, u, ndod, dod)
            
            q = sqrt(dot_product(eigenvec, u))
            eigenvec = eigenvec/q; u = u/q
            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, v)
            call reset_dirichletbound1(nr_total, v, ndod, dod)
            
            norm = norm2(g)
            deallocate(ll, qq)
        end do

    end subroutine LOBPCGSTAB

end module solverheat3