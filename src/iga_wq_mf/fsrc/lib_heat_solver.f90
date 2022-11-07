module heat_solver

    use heat_spmf
    type cgsolver
        logical :: withscaling = .false., withdiag = .false.
        integer :: matrixfreetype = 1
        double precision, dimension(:), allocatable :: factor
        double precision, dimension(:), pointer :: diag
    end type cgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        type(thermomat), pointer :: mat
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

            call mf_wq_get_cu_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_wq_get_ku_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.3) then

            call mf_wq_get_kcu_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else 
            stop 'function not defined'
        end if

    end subroutine matrixfree_spMdV

    subroutine setup_eigendiag(solv, nr_total, diag)

        use omp_lib
        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        integer, intent(in) :: nr_total
        double precision, target, intent(in) :: diag
        dimension :: diag(nr_total)

        solv%withdiag = .true.
        solv%diag => diag
        
    end subroutine setup_eigendiag

    subroutine setup_FDscaling(solv, nr_total, factor_up, factor_down)

        use omp_lib
        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: factor_up, factor_down
        dimension :: factor_up(nr_total), factor_down(nr_total)

        ! Local data
        ! ----------
        integer :: i, nb_tasks
        double precision, target :: factor
        dimension :: factor(nr_total)
        
        solv%withscaling = .true.
        
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
        do i = 1, nr_total
            factor(i) = sqrt(factor_up(i)/factor_down(i))
        end do  
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

        allocate(solv%factor(nr_total))
        solv%factor = factor

    end subroutine setup_FDscaling

    subroutine fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver), pointer :: solv
        integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
        double precision, intent(in) :: U_u, U_v, U_w, array_in
        dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: i, nb_tasks
        double precision :: array_temp, dummy
        dimension :: array_temp(nr_total), dummy(nr_total)

        dummy = array_in
        if (solv%withscaling) then
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
            do i = 1, nr_total
                dummy(i) = solv%factor(i)*dummy(i) 
            end do  
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
        end if  
    
        ! Compute (Uw x Uv x Uu)'.array_in
        call sumproduct3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        transpose(U_u), transpose(U_v), transpose(U_w), dummy, array_temp)
        
        if (solv%withdiag) then
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
            do i = 1, nr_total
                array_temp(i) = array_temp(i)/solv%diag(i)
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        end if
    
        ! Compute (Uw x Uv x Uu).array_temp
        call sumproduct3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, dummy)

        if (solv%withscaling) then
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
            do i = 1, nr_total
                dummy(i) = solv%factor(i)*dummy(i) 
            end do  
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
        end if
        array_out = dummy
        
    end subroutine fast_diagonalization

    subroutine CG(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        type(thermomat), pointer :: mat
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

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, normb
        double precision :: r, p, Ap
        dimension :: r(nr_total), p(nr_total), Ap(nr_total)
        integer :: iter

        x = 0.d0; r = b; normb = maxval(abs(r))
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return 

        rsold = dot_product(r, r); p = r
        
        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            
            alpha = rsold/dot_product(p, Ap)
            x = x + alpha * p
            r = r - alpha * Ap

            resPCG(iter+1) = maxval(abs(r))/normb
            if (resPCG(iter+1).le.threshold) exit
        
            rsnew = dot_product(r, r)
            p = r + rsnew/rsold * p
            rsold = rsnew
        end do

    end subroutine CG

    subroutine PCG(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        type(thermomat), pointer :: mat
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

        double precision, intent(in) :: U_u, U_v, U_w
        dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, normb
        double precision :: r, p, Ap, z
        dimension :: r(nr_total), p(nr_total), Ap(nr_total), z(nr_total)
        integer :: iter

        x = 0.d0; r = b; normb = maxval(abs(r))
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        call fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, r, z)
        rsold = dot_product(r, z); p = z
        
        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)

            alpha = rsold/dot_product(p, Ap)
            x = x + alpha * p
            r = r - alpha * Ap

            resPCG(iter+1) = maxval(abs(r))/normb
            if (resPCG(iter+1).le.threshold) exit       

            call fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, r, z)
            rsnew = dot_product(r, z)

            p = z + rsnew/rsold * p
            rsold = rsnew
        end do

    end subroutine PCG

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        type(thermomat), pointer :: mat
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

        x = 0.d0; r = b; rhat = r; p = r
        rsold = dot_product(r, rhat); normb = maxval(abs(r))
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
            omega = dot_product(As, s)/dot_product(As, As)
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = maxval(abs(r))/normb
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
                        data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver), pointer :: solv
        type(thermomat), pointer :: mat
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

        double precision, intent(in) :: U_u, U_v, U_w
        dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w)

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

        x = 0.d0; r = b; rhat = r; p = r
        rsold = dot_product(r, rhat); normb = maxval(abs(r))
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, p, ptilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, s, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = maxval(abs(r))/normb
            if (resPCG(iter+1).le.threshold) exit
    
            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module heat_solver