module stsolverheat3

    use matrixfreestheat
    use datastructure
    type stcgsolver
        integer :: dimen = 4
        type(structure) :: temp_struct
    end type stcgsolver

contains

    subroutine matrixfree_spMdV(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t, indj_T_u, indj_T_v, indj_T_w, indj_T_t
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
        
        ! Local data
        ! ----------
        double precision :: tmp(nr_total)

        call mf_partialt_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_w, data_W_t, array_in, array_out)

        call mf_gradx_u_gradx_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_w, data_W_t, array_in, tmp)

        array_out = array_out + tmp

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, &
                data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t, table)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), data_B_t(nnz_t, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2)
        
        ! Local data
        ! ----------
        double precision :: dummymean(solv%dimen)

        solv%temp_struct%isspacetime = .true.; dummymean = 1.d0
        call init_4datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t)
        call update_datastructure(solv%temp_struct, solv%dimen, table)
        call space_eigendecomposition(solv%temp_struct, solv%dimen, dummymean)
        !! diageigvalues has not been created
    
    end subroutine initializefastdiag

    subroutine solve_schurtriangular__(solv, nr, coefs, b, x)
        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr
        double precision, intent(in) :: coefs(2)
        double precision, intent(in) :: b(nr)
        double precision, intent(out) :: x(nr)
        
        ! Local data
        ! ----------
        double complex :: uptrg(nr, nr), Vb(nr), sol(nr), Usol(nr)

        Vb = matmul(transpose(conjg(solv%temp_struct%Lschur_tm)), dcmplx(b))
        uptrg = dcmplx(coefs(1))*solv%temp_struct%Luptr_tm + dcmplx(coefs(2))*solv%temp_struct%Ruptr_tm
        call solve_complex_uppertriangular_system(nr, uptrg, Vb, sol)
        Usol = matmul(solv%temp_struct%Rschur_tm, sol)
        x = realpart(Usol)

    end subroutine solve_schurtriangular__

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w, nr_t, i, j, k, l
        double precision, allocatable, dimension(:, :) :: identity
        double precision, allocatable, dimension(:) :: tmp, tmp1, tmp2, tmp3, btmp, sol
        integer, allocatable, dimension(:, :, :, :) :: dof
        integer, allocatable, dimension(:) :: indices 
        double precision :: eigval
        integer :: nrows(solv%dimen)

        array_out = 0.d0

        nrows = solv%temp_struct%nrows
        nr_u = nrows(1); nr_v = nrows(2); nr_w = nrows(3); nr_t = nrows(4)
        dof  = reshape(solv%temp_struct%dof, shape=(/nr_u, nr_v, nr_w, nr_t/))

        ! Compute (Ut x Uw x Uv x Uu)'.array_in
        allocate(tmp(nr_u*nr_v*nr_w*nr_t), identity(nr_t, nr_t))

        call sumfacto4d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        transpose(solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                        identity, array_in(solv%temp_struct%dof), tmp)

        allocate(tmp1(nr_total))
        tmp1 = 0.d0; tmp1(solv%temp_struct%dof) = tmp
        deallocate(tmp)

        allocate(btmp(nr_t), sol(nr_t), indices(nr_t), tmp2(nr_total))
        tmp2 = 0.d0
        do k = 1, nr_w
            do j = 1, nr_v
                do i = 1, nr_u
                    l = i + (j - 1)*nr_u + (k - 1)*nr_u*nr_v
                    eigval = solv%temp_struct%diageigval_sp(l)
                    indices = dof(i, j, k, :)
                    btmp = tmp1(indices)
                    call solve_schurtriangular__(solv, nr_t, (/1.d0, eigval/), btmp, sol)
                    tmp2(indices) = sol
                end do
            end do
        end do  
        deallocate(tmp1)

        ! Compute (Ut x Uw x Uv x Uu).array_tmp
        allocate(tmp3(nr_u*nr_v*nr_w*nr_t))
        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), &
                        identity, tmp2(solv%temp_struct%dof), tmp3)
        array_out(solv%temp_struct%dof) = tmp3
        deallocate(tmp2, tmp3)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(nnz)

        call set2zero(solv%temp_struct, nnz, array)

    end subroutine clear_dirichlet

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_w, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

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
        call clear_dirichlet(solv, nr_total, r) 
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, p, Ap)

            call clear_dirichlet(solv, nr_total, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, s, As)
            call clear_dirichlet(solv, nr_total, As)
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

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_w, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

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
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_spMdV(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
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

end module stsolverheat3