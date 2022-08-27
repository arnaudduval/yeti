! ==========================
! module: Penalty or Lagrange multipliers method to solve Dirichlet problem
! author: Joaquin Cornejo
! 
! LAGRANGE:
! In this module, one could find functions to solve Au = b using Lagrange multipliers
! That is, we want to solve Mx = F, with Bx = g (Dirichlet condition). In matrix representation:
! [ M  B'    [x  = [F
!   B  0 ]    y]    g]
! or  A       u  =  b
! with y, the lagrange multiplier vector. 
! As preconditioner we use: 
! P = [ MM  0
!       0   SS]
! Where MM is an approximation of M and SS is an approximation of S 
! S is the Schur complement, S = B M^-1 B'
! In practice, to compute MM we use fast diagonalisation and 
! to compute SS we use the inverse of the diagonal of MM (in Dirichlet nodes)
!
! PENALTY:
! We want to solve Mx = F, with Bx = g (Dirichlet condition). In matrix representation:
! [ Knn  Knd      [xn  = [Fn
!   Kdn  Kdd+P]    xd]    Fd+Pxd]
! or  A           u  =  b
! ==========================

subroutine create_block_B(nr, nc, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT)
    !! Creates B matrix of the new matrix. 
    !! where the new matrix is [A, BT; B, 0]. 
    !! The block B is a matrix of size nrxnc
    !! Returns B and BT in CSR format
    !! In this special matrix, the number of non zero values is equal to the number of rows 

    implicit none
    ! Input/output data
    ! -----------------
    integer, intent(in) :: nr, nc, dod
    dimension :: dod(nr)

    integer, intent(out) :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(nr+1), indj_B(nr), indi_BT(nc+1), indj_BT(nr)
    double precision, intent(out) :: B, BT
    dimension :: B(nr), BT(nr)
    
    ! Local data
    ! -----------------
    integer :: i, indi_coo, indj_coo
    dimension :: indi_coo(nr), indj_coo(nr)
    double precision :: data_coo
    dimension :: data_coo(nr, 1)

    ! Get COO format
    do i = 1, nr
        indi_coo(i) = i
        indj_coo(i) = dod(i)
        data_coo(i, 1) = 1.d0 
    end do

    ! Get B in CSR format
    call coo2csr(1, nr, nr, data_coo, indi_coo, indj_coo, B, indj_B, indi_B)

    ! Get B transpose in CSR format
    call coo2csr(1, nc, nr, data_coo, indj_coo, indi_coo, BT, indj_BT, indi_BT)

end subroutine create_block_B

! STEADY HEAT TRANSFER CASE: M = K + B' P B and F = f + B' P g
! -------------------------------------------------------------
! With Lagrange multipliers method

subroutine mf_wq_get_Au_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_B, indj_B, indi_BT, indj_BT, B, BT, &
                                penalty, ndod, array_in, array_out)
    !! Computes A.u in steady heat 3D case with lagrange multipliers, where 
    !! A = [ M  B'   
    !!       B  0 ] with M = K + B' P B and P a penalty diagonal matrix = p Identity
    !! Array input is composed of 2 blocks = [x y], then the output will be [M.x + B'.y, B.x]
    !! This, if we replace, [K.x + B'.(p B.x + y), B.x] 
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    double precision, intent(in) :: penalty
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)

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

    integer, intent(in) :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nr_total+1), indj_BT(ndod)
    double precision :: B, BT
    dimension :: B(ndod), BT(ndod)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total+ndod)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total+ndod)

    ! Local data 
    ! ------------------
    double precision :: x, y
    dimension :: x(nr_total), y(ndod)
    
    double precision :: Kx, Bx, L, BTL
    dimension :: Kx(nr_total), Bx(ndod), L(ndod), BTL(nr_total)

    ! Assign blocks
    x = array_in(1:nr_total)
    y = array_in(nr_total+1:nr_total+ndod)

    ! Compute K.x
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, x, Kx)

    ! Compute B.x                   
    call spMdotdV(ndod, nr_total, ndod, indi_B, indj_B, B, x, Bx)

    ! Compute L = p B.x + y
    L = penalty*Bx + y
    
    ! Compute K.x + B'.L
    call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, L, BTL)

    ! Set values
    array_out(1:nr_total) = Kx + BTL
    array_out(nr_total+1:nr_total+ndod) = Bx

end subroutine mf_wq_get_Au_shlm_3d

subroutine fd_shlm_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, pardiag, ndod, dod, array_in, array_out)
    !! Solves Z.s = r (Z is the preconditioner) in steady heat 3D case with lagrange multipliers, where 
    !!     [ MM  0    [s1  = [r1
    !!       0   SS ]  s2]    r2]
    !! To keep it simple, we solve independently s1 = MM^-1 r1 and we set s2 = 0
    !! MM is an approximation of M = K + B' P B and P a penalty diagonal matrix. To compute MM^-1 we use fast diagonalization
    !! SS is an approximation of Schur complement S = B MM^-1 B'. With B a block matrix of 0 and 1, S is really a block matrix
    !! of Dirichlet nodes. Then we approximate S as the inverse of diagonal of MM (in Dirichlet nodes)
    !! Indices must be in CSR format
    !! ------------------------------To review SS

    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: ndod, nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag, pardiag
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), eigendiag(nr_total), pardiag(nr_total)

    integer, intent(in) :: dod
    dimension :: dod(ndod)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total+ndod)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total+ndod)

    ! Local data
    ! ---------------
    integer :: i
    double precision :: r1, r2, s1, s2
    dimension :: r1(nr_total), r2(ndod), s1(nr_total), s2(ndod)

    ! Assign blocks 
    r1 = array_in(1:nr_total)
    r2 = array_in(nr_total+1:nr_total+ndod)

    ! Solve firts block
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, r1, s1)

    ! Solve second block
    do i = 1, ndod
        s2(i) = r2(i)*pardiag(dod(i)) ! to review
    end do

    ! Assign values
    array_out(1:nr_total) = s1
    array_out(nr_total+1:nr_total+ndod) = s2

end subroutine fd_shlm_3d

subroutine mf_wq_solve_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, x, energy)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition), where M = M = K + B' P B and P a penalty diagonal matrix
    !! and F = f + B' P g. Using Lagrange multipliers, the linear system is:
    !! [ M  B'    [x  = [F
    !!   B  0 ]    y]    g]
    !! CSR FORMAT
                        
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d=3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: table, dod
    dimension :: table(d, 2), dod(ndod)
    integer, intent(in) :: nbIter
    double precision, intent(in) :: epsilon, f, g
    dimension :: f(nr_total), g(ndod)
    
    double precision, intent(out) :: x, energy
    dimension :: x(nr_total), energy(nbIter)

    ! Local data
    ! ------------------
    double precision :: dummy_tol
    ! Pre / Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, xt, Ax
    dimension ::    r(nr_total+ndod), rhat(nr_total+ndod), p(nr_total+ndod), s(nr_total+ndod), ptilde(nr_total+ndod), &
                    Aptilde(nr_total+ndod), Astilde(nr_total+ndod), stilde(nr_total+ndod), xt(nr_total+ndod), Ax(nr_total)
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    ! Block B
    integer :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nr_total+1), indj_BT(ndod)
    double precision :: B, BT, BTPg
    dimension :: B(ndod), BT(ndod), BTPg(nr_total)
    double precision :: penalty 

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! --------------------------------------------
    ! DIAGONAL DECOMPOSITION
    ! -------------------------------------------- 
    ! Initialize
    dummy_tol = epsilon
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    
    do iter = 1, 2
        call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                    Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
    end do

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), table(1, :), D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), table(2, :), D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), table(3, :), D_w, U_w, Kdiag_w, Mdiag_w)   
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    allocate(Deigen(nr_total))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
    deallocate(I_u, I_v, I_w)

    ! Find diagonal of preconditioner
    allocate(Dparametric(nr_total))
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                Kdiag_u, Kdiag_v, Kdiag_w, Dparametric)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    ! Find diagonal of real matrix (K in this case)
    allocate(Dphysical(nr_total))
    call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)

    ! Define penalty
    penalty = max(100.d0, 100.d0*maxval(Dphysical))

    if (nbIter.gt.0) then
        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        ! The system to solve is A u = b
        ! Where b = [f + B' P g
        !            g ]
        ! So, to solve this system we initialize r = b - A u, if u = 0, r = b

        ! Create block B
        call create_block_B(ndod, nr_total, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT) 

        ! Define r
        call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, penalty*g, BTPg)
        r(1:nr_total) = f + BTPg
        r(nr_total+1:nr_total+ndod) = g

        rhat = r; p = r
        rsold = dot_product(r, rhat)
        energy(1) = 0.d0 ! Energy : 0.5 u' A u - u' f  

        do iter = 1, nbIter

            call fd_shlm_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, ndod, dod, p, ptilde)
            call mf_wq_get_Au_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_shlm_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, ndod, dod, s, stilde)
            call mf_wq_get_Au_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            xt = xt + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            ! Save result
            x = xt(1:nr_total)

            call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)
            
            energy(iter) = 0.5 * dot_product(x, Ax) - dot_product(x, f) ! Energy : 0.5 u' A u - u' f       
            ! if (energy(iter+1).le.epsilon) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end if

end subroutine mf_wq_solve_shlm_3d

! With penalty method

subroutine mf_wq_get_Au_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_B, indj_B, indi_BT, indj_BT, B, BT, &
                                penalty, ndod, array_in, array_out)
    !! Computes A.u in steady heat 3D case with penalty method, where 
    !! A = K + B' P B, with P a penalty diagonal matrix = p Identity
    !! The output will be Au = Ku + p B' B u 
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    double precision, intent(in) :: penalty
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)

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

    integer, intent(in) :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nr_total+1), indj_BT(ndod)
    double precision :: B, BT
    dimension :: B(ndod), BT(ndod)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ------------------   
    double precision :: Ku, Bu, BTBu
    dimension :: Ku(nr_total), Bu(ndod), BTBu(nr_total)

    ! Compute K.u
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, Ku)

    ! Compute B.u                   
    call spMdotdV(ndod, nr_total, ndod, indi_B, indj_B, B, array_in, Bu)
    
    ! Compute B' B u 
    call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, Bu, BTBu)

    ! Set values
    array_out = Ku + penalty*BTBu

end subroutine mf_wq_get_Au_shp_3d

subroutine fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, array_in, array_out)
    !! Solves MM.s = r (MM is the preconditioner) in steady heat 3D case with lagrange multipliers, where 
    !! MM is an approximation of M = K + B' P B and P a penalty diagonal matrix. 
    !! With B a block matrix of 0 and 1. To compute MM^-1 we use fast diagonalization.
    !! Indices must be in CSR format

    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), eigendiag(nr_total)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! By now, we test this approximation. It could change later
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, array_in, array_out)

end subroutine fd_shp_3d

subroutine mf_wq_solve_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, x, energy)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition), where M = M = K + B' P B and P a penalty diagonal matrix
    !! and F = f + B' P g. We use penalty method
    !! CSR FORMAT
                        
    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d=3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
    double precision, intent(in) :: coefs
    dimension :: coefs(d, d, nc_total)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: table, dod
    dimension :: table(d, 2), dod(ndod)
    integer, intent(in) :: nbIter
    double precision, intent(in) :: epsilon, f, g
    dimension :: f(nr_total), g(ndod)
    
    double precision, intent(out) :: x, energy
    dimension :: x(nr_total), energy(nbIter)

    ! Local data
    ! ------------------
    double precision :: dummy_tol
    ! Pre / Conjugate gradient algoritm
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, dummy, Ax
    dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), ptilde(nr_total), &
                    Aptilde(nr_total), Astilde(nr_total), stilde(nr_total), dummy(nr_total), Ax(nr_total)
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    ! Block B
    integer :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nr_total+1), indj_BT(ndod)
    double precision :: B, BT, BTPg
    dimension :: B(ndod), BT(ndod), BTPg(nr_total)
    double precision :: penalty 

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize B transpose in CSR format
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! --------------------------------------------
    ! DIAGONAL DECOMPOSITION
    ! -------------------------------------------- 
    ! Initialize
    dummy_tol = epsilon
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0
    
    do iter = 1, 2
        call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                    Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
    end do

    ! --------------------------------------------
    ! EIGEN DECOMPOSITION
    ! -------------------------------------------- 
    allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), table(1, :), D_u, U_u, Kdiag_u, Mdiag_u)

    allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), table(2, :), D_v, U_v, Kdiag_v, Mdiag_v)    

    allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), table(3, :), D_w, U_w, Kdiag_w, Mdiag_w)   
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

    ! Find diagonal of eigen values
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    allocate(Deigen(nr_total))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
    deallocate(I_u, I_v, I_w)

    ! Find diagonal of preconditioner
    allocate(Dparametric(nr_total))
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                Kdiag_u, Kdiag_v, Kdiag_w, Dparametric)
    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    ! Find diagonal of real matrix (K in this case)
    allocate(Dphysical(nr_total))
    call wq_find_conductivity_diagonal_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)

    ! Define penalty
    penalty = max(100.d0, 10000.d0*maxval(Dphysical))
    Dphysical(dod) = Dphysical(dod) + penalty 

    if (nbIter.gt.0) then
        ! -------------------------------------------
        ! Preconditioned Conjugate Gradient algorithm
        ! -------------------------------------------
        ! The system to solve is A u = b
        ! Where b = f + B' P g
        ! So, to solve this system we initialize r = b - A u, if u = 0, r = b

        ! Create block B
        call create_block_B(ndod, nr_total, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT) 

        ! Define r
        call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, penalty*g, BTPg)
        r = f + BTPg

        rhat = r; p = r
        rsold = dot_product(r, rhat)
        energy(1) = 0.d0 ! Energy : 0.5 u' A u - u' f

        do iter = 1, nbIter

            dummy = p
            call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy) 
            call fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, ptilde)
            call fd_sqr_scaling(nr_total, Dparametric, Dphysical, ptilde)  

            call mf_wq_get_Au_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            dummy = s
            call fd_sqr_scaling(nr_total, Dparametric, Dphysical, dummy)  
            call fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, dummy, stilde)
            call fd_sqr_scaling(nr_total, Dparametric, Dphysical, stilde) 

            call mf_wq_get_Au_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde 
            r = s - omega*Astilde   

            call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)
            
            energy(iter) = 0.5 * dot_product(x, Ax) - dot_product(x, f) ! Energy : 0.5 u' A u - u' f       
            ! if (RelRes(iter+1).le.epsilon) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end if

end subroutine mf_wq_solve_shp_3d