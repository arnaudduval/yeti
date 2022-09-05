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

subroutine create_block_B(nc, ndod, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT)
    !! Creates B matrix of the new matrix. 
    !! where the new matrix is [A, BT; B, 0]. 
    !! The block B is a matrix of size nrxnc
    !! Returns B and BT in CSR format
    !! In this special matrix, the number of non zero values is equal to the number of rows (nr = ndod)

    implicit none
    ! Input/output data
    ! -----------------
    integer, intent(in) :: nc, ndod, dod
    dimension :: dod(ndod)

    integer, intent(out) :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nc+1), indj_BT(ndod)
    double precision, intent(out) :: B, BT
    dimension :: B(ndod), BT(ndod)
    
    ! Local data
    ! -----------------
    integer :: i, nr, indi_coo, indj_coo
    dimension :: indi_coo(ndod), indj_coo(ndod)
    double precision :: data_coo
    dimension :: data_coo(ndod, 1)

    ! Initialize
    nr = ndod

    ! Get COO format
    do i = 1, nr ! nr = ndod
        indi_coo(i) = i
        indj_coo(i) = dod(i)
        data_coo(i, 1) = 1.d0 
    end do

    ! Get B in CSR format
    call coo2csr(1, nr, ndod, data_coo, indi_coo, indj_coo, B, indj_B, indi_B)

    ! Get B transpose in CSR format
    call coo2csr(1, nc, ndod, data_coo, indj_coo, indi_coo, BT, indj_BT, indi_BT)

end subroutine create_block_B

subroutine create_block_C(nr, ndod, dod, indi_C, indj_C, C, indi_CT, indj_CT, CT)
    !! Creates the block matrix C, a matrix of size nrxnc. C M C' = Mnn
    !! Returns C and CT in CSR format
    !! In this special matrix, the number of non zero values is equal to the number of rows (nc = nr + ndod, nr = ndof)

    implicit none
    ! Input/output data
    ! -----------------
    integer, intent(in) :: nr, ndod, dod
    dimension :: dod(ndod)

    integer, intent(out) :: indi_C, indj_C, indi_CT, indj_CT
    dimension :: indi_C(nr+1), indj_C(nr), indi_CT(nr+ndod+1), indj_CT(nr)
    double precision, intent(out) :: C, CT
    dimension :: C(nr), CT(nr)
    
    ! Local data
    ! -----------------
    integer :: i, j, nc, indi_coo, indj_coo
    dimension :: indi_coo(nr), indj_coo(nr)
    integer, allocatable, dimension(:) :: dof
    double precision :: data_coo
    dimension :: data_coo(nr, 1)

    ! Initialize 
    nc = nr + ndod
    if (any(dod.ge.nc+1)) stop 'Problem creating C'
    
    ! Get dof as complement of dod
    allocate(dof(nr))
    dof = 1; i = 1; j = 1
    do while ((j.le.nc).and.(i.le.nr))
        if (any(dod.eq.j)) then
            continue
        else
            dof(i) = j
            i = i + 1 
        end if
        j = j + 1
    end do

    ! Get COO format
    do i = 1, nr ! ndof = nr
        indi_coo(i) = i
        indj_coo(i) = dof(i)
        data_coo(i, 1) = 1.d0 
    end do

    ! Get C in CSR format
    call coo2csr(1, nr, nr, data_coo, indi_coo, indj_coo, C, indj_C, indi_C)

    ! Get C transpose in CSR format
    call coo2csr(1, nc, nr, data_coo, indj_coo, indi_coo, CT, indj_CT, indi_CT)

end subroutine create_block_C

! STEADY HEAT TRANSFER CASE:
! -------------------------------------------------------------
! With Lagrange multipliers method
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! To be modified

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
    x = array_in(:nr_total)
    y = array_in(nr_total+1:)

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
    array_out(:nr_total) = Kx + BTL
    array_out(nr_total+1:) = Bx

end subroutine mf_wq_get_Au_shlm_3d

subroutine fd_shlm_3d(nr_total, ndod, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, pardiag, phydiag, dod, array_in, array_out)
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
    integer, intent(in) :: nr_total, ndod, nr_u, nr_v, nr_w, dod
    dimension :: dod(ndod)
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag, pardiag, phydiag
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigendiag(nr_total), pardiag(nr_total), phydiag(nr_total) 

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total+ndod)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total+ndod)

    ! Local data
    ! ---------------
    integer :: i, j
    double precision :: r1, r2, s1, s2
    dimension :: r1(nr_total), r2(ndod), s1(nr_total), s2(ndod)

    ! Assign blocks 
    r1 = array_in(:nr_total)
    r2 = array_in(nr_total+1:)

    ! Solve firts block
    s1 = 0.d0
    call fd_sqr_scaling(nr_total, pardiag, phydiag, r1)
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, r1, s1)
    call fd_sqr_scaling(nr_total, pardiag, phydiag, s1)

    ! Solve second block
    s2 = 0.d0
    do i = 1, ndod
        j = dod(i)
        s2(i) = r2(i)*pardiag(j)/phydiag(j) ! to review: Identity
    end do

    ! Assign values
    array_out(:nr_total) = s1
    array_out(nr_total+1:) = s2

end subroutine fd_shlm_3d

subroutine mf_wq_solve_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, x, energy)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition), where M = K + B' P B and P a penalty diagonal matrix
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
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, xt, Kx
    dimension ::    r(nr_total+ndod), rhat(nr_total+ndod), p(nr_total+ndod), s(nr_total+ndod), ptilde(nr_total+ndod), &
                    Aptilde(nr_total+ndod), Astilde(nr_total+ndod), stilde(nr_total+ndod), xt(nr_total+ndod), Kx(nr_total)
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
        ! Where b = [f + B' P g
        !            g ]
        ! And M = K + B' P B
        ! So, A = [ M  B'   
        !           B  0 ] 
        ! To solve this system we initialize r = b - A u and u = [x, y]

        ! Initialize
        xt = 0.d0; xt(dod) = g
        x = xt(:nr_total)
        call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Kx)

        ! Create block B
        call create_block_B(nr_total, ndod, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT) 

        ! Define r
        r = 0.d0
        call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, penalty*g, BTPg)
        r(:nr_total) = f - Kx

        rhat = r; p = r
        rsold = dot_product(r, rhat)
        energy(1) = 0.d0 ! Energy : 0.5 u' A u - u' f  
    
        do iter = 1, nbIter

            call fd_shlm_3d(nr_total, ndod, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, dod, p, ptilde)
            call mf_wq_get_Au_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_shlm_3d(nr_total, ndod, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, dod, s, stilde)
            call mf_wq_get_Au_shlm_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, stilde, Astilde)

            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            xt = xt + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            ! Save result
            x = xt(:nr_total)

            call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Kx)

            energy(iter) = 0.5 * dot_product(x, Kx) - dot_product(x, f) ! Energy : 0.5 u' A u - u' f       
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

subroutine fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, pardiag, phydiag, array_in, array_out)
    !! Solves MM.s = r (MM is the preconditioner) in steady heat 3D case with lagrange multipliers, where 
    !! MM is an approximation of M = K + B' P B and P a penalty diagonal matrix. 
    !! With B a block matrix of 0 and 1. To compute MM^-1 we use fast diagonalization.
    !! Indices must be in CSR format

    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag, pardiag, phydiag
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigendiag(nr_total), pardiag(nr_total), phydiag(nr_total)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! -------------
    double precision :: array_temp(nr_total)

    ! Initialize
    array_temp = array_in

    ! Scaling
    call fd_sqr_scaling(nr_total, pardiag, phydiag, array_temp)

    ! By now, we test this approximation. It could change later
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, array_temp, array_out)

    ! Scaling
    call fd_sqr_scaling(nr_total, pardiag, phydiag, array_out)

end subroutine fd_shp_3d

subroutine mf_wq_solve_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, x, energy)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition), where M = K + B' P B and P a penalty diagonal matrix
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
    double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, Ax
    dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), ptilde(nr_total), &
                    Aptilde(nr_total), Astilde(nr_total), stilde(nr_total), Ax(nr_total)
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
        ! Where b = f + B' P g and A = K + B' P B
        ! So, to solve this system we initialize r = b - A u, and B u = g
        ! Then r = f - A u

        ! Initialize
        x = 0.d0; x(dod) = g
        call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)

        ! Create block B
        call create_block_B(nr_total, ndod, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT) 

        ! Define r
        call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, penalty*g, BTPg)
        r = f - Ax

        rhat = r; p = r
        rsold = dot_product(r, rhat)
        energy(1) = 0.d0 ! Energy : 0.5 u' A u - u' f

        do iter = 1, nbIter

            call fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, p, ptilde)
            call mf_wq_get_Au_shp_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_B, indj_B, indi_BT, indj_BT, B, BT, penalty, ndod, ptilde, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde

            call fd_shp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, s, stilde)
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

! With substitution method

subroutine mf_wq_get_Au_shs_3d(Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, indi_C, indj_C, indi_CT, indj_CT, C, CT, &
                                ndof, array_in, array_out)
    !! Computes Ann.u in steady heat 3D case with substitution method, where 
    !! But A is given as [Ann, And; Adn, Add]. So Ann u =  C A C' u, where C is a zeros and ones matrix.
    !! Indices must be in CSR format

    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, parameter :: d = 3 
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
    double precision, intent(in) :: Kcoefs
    dimension :: Kcoefs(d, d, nc_total)

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

    integer, intent(in) :: indi_C, indj_C, indi_CT, indj_CT
    dimension :: indi_C(ndof+1), indj_C(ndof), indi_CT(nr_total+1), indj_CT(ndof)
    double precision :: C, CT
    dimension :: C(ndof), CT(ndof)

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data 
    ! ------------------   
    double precision :: KBTu, BTu
    dimension :: KBTu(nr_total), BTu(nr_total)

    ! Compute B'.u                   
    call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, array_in, BTu)

    ! Compute K.u
    call mf_wq_get_ku_3d(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, BTu, KBTu)

    ! Compute B.u                   
    call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, KBTu, array_out)

end subroutine mf_wq_get_Au_shs_3d

subroutine fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, pardiag, phydiag, method, &
                    ndof, indi_C, indj_C, indi_CT, indj_CT, C, CT, array_in, array_out)
    !! Solves MM.s = r (MM is the preconditioner) in steady heat 3D (or transient) case with substitution method, where 
    !! MM is an approximation of M = K.
    !! Indices must be in CSR format

    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w, ndof
    double precision, intent(in) :: U_u, U_v, U_w, eigendiag, pardiag, phydiag
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigendiag(nr_total), pardiag(nr_total), phydiag(nr_total)

    integer, intent(in) :: indi_C, indj_C, indi_CT, indj_CT
    dimension :: indi_C(ndof+1), indj_C(ndof), indi_CT(nr_total+1), indj_CT(ndof)
    double precision :: C, CT
    dimension :: C(ndof), CT(ndof)
    character(len=10), intent(in) :: method

    double precision, intent(in) :: array_in
    dimension :: array_in(ndof)

    double precision, intent(out) :: array_out
    dimension :: array_out(ndof)

    ! Local data
    ! -------------
    double precision :: BTu, KBTu
    dimension :: BTu(nr_total), KBTu(nr_total)

    ! Compute B'.u                   
    call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, array_in, BTu)

    ! Scaling
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, BTu)
    end if

    ! By now, we test this approximation. It could change later
    call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigendiag, BTu, KBTu)

    ! Scaling
    if ((method.eq.'JMS').or.(method.eq.'TDS')) then
        call fd_sqr_scaling(nr_total, pardiag, phydiag, KBTu)
    end if

    ! Compute B.u                   
    call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, KBTu, array_out)

end subroutine fd_tshs_3d

subroutine mf_wq_solve_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, ndod, dod, f, g, nbIter, epsilon, method, nnz_cond, cond, JJ, x, energy)
    !! Precontionned bi-conjugate gradient to solve steady heat problems
    !! We want to solve M x = F, with Bx = g (Dirichlet condition). Then, we use substitution method, so
    !! Mnn xn = Fn - Mnd xd and xd = g 
    !! CSR FORMAT

    use tensor_methods
    implicit none 
    ! Input / output data
    ! ---------------------
    integer, parameter :: d = 3
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
    character(len=10), intent(in) :: method
    integer, intent(in) :: nbIter, nnz_cond
    double precision, intent(in) :: cond
    dimension :: cond(3, 3, nnz_cond)
    double precision, intent(in) :: epsilon, f, g, JJ
    dimension :: f(nr_total), g(ndod), JJ(3, 3, nc_total)
    
    double precision, intent(out) :: x, energy
    dimension :: x(nr_total), energy(nbIter)

    ! Local data
    ! ------------------
    double precision :: dummy_tol
    ! Pre / Conjugate gradient algoritm
    double precision :: Lu, Lv, Lw, lamb_u, lamb_v, lamb_w
    double precision :: rsold, rsnew, alpha, omega, beta
    double precision :: r, rhat, p, Ap, s, As, ptilde, Aptilde, stilde, Astilde, xn, Ax, fAx
    dimension ::    r(nr_total-ndod), rhat(nr_total-ndod), p(nr_total-ndod), Ap(nr_total-ndod), &
                    As(nr_total-ndod), s(nr_total-ndod), ptilde(nr_total-ndod), Aptilde(nr_total-ndod), &
                    Astilde(nr_total-ndod), stilde(nr_total-ndod), xn(nr_total-ndod), Ax(nr_total), fAx(nr_total)
    integer :: iter

    ! Fast diagonalization
    double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
    double precision, dimension(:), allocatable :: Dparametric, Dphysical, Deigen
    double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
    double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

    ! Block C
    integer :: ndof
    integer, allocatable, dimension(:) :: indi_C, indj_C, indi_CT, indj_CT
    double precision, allocatable, dimension(:) :: C, CT

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    ! Initialize
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    ! Initiate variables
    ndof = nr_total - ndod
    x = 0.d0; xn = 0.d0; energy = 0.d0; dummy_tol = epsilon
    x(dod) = g

    if (any(dod.le.0)) stop 'Indices must be greater than 0'

    ! Create block C
    allocate(indi_C(ndof+1), indj_C(ndof), indi_CT(nr_total+1), indj_CT(ndof), C(ndof), CT(ndof))
    call create_block_C(ndof, ndod, dod, indi_C, indj_C, C, indi_CT, indj_CT, CT)     

    if (method.eq.'WP') then 
        if (nbIter.gt.0) then
            ! ----------------------------
            ! Conjugate Gradient algorithm
            ! ----------------------------
            ! The system to solve is Knn xn = bn
            ! Where bn = fn - Knd xd 
            ! So, to solve this system we initialize r = bn - Knn xn, if xn = 0, r = bn
            ! Initialize
            call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)

            fAx = f - Ax 
            call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, fAx, r)
            rhat = r; p = r
            rsold = dot_product(r, rhat)
            energy(1) = 1.d0

            do iter = 1, nbIter
                call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, p, Ap)
                alpha = rsold/dot_product(Ap, rhat)
                s = r - alpha*Ap

                call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, s, As)
                omega = dot_product(As, s)/dot_product(As, As)
                xn = xn + alpha*p + omega*s
                r = s - omega*As

                call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, xn, x)
                x(dod) = g
                call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)
                
                energy(iter) = 0.5 * dot_product(x, Ax) - dot_product(x, f) ! Energy : 0.5 u' A u - u' f          
                ! if (energy(iter+1).le.epsilon) exit

                rsnew = dot_product(r, rhat)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Ap)
                rsold = rsnew
            end do

        end if

    else  
        ! Initialize
        allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
        Mcoef_u = 1.d0; Kcoef_u = 1.d0
        Mcoef_v = 1.d0; Kcoef_v = 1.d0
        Mcoef_w = 1.d0; Kcoef_w = 1.d0

        if ((method.eq.'TDS').or.(method.eq.'TDC')) then 
            ! --------------------------------------------
            ! DIAGONAL DECOMPOSITION
            ! --------------------------------------------            
            do iter = 1, 2
                call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, coefs, &
                                            Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
            end do

        else if ((method.eq.'JMS').or.(method.eq.'JMC')) then 
            ! --------------------------------------------
            ! MY METHOD
            ! --------------------------------------------
            ! Find dimensions and conductivity
            call jacobien_mean_3d(nc_u, nc_v, nc_w, nc_total, JJ, Lu, Lv, Lw)
            call conductivity_mean_3d(nc_u, nc_v, nc_w, nnz_cond, cond, lamb_u, lamb_v, lamb_w)
            
            Kcoef_u = lamb_u*Lv*Lw/Lu
            Kcoef_v = lamb_v*Lw*Lu/Lv
            Kcoef_w = lamb_w*Lu*Lv/Lw

        end if

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

        if ((method.eq.'TDS').or.(method.eq.'JMS')) then
            ! --------------------------------------------
            ! SCALING
            ! --------------------------------------------
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
        end if

        if (nbIter.gt.0) then
            ! -------------------------------------------
            ! Preconditioned Conjugate Gradient algorithm
            ! -------------------------------------------
            ! The system to solve is Knn xn = bn
            ! Where bn = fn - Knd xd 
            ! So, to solve this system we initialize r = bn - Knn xn, if xn = 0, r = bn
            ! Initialize
            call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, x, Ax)

            fAx = f - Ax 
            call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, fAx, r)
            rhat = r; p = r
            rsold = dot_product(r, rhat)
            energy(1) = 1.d0

            do iter = 1, nbIter

                call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                                ndof, indi_C, indj_C, indi_CT, indj_CT, C, CT, p, ptilde)

                call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, ptilde, Aptilde)

                alpha = rsold/dot_product(Aptilde, rhat)
                s = r - alpha*Aptilde

                call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
                                ndof, indi_C, indj_C, indi_CT, indj_CT, C, CT, s, stilde)

                call mf_wq_get_Au_shs_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                                    indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, stilde, Astilde)

                omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
                xn = xn + alpha*ptilde + omega*stilde
                r = s - omega*Astilde    
                
                call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, xn, x)
                x(dod) = g
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
    end if

end subroutine mf_wq_solve_shs_3d

! ! TRANSIENT HEAT TRANSFER CASE: 
! ! -------------------------------------------------------------
! subroutine mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                                 indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
!                                 data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                                 data_W_u, data_W_v, data_W_w, indi_C, indj_C, indi_CT, indj_CT, C, CT, &
!                                 ndof, alpha, dt, array_in, array_out)
!     !! Computes Ann.u in steady heat 3D case with substitution method, where 
!     !! But A is given as [Ann, And; Adn, Add]. So Ann u =  C A C' u, where C is a zeros and ones matrix.
!     !! Indices must be in CSR format

!     use tensor_methods
!     implicit none 
!     ! Input / output 
!     ! -------------------
!     integer, parameter :: d = 3 
!     integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndof
!     double precision, intent(in) :: Kcoefs, Ccoefs
!     dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)

!     integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
!     dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
!     integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
!     dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
!     double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
!     dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

!     integer, intent(in) :: indi_u, indi_v, indi_w
!     dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
!     integer, intent(in) :: indj_u, indj_v, indj_w
!     dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
!     double precision, intent(in) :: data_W_u, data_W_v, data_W_w
!     dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

!     integer, intent(in) :: indi_C, indj_C, indi_CT, indj_CT
!     dimension :: indi_C(ndof+1), indj_C(ndof), indi_CT(nr_total+1), indj_CT(ndof)
!     double precision :: C, CT
!     dimension :: C(ndof), CT(ndof)

!     double precision, intent(in) :: alpha, dt, array_in
!     dimension :: array_in(ndof)

!     double precision, intent(out) :: array_out
!     dimension :: array_out(ndof)

!     ! Local data 
!     ! ------------------   
!     double precision :: CBTu, KBTu, BTu
!     dimension :: CBTu(nr_total), KBTu(nr_total), BTu(nr_total)

!     ! Compute B'.u                   
!     call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, array_in, BTu)

!     ! Compute K.u
!     call mf_wq_get_ku_3d(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                         indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
!                         data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                         data_W_u, data_W_v, data_W_w, BTu, KBTu)

!     ! Compute C.u
!     call mf_wq_get_cu_3d(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                         indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
!                         indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, BTu, CBTu)

!     ! Compute alpha dt Ku + Cu
!     CBTu = CBTu + alpha*dt*KBTu

!     ! Compute B.u                   
!     call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, CBTu, array_out)

! end subroutine mf_wq_get_Au_ths_3d

! subroutine mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                             nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                             data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
!                             newmark, dt, table, ndod, dod, f, g, nbIter, epsilon, x, energy)

!     use tensor_methods
!     implicit none 
!     ! Input / output data
!     ! ---------------------
!     integer, parameter :: d = 3
!     integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ndod
!     double precision, intent(in) :: Kcoefs, Ccoefs
!     dimension :: Kcoefs(d, d, nc_total), Ccoefs(nc_total)
!     integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
!     dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
!                     indi_v(nr_v+1), indj_v(nnz_v), &
!                     indi_w(nr_w+1), indj_w(nnz_w)
!     double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
!     dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
!                     data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
!                     data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

!     double precision, intent(in) :: newmark, dt
!     integer, intent(in) :: table, dod
!     dimension :: table(d, 2), dod(ndod)
!     integer, intent(in) :: nbIter
!     double precision, intent(in) :: epsilon, f, g
!     dimension :: f(nr_total), g(ndod)
    
!     double precision, intent(out) :: x, energy
!     dimension :: x(nr_total), energy(nbIter)

!     ! Local data
!     ! ------------------
!     character (len=10) :: method = 'TDS'
!     double precision :: dummy_tol
!     ! Pre / Conjugate gradient algoritm
!     double precision :: rsold, rsnew, alpha, omega, beta
!     double precision :: r, rhat, p, s, ptilde, Aptilde, stilde, Astilde, xn, Ax, fAx
!     dimension ::    r(nr_total-ndod), rhat(nr_total-ndod), p(nr_total-ndod), &
!                     s(nr_total-ndod), ptilde(nr_total-ndod), Aptilde(nr_total-ndod), &
!                     Astilde(nr_total-ndod), stilde(nr_total-ndod), xn(nr_total-ndod), Ax(nr_total), fAx(nr_total)
!     integer :: iter

!     ! Fast diagonalization
!     double precision, dimension(:), allocatable :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
!     double precision, dimension(:), allocatable :: Kdiag_u, Kdiag_v, Kdiag_w, Mdiag_u, Mdiag_v, Mdiag_w
!     double precision, dimension(:), allocatable :: Dparametric, Dtemp, Dphysical, Deigen
!     double precision, dimension(:, :), allocatable :: U_u, U_v, U_w
!     double precision, dimension(:), allocatable :: D_u, D_v, D_w, I_u, I_v, I_w

!     ! Block C
!     integer :: ndof
!     integer, allocatable, dimension(:) :: indi_C, indj_C, indi_CT, indj_CT
!     double precision, allocatable, dimension(:) :: C, CT

!     ! Csr format
!     integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
!     dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
!                     indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
!     double precision :: data_BT_u, data_BT_v, data_BT_w
!     dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

!     ! Initialize
!     call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
!     call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
!     call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

!     ! Initiate variables
!     ndof = nr_total - ndod
!     x = 0.d0; xn = 0.d0; energy = 0.d0; dummy_tol = epsilon
!     x(dod) = g

!     if (any(dod.le.0)) stop 'Indices must be greater than 0'

!     ! Create block C
!     allocate(indi_C(ndof+1), indj_C(ndof), indi_CT(nr_total+1), indj_CT(ndof), C(ndof), CT(ndof))
!     call create_block_C(ndof, ndod, dod, indi_C, indj_C, C, indi_CT, indj_CT, CT)     

!     ! Initialize
!     allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w))            
!     Mcoef_u = 1.d0; Kcoef_u = 1.d0
!     Mcoef_v = 1.d0; Kcoef_v = 1.d0
!     Mcoef_w = 1.d0; Kcoef_w = 1.d0

!     ! --------------------------------------------
!     ! DIAGONAL DECOMPOSITION
!     ! --------------------------------------------            
!     do iter = 1, 2
!         call tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, Kcoefs, &
!                                     Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)
!     end do

!     ! --------------------------------------------
!     ! EIGEN DECOMPOSITION
!     ! -------------------------------------------- 
!     allocate(U_u(nr_u, nr_u), D_u(nr_u), U_v(nr_v, nr_v), D_v(nr_v), U_w(nr_w, nr_w), D_w(nr_w))
!     allocate(Kdiag_u(nr_u), Mdiag_u(nr_u))
!     call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
!                             data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
!                             data_W_u(:, 4), table(1, :), D_u, U_u, Kdiag_u, Mdiag_u)

!     allocate(Kdiag_v(nr_v), Mdiag_v(nr_v))
!     call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
!                             data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
!                             data_W_v(:, 4), table(2, :), D_v, U_v, Kdiag_v, Mdiag_v)    

!     allocate(Kdiag_w(nr_w), Mdiag_w(nr_w))
!     call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
!                             data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
!                             data_W_w(:, 4), table(3, :), D_w, U_w, Kdiag_w, Mdiag_w)   
!     deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

!     ! Find diagonal of eigen values
!     allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
!     allocate(Deigen(nr_total))
!     I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
!     call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, Deigen)
!     ! Update eigen diagonal
!     Deigen = 1.d0 + newmark*dt*Deigen
!     deallocate(I_u, I_v, I_w)

!     ! --------------------------------------------
!     ! SCALING
!     ! --------------------------------------------
!     ! Find diagonal of preconditioner
!     allocate(Dparametric(nr_total), Dtemp(nr_total))
!     call kron_product_3vec(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, Dtemp, 1.d0)
!     call find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
!                                 Kdiag_u, Kdiag_v, Kdiag_w, Dparametric)
!     ! Update diagonal
!     Dparametric = Dtemp + newmark*dt*Dparametric
!     deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

!     ! Find diagonal of real matrix (K in this case)
!     allocate(Dphysical(nr_total))
!     call wq_find_capacity_diagonal_3d(Ccoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                             nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                             data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dtemp)

!     call wq_find_conductivity_diagonal_3D(Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                             nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                             data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, Dphysical)
!     ! Update diagonal
!     Dphysical = Dtemp + newmark*dt*Dphysical

!     if (nbIter.gt.0) then
!         ! -------------------------------------------
!         ! Preconditioned Conjugate Gradient algorithm
!         ! -------------------------------------------
!         ! The system to solve is Knn xn = bn
!         ! Where bn = fn - Knd xd 
!         ! So, to solve this system we initialize r = bn - Knn xn, if xn = 0, r = bn
!         ! Initialize
!         call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                             indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
!                             indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
!                             1.d0, newmark*dt, x, Ax)

!         fAx = f - Ax 
!         call spMdotdV(ndof, nr_total, ndof, indi_C, indj_C, C, fAx, r)
!         rhat = r; p = r
!         rsold = dot_product(r, rhat)
!         energy(1) = 1.d0

!         do iter = 1, nbIter

!             call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
!                             ndof, indi_C, indj_C, indi_CT, indj_CT, C, CT, p, ptilde)

!             call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, & 
!                                 nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                                 indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
!                                 indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
!                                 indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, newmark, dt, ptilde, Aptilde)

!             alpha = rsold/dot_product(Aptilde, rhat)
!             s = r - alpha*Aptilde

!             call fd_tshs_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, Dparametric, Dphysical, method, &
!                             ndof, indi_C, indj_C, indi_CT, indj_CT, C, CT, s, stilde)

!             call mf_wq_get_Au_ths_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, &
!                                 nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                                 indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
!                                 indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
!                                 indi_C, indj_C, indi_CT, indj_CT, C, CT, ndof, newmark, dt, stilde, Astilde)

!             omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
!             xn = xn + alpha*ptilde + omega*stilde
!             r = s - omega*Astilde    
            
!             call spMdotdV(nr_total, ndof, ndof, indi_CT, indj_CT, CT, xn, x)
!             x(dod) = g
!             call mf_wq_get_kcu_3d(Ccoefs, Kcoefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                             indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
!                             indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
!                             1.d0, newmark*dt, x, Ax)
            
!             energy(iter) = 0.5 * dot_product(x, Ax) - dot_product(x, f) ! Energy : 0.5 u' A u - u' f          
!             ! if (energy(iter+1).le.epsilon) exit

!             rsnew = dot_product(r, rhat)
!             beta = (alpha/omega)*(rsnew/rsold)
!             p = r + beta*(p - omega*Aptilde)
!             rsold = rsnew

!         end do

!     end if

! end subroutine mf_wq_solve_ths_linear_3d

! subroutine mf_wq_ths_nonlinear_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                         nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                         data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, nbpts, table_cond, table_cap, &
!                         newmark, table, ndod, dod, invJ, detJ, sizeF, time_list, Fext, temperature, speed)
!     use tensor_methods
!     use heat_transfer
!     implicit none 
!     ! Input / output data
!     ! ---------------------
!     double precision, parameter :: tol = 1.d-8
!     ! Geometry
!     integer, parameter :: nbIterRaphson = 30, nbIterSolver = 300, d = 3
!     integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, sizeF, nbpts
!     integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
!     dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
!                     indi_v(nr_v+1), indj_v(nnz_v), &
!                     indi_w(nr_w+1), indj_w(nnz_w)
!     double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
!     dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
!                     data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
!                     data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
!     ! Physics
!     double precision, intent(in) :: table_cond, table_cap, newmark
!     dimension :: table_cond(nbpts, 2), table_cap(nbpts, 2)
!     integer, intent(in) :: ndod, table, dod
!     dimension :: table(d, 2), dod(ndod)
!     double precision, intent(in) :: invJ, detJ, Fext, time_list
!     dimension :: invJ(d, d, nc_total), detJ(nc_total), Fext(nr_total, sizeF+1), time_list(sizeF+1) 
    
!     double precision, intent(out) :: temperature, speed
!     dimension :: temperature(nr_total, sizeF+1), speed(nr_total, sizeF+1)

!     ! Local data
!     ! -----------   
!     double precision :: g(ndod), energy(nbIterSolver)
!     double precision :: delta_time, TTn0, VVn0, ddVV, TTn1, TTn1_t0, VVn1, F, TT_interp
!     dimension ::    TTn0(nr_total), VVn0(nr_total), ddVV(nr_total), TTn1(nr_total), TTn1_t0(nr_total), &
!                     VVn1(nr_total), F(nr_total), TT_interp(nc_total)
!     double precision :: KK, CC, Kcoefs, Ccoefs, Fint, CdT, KT, ddFF
!     dimension ::    KK(nc_total), CC(nc_total), Kcoefs(d, d, nc_total), Ccoefs(nc_total), &
!                     Fint(nr_total), CdT(nr_total), KT(nr_total), ddFF(nr_total)
!     double precision :: relerror, prod1, prod2
!     integer :: i, j

!     ! Initialize
!     g = 0.d0 !!!!!! Pour l'instant on considère Dirichlet = 0

!     ! --------------------------------------------
!     ! SOLVE
!     ! -------------------------------------------- 
!     ! Initialize
!     temperature = 0.d0; speed = 0.d0
    
!     do i = 2, sizeF+1
!         ! Get delta time
!         delta_time = time_list(i) - time_list(i-1)

!         ! Initialize 
!         TTn0 = temperature(:, i-1)
!         VVn0 = speed(:, i-1)   
!         ddVV = 0.d0

!         ! Prediction of new step
!         TTn1 = TTn0 + delta_time*(1-newmark)*VVn0
!         TTn1_t0 = TTn1
!         VVn1 = 0.d0

!         ! Get force of new step
!         F = Fext(:, i)
!         prod2 = dot_product(F, F)

!         ! Newton Raphson
!         do j = 1, nbIterRaphson
!             print*, 'Step: ', i-1, ' Iteration: ', j-1

!             ! Compute temperature (at each quadrature point) 
!             call interpolate_fieldphy_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                                 indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                                 data_B_u, data_B_v, data_B_w, 1, TTn1, TT_interp)

!             ! Interpolate capacity and conductivity at each quadrature point 
!             call compute_heat_properties(nbpts, table_cond, table_cap, nc_total, TT_interp, KK, CC)

!             ! Compute coefficients to compute tangent matrix
!             call compute_heat_coefficients(nc_total, KK, CC, invJ, detJ, Kcoefs, Ccoefs)
            
!             ! Compute Fint = C dT + K T 
!             call mf_wq_get_cu_3d_csr(Ccoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                                 nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                                 data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, VVn1, CdT)

!             call mf_wq_get_ku_3d_csr(Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                                 nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                                 data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, TTn1, KT)

!             Fint = CdT + KT

!             ! Compute residue
!             ddFF = F - Fint
!             ddFF(dod) = 0.d0 
!             prod1 = dot_product(ddFF, ddFF)
!             relerror = sqrt(prod1/prod2)
!             print*, "Raphson with error: ", relerror
!             if (isnan(relerror)) stop
            
!             ! Verify
!             if (relerror.le.1e-6) then 
!                 exit
!             else

!                 ! Solve by iterations 
!                 call mf_wq_solve_ths_linear_3d(Ccoefs, Kcoefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                                 nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                                 data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
!                                 newmark, delta_time, table, ndod, dod, ddFF, g, nbIterSolver, tol, ddVV, energy)

!                 ! Update values
!                 VVn1 = VVn1 + ddVV
!                 TTn1 = TTn1_t0 + newmark*delta_time*VVn1

!             end if                
!         end do
        
!         ! Set values
!         speed(:, i) = VVn1
!         temperature(:, i) = TTn1
                
!     end do

! end subroutine mf_wq_ths_nonlinear_3d