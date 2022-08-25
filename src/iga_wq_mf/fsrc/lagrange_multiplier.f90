! ==========================
! module: Lagrange multipliers to solve Dirichlet problem
! author: Joaquin Cornejo
! 
! In this module, one could find functions to solve Ax = f using Lagrange multipliers
! That is, instead of having Mx = f, with Bx = g (Dirichlet condition), we solve
! [ M  B'    [x  = [f
!   B  0 ]    y]    g]
! with y, the lagrange multiplier vector. This equation is the kind A u = b
! By the moment, these function has not been proven 
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

! STEADY HEAT TRANSFER CASE: M = K + B' P B
subroutine mf_wq_get_Au_steady_lagrange_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_W_u, data_W_v, data_W_w, penalty, ndod, dod, array_input, array_output)
    !! Computes A.u in 3D case, where 
    !! A = [ M  B'   
    !!       B  0 ] with M = K + B' P B and P a penalty diagonal matrix = p Identity
    !! Array input is composed of 2 blocks = [T lambda], then the output will be [M.T + B'.lambda, B.T]
    !! This, if we replace, [K.T + B'.(p B.T + lambda), B.T] 
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

    integer, intent(in) :: dod
    dimension :: dod(ndod)

    double precision, intent(in) :: array_input
    dimension :: array_input(nr_total+ndod)

    double precision, intent(out) :: array_output
    dimension :: array_output(nr_total+ndod)

    ! Local data 
    ! ------------------
    integer :: indi_B, indj_B, indi_BT, indj_BT
    dimension :: indi_B(ndod+1), indj_B(ndod), indi_BT(nr_total+1), indj_BT(ndod)
    double precision :: B, BT
    dimension ::  B(ndod), BT(ndod)
    double precision :: KKTH, BBTH, L, BTL
    dimension :: KKTH(nr_total), BBTH(ndod), L(ndod), BTL(nr_total)

    ! Computes the first part of the output vector
    ! --------------------------------------------
    ! Compute K.Th
    call mf_wq_get_ku_3d(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_input(1:nr_total), KKTH)

    ! Compute B.Th
    call create_block_B(ndod, nr_total, dod, indi_B, indj_B, B, indi_BT, indj_BT, BT)                    
    call spMdotdV(ndod, nr_total, ndod, indi_B, indj_B, B, array_input(1:nr_total), BBTH)

    ! Compute p B.Th + lambda
    L = penalty*BBTH + array_input(nr_total+1:nr_total+ndod)
    
    ! Compute K.Th + B'.L
    call spMdotdV(nr_total, ndod, ndod, indi_BT, indj_BT, BT, L, BTL)

    ! Set the first part
    array_output(1:nr_total) = KKTH + BTL

    ! Computes the second part of the output vector
    ! --------------------------------------------
    array_output(nr_total+1:nr_total+ndod) = BBTH

end subroutine mf_wq_get_Au_steady_lagrange_3d

! subroutine fd_lagrange_steady_heat()
! end subroutine fd_lagrange_steady_heat

! subroutine mf_wq_solve_steady_lagrange_3d()
! end subroutine mf_wq_solve_steady_lagrange_3d