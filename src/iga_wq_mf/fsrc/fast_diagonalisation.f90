! ====================================================
! module :: Fast diagonalization 
! author :: Joaquin Cornejo
! ====================================================

subroutine eigen_decomposition_py(nr, nc, nnz, indi, indj, data_B0, data_W0, data_B1, data_W1, &
                                Mcoef, Kcoef, robin_condition, eigenvalues, eigenvectors)
    !! Eigen decomposition generalized KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B0, data_W0, data_B1, data_W1
    dimension :: data_B0(nnz), data_W0(nnz), data_B1(nnz), data_W1(nnz)
    double precision, intent(in) :: Mcoef, Kcoef
    dimension :: Mcoef(nc), Kcoef(nc)
    integer, intent(in) :: robin_condition
    dimension :: robin_condition(2)
            
    double precision, intent(out) :: eigenvalues, eigenvectors
    dimension :: eigenvalues(nr), eigenvectors(nr, nr)

    ! Local data
    ! ----------
    double precision :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    call eigen_decomposition(nr, nc, Mcoef, Kcoef, nnz, indi, indj, &
                            data_B0, data_W0, data_B1, data_W1, robin_condition, &
                            eigenvalues, eigenvectors, Kdiag, Mdiag)

end subroutine eigen_decomposition_py

subroutine fd_steady_heat_3d_py(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to steady heat problems
    !! by G. Sanaglli and M. Tani
    
    use solverheat
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    eigen_diag(nr_total), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call setup_preconditionerdiag(solv, nr_total, eigen_diag)
    call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)

end subroutine fd_steady_heat_3d_py

subroutine fd_interpolation_3d_py(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applieg in control points interpolation problems
    !! by G. Sanaglli and M. Tani
    
    use solverheat
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, array_in
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv

    allocate(solv)
    call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)

end subroutine fd_interpolation_3d_py