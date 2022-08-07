! ====================================================
! module :: Fast diagonalization 
! author :: Joaquin Cornejo
! ====================================================

subroutine eigen_decomposition_py(nr, nc, nnz, indi, indj, data_B0, data_W0, data_B1, data_W1, &
                                Method, Mcoef, Kcoef, t_robin, eigenvalues, eigenvectors)
    !! Eigen decomposition generalized KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    use tensor_methods
    implicit none 
    ! Input / output 
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B0, data_W0, data_B1, data_W1
    dimension :: data_B0(nnz), data_W0(nnz), data_B1(nnz), data_W1(nnz)
    character(len=10), intent(in) :: Method
    double precision, intent(in) :: Mcoef, Kcoef
    dimension :: Mcoef(nc), Kcoef(nc)
    integer, intent(in) :: t_robin
    dimension :: t_robin(2)
            
    double precision, intent(out) :: eigenvalues, eigenvectors
    dimension :: eigenvalues(nr), eigenvectors(nr, nr)

    ! Local data
    ! -----------------
    double precision :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    call eigen_decomposition(nr, nc, Mcoef, Kcoef, nnz, indi, indj, &
                            data_B0, data_W0, data_B1, data_W1, Method, t_robin, &
                            eigenvalues, eigenvectors, Kdiag, Mdiag)

end subroutine eigen_decomposition_py

subroutine fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, &
                                diagonal, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to steady heat equations
    !! by G. Sanaglli and M. Tani
    
    use tensor_methods
    use omp_lib
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    diagonal(nr_total), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! -------------
    integer :: i, nb_tasks
    double precision :: array_temp
    dimension :: array_temp(nr_total)

    ! First part 
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                    transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

    ! Second part 
    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
    do i = 1, nr_total
        array_temp(i) = array_temp(i)/diagonal(i)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! Third part
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                    U_u, U_v, U_w, array_temp, array_out)
    
end subroutine fd_steady_heat_3d

subroutine fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! by G. Sanaglli and M. Tani
    
    use tensor_methods
    use omp_lib
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, array_in
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! -------------
    double precision :: array_temp
    dimension :: array_temp(nr_total)

    ! First part 
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                    transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

    ! Second part
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, array_out)

end subroutine fd_interpolation_3d

subroutine fd_transient_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, diagonal, array_in, dt, theta, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to transient heat equations
    !! by G. Sanaglli and M. Tani
    
    use tensor_methods
    use omp_lib
    implicit none
    ! Input / output  data 
    !---------------------
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
    dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                    diagonal(nr_total), array_in(nr_total)
    double precision, intent(in) :: dt, theta

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! -------------
    integer :: i, nb_tasks
    double precision :: array_temp, diagonal_new
    dimension :: array_temp(nr_total), diagonal_new(nr_total)

    ! First part 
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

    ! Second part 
    diagonal_new = 1.d0/dt
    if (theta.gt.0.d0) then 
        diagonal_new = diagonal_new + theta*diagonal
    end if
    !$OMP PARALLEL 
    nb_tasks = omp_get_num_threads()
    !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
    do i = 1, nr_total
        array_temp(i) = array_temp(i)/diagonal_new(i)
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    ! Third part
    call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, array_out)

end subroutine fd_transient_heat_3d

subroutine fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, diagonal, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    use omp_lib
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: d = 3
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
    dimension ::    U_u(nr_u, nr_u, d), U_v(nr_v, nr_v, d), U_w(nr_w, nr_w, d), &
                    diagonal(d, nr_total), array_in(d, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(d, nr_total)

    ! Local data
    ! -------------
    integer :: i
    double precision :: array_temp
    dimension :: array_temp(nr_total)

    do i = 1, d 
        call fd_steady_heat_3d(nr_total, nr_u, nr_v, nr_w, U_u(:, :, i), U_v(:, :, i), U_w(:, :, i), &
                                diagonal(i, :), array_in(i, :), array_temp)
        array_out(i, :) = array_temp
    end do
    
end subroutine fd_elasticity_3d