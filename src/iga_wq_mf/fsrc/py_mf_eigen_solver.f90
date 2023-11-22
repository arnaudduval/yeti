subroutine solver_helmholtz_lobpcg_2d(nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, table, &
                            invJ, detJ, ishigher, nbIterPCG, threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver2
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    logical, intent(in) :: table
    dimension :: table(dimen, 2)

    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold
    
    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(nr_u*nr_v)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: i, nc_list(dimen)
    double precision :: Kprop(dimen, dimen, nc_total), Cprop(nc_total), mean(dimen) = 1.d0

    ! Csr format
    integer :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v)
    double precision :: data_BT_u, data_BT_v
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)

    Cprop = 1.d0; Kprop = 0.d0
    do i = 1, dimen
        Kprop(i, i, :) = 1.d0
    end do

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, Kprop)
    call setup_capacityprop(mat, nc_total, Cprop)
    nc_list = (/nc_u, nc_v/)

    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                            data_W_u, data_W_v, table, mean)

    call LOBPCGSTAB(solv, mat, nr_u*nr_v, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, ishigher, nbIterPCG, threshold, eigenvec, eigenval)

end subroutine solver_helmholtz_lobpcg_2d

subroutine solver_helmholtz_lobpcg_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, invJ, detJ, ishigher, nbIterPCG, threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver3
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: nbIterPCG
    double precision, intent(in) :: threshold

    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    integer :: i, nc_list(dimen)
    double precision :: Kprop(dimen, dimen, nc_total), Cprop(nc_total), mean(dimen) = 1.d0

    ! Csr format
    integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
    dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                    indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
    
    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    Cprop = 1.d0; Kprop = 0.d0
    do i = 1, dimen
        Kprop(i, i, :) = 1.d0
    end do

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, Kprop)
    call setup_capacityprop(mat, nc_total, Cprop)
    nc_list = (/nc_u, nc_v, nc_w/)

    
    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                            data_W_u, data_W_v, data_W_w, table, mean)

    call LOBPCGSTAB(solv, mat, nr_u*nr_v*nr_w, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, data_BT_u, data_BT_v, data_BT_w, &
                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, &
                ishigher, nbIterPCG, threshold, eigenvec, eigenval)

end subroutine solver_helmholtz_lobpcg_3d

subroutine fastdiagonalization_2d(nr_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use heatsolver2
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)
    
    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver) :: solv
    double precision :: meanval(dimen) = 1.d0
    logical :: table(dimen, 2) = .true.

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                            data_W_u, data_W_v, table, meanval)

    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine fastdiagonalization_2d

subroutine fastdiagonalization_3d(nr_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use heatsolver3
    use datastructure
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)
    
    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(cgsolver) :: solv
    double precision :: meanval(dimen) = 1.d0
    logical :: table(dimen, 2) = .true.

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                            data_W_u, data_W_v, data_W_w, table, meanval)

    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine fastdiagonalization_3d