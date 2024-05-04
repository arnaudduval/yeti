subroutine solver_eig_heat_2d(nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, table, &
                            invJ, detJ, Cprop, Kprop, ishigher, iterations, &
                            threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver
    use structured_data
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

    double precision, intent(in) :: invJ, detJ, Kprop, Cprop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Kprop(dimen, dimen, nc_total), Cprop(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold
    
    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(nr_u*nr_v)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nr_total, nc_list(dimen)

    call init_2basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v,&
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v,  &
                        data_W_u, data_W_v)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsr2csc(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    call getcsr2csc(globsyst)

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, Kprop)
    call setup_capacityprop(mat, nc_total, Cprop)

    nc_list = (/nc_u, nc_v/)
    nr_total = nr_u*nr_v
    solv%withdiag = .true.
    call compute_variablesmean(mat, nc_list)
    call setup_meancoefs(redsyst, dimen, mat%Kmean(1:dimen))

    call space_eigendecomposition(redsyst)
    call initialize_solver(solv, globsyst, redsyst)
    call RQMIN(solv, mat, nr_total, ishigher, iterations, threshold, eigenvec, eigenval)

end subroutine solver_eig_heat_2d

subroutine solver_eig_heat_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, invJ, detJ, Cprop, Kprop, ishigher, iterations, threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver
    use structured_data
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

    double precision, intent(in) :: invJ, detJ, Kprop, Cprop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Kprop(dimen, dimen, nc_total), Cprop(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold

    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nr_total, nc_list(dimen)
    
    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w,&
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsr2csc(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    call getcsr2csc(globsyst)

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, Kprop)
    call setup_capacityprop(mat, nc_total, Cprop)

    nc_list = (/nc_u, nc_v, nc_w/)
    nr_total = nr_u*nr_v*nr_w
    solv%withdiag = .true.
    call compute_variablesmean(mat, nc_list)
    call setup_meancoefs(redsyst, dimen, mat%Kmean(1:dimen))

    call space_eigendecomposition(redsyst)
    call initialize_solver(solv, globsyst, redsyst)
    call RQMIN(solv, mat, nr_total, ishigher, iterations, threshold, eigenvec, eigenval)

end subroutine solver_eig_heat_3d

subroutine solver_eig_elasticity_2d(nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            table, invJ, detJ, nbmechArgs, mechArgs, Mprop, ishigher, iterations, &
                            threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeplasticity
    use plasticitysolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, nbmechArgs
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    logical, intent(in) :: table
    dimension :: table(dimen, 2, dimen) 

    double precision, intent(in) :: invJ, detJ, mechArgs, Mprop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), mechArgs(nbmechArgs, nc_total), Mprop(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold

    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(dimen, nr_u*nr_v)

    ! Local data
    ! ----------
    type(mecamat) :: mat
    type(cgsolver) :: solv
    type(basis_data) :: globsyst
    type(reduced_system), target :: redsyst(dimen)
    integer :: i, nr_total, nc_list(dimen)

    call init_2basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v,&
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v,  &
                        data_W_u, data_W_v)
    do i = 1, dimen
        call copybasisdata(globsyst, redsyst(i)%basisdata)
        call update_reducedsystem(redsyst(i), dimen, table(:, :, i))
        call getcsrc2dense(redsyst(i)%basisdata)
        call getcsr2csc(redsyst(i)%basisdata)
    end do
    call getcsrc2dense(globsyst)
    call getcsr2csc(globsyst)

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_mechanicalArguments(mat, nbmechArgs, mechArgs)
    call setup_massprop(mat, nc_total, Mprop)

    nc_list = (/nc_u, nc_v/)
    nr_total = nr_u*nr_v
    solv%withdiag = .true.
    call compute_variablesmean(mat, nc_list)
    do i = 1, dimen
        call setup_meancoefs(redsyst(i), dimen, mat%Smean(i, 1:dimen))
    end do

    do i = 1, dimen
        call space_eigendecomposition(redsyst(i))
    end do
    call initialize_solver(solv, globsyst, redsyst)
    call RQMIN(solv, mat, nr_total, ishigher, iterations, threshold, eigenvec, eigenval)

end subroutine solver_eig_elasticity_2d

subroutine solver_eig_elasticity_3d(nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, invJ, detJ, nbmechArgs, mechArgs, Mprop, ishigher, iterations, &
                            threshold, eigenval, eigenvec)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use matrixfreeplasticity
    use plasticitysolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, nbmechArgs
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    logical, intent(in) :: table
    dimension :: table(dimen, 2, dimen) 

    double precision, intent(in) :: invJ, detJ, mechArgs, Mprop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), mechArgs(nbmechArgs, nc_total), Mprop(nc_total)
    logical, intent(in) :: ishigher
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold

    double precision, intent(out) :: eigenvec, eigenval
    dimension :: eigenvec(dimen, nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    type(mecamat) :: mat
    type(cgsolver) :: solv
    type(basis_data) :: globsyst
    type(reduced_system), target :: redsyst(dimen)
    integer :: i, nr_total, nc_list(dimen)

    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, &
                        data_B_w, data_W_u, data_W_v, data_W_w)
    do i = 1, dimen
        call copybasisdata(globsyst, redsyst(i)%basisdata)
        call update_reducedsystem(redsyst(i), dimen, table(:, :, i))
        call getcsrc2dense(redsyst(i)%basisdata)
        call getcsr2csc(redsyst(i)%basisdata)
        call space_eigendecomposition(redsyst(i))
    end do
    call getcsrc2dense(globsyst)
    call getcsr2csc(globsyst)

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_jacobienjacobien(mat)
    call setup_mechanicalArguments(mat, nbmechArgs, mechArgs)
    call setup_massprop(mat, nc_total, Mprop)

    nc_list = (/nc_u, nc_v, nc_w/)
    nr_total = nr_u*nr_v
    solv%withdiag = .true.
    call compute_variablesmean(mat, nc_list)
    do i = 1, dimen
        call setup_meancoefs(redsyst(i), dimen, mat%Smean(i, 1:dimen))
    end do

    do i = 1, dimen
        call space_eigendecomposition(redsyst(i))
    end do
    call initialize_solver(solv, globsyst, redsyst)
    call RQMIN(solv, mat, nr_total, ishigher, iterations, threshold, eigenvec, eigenval)

end subroutine solver_eig_elasticity_3d

subroutine fastdiagonalization_2d(nr_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use heatsolver
    use structured_data
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
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    logical :: table(dimen, 2) = .true.

    call init_2basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                        data_W_u, data_W_v)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    call space_eigendecomposition(redsyst)

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
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

    use heatsolver
    use structured_data
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
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    logical :: table(dimen, 2) = .true.

    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, &
                        data_B_w, data_W_u, data_W_v, data_W_w)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    call space_eigendecomposition(redsyst)

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine fastdiagonalization_3d

subroutine elfastdiagonalization_2d(nr_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use plasticitysolver
    use structured_data
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
    dimension :: array_in(dimen, nr_total)
    
    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst(dimen)
    logical :: table(dimen, 2, dimen) = .true.
    integer :: i

    call init_2basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                        data_W_u, data_W_v)
    do i = 1, dimen
        call copybasisdata(globsyst, redsyst(i)%basisdata)
        call update_reducedsystem(redsyst(i), dimen, table(:, :, i))
        call getcsrc2dense(redsyst(i)%basisdata)
    end do
    call getcsrc2dense(globsyst)
    do i = 1, dimen
        call space_eigendecomposition(redsyst(i))
    end do

    if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine elfastdiagonalization_2d

subroutine elfastdiagonalization_3d(nr_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use plasticitysolver
    use structured_data
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
    dimension :: array_in(dimen, nr_total)
    
    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst(dimen)
    logical :: table(dimen, 2, dimen) = .true.
    integer :: i

    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w)
    do i = 1, dimen
        call copybasisdata(globsyst, redsyst(i)%basisdata)
        call update_reducedsystem(redsyst(i), dimen, table(:, :, i))
        call getcsrc2dense(redsyst(i)%basisdata)
    end do
    call getcsrc2dense(globsyst)
    do i = 1, dimen
        call space_eigendecomposition(redsyst(i))
    end do

    if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine elfastdiagonalization_3d

subroutine sptfastdiagonalization_2d(nr_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
        nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
        indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, &
        data_W_t, array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    logical :: table(dimen_sp+1, 2) = .true.
    integer :: dimen

    dimen = dimen_sp + 1; table(dimen, 2) = .false.
    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_B_u, data_B_v, data_B_t, &
                        data_W_u, data_W_v, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.
    
    call space_eigendecomposition(redsyst)
    call time_schurdecomposition(redsyst)

    if (nr_total.ne.nr_u*nr_v*nr_t) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine sptfastdiagonalization_2d

subroutine sptfastdiagonalization_3d(nr_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
        indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, &
        data_W_w, data_W_t, array_in, array_out)
    !! (Preconditioned) Conjugate gradient algorithm to solver linear heat problems 
    !! It solves Ann xn = bn, where Ann is Knn (steady heat problem) and bn = Fn - And xd
    !! bn is compute beforehand (In python).
    !! IN CSR FORMAT

    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)

    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    logical :: table(dimen_sp+1, 2) = .true.
    integer :: dimen

    dimen = dimen_sp + 1; table(dimen, 2) = .false.
    call init_4basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                        data_W_u, data_W_v, data_W_w, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.

    call space_eigendecomposition(redsyst)
    call time_schurdecomposition(redsyst)

    if (nr_total.ne.nr_u*nr_v*nr_w*nr_t) stop 'Size problem'
    call initialize_solver(solv, globsyst, redsyst)
    call applyfastdiag(solv, nr_total, array_in, array_out)

end subroutine sptfastdiagonalization_3d