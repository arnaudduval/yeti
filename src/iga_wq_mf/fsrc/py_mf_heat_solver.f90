! ====================================================
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! It is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

subroutine mf_capacity_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, isLumped, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v    
    integer, intent(in) :: indi_u, indi_v
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
    integer, intent(in) ::  indj_u, indj_v
    dimension :: indj_u(nnz_u), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    logical, intent(in) :: isLumped
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_2basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, data_W_u, data_W_v)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    mat%isLumped = isLumped
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, prop)
    call mf_u_v(mat, basisdata, nr_total, array_in, array_out)

end subroutine mf_capacity_2d

subroutine mf_capacity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isLumped, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes C.u where C is capacity matrix in 3D
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w    
    integer, intent(in) :: indi_u, indi_v, indi_w
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    logical, intent(in) :: isLumped
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w,&
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    mat%isLumped = isLumped
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, prop)
    call mf_u_v(mat, basisdata, nr_total, array_in, array_out)
end subroutine mf_capacity_3d

subroutine mf_conductivity_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat 
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_2basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v,&
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                        data_W_u, data_W_v)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_gradu_gradv(mat, basisdata, nr_total, array_in, array_out)
    
end subroutine mf_conductivity_2d

subroutine mf_conductivity_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w,&
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_gradu_gradv(mat, basisdata, nr_total, array_in, array_out)
    
end subroutine mf_conductivity_3d

subroutine mf_thmchcoupled_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                        data_B_u, data_B_v, data_W_u, data_W_v, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_2basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                        data_W_u, data_W_v)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_thmchcoupledprop(mat, nc_total, prop)
    call mf_gradu_tv(mat, basisdata, nr_total, array_in, array_out)
    
end subroutine mf_thmchcoupled_2d

subroutine mf_thmchcoupled_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                        invJ, detJ, prop, array_in, array_out)
    !! Computes K.u where K is conductivity matrix in 3D 
    !! This function is adapted to python
    !! IN CSR FORMAT

    use matrixfreeheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(basis_data) :: basisdata
    call init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w,&
                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, &
                        data_W_u, data_W_v, data_W_w)
    call getcsrc2dense(basisdata)
    call getcsr2csc(basisdata)
    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_thmchcoupledprop(mat, nc_total, prop)
    call mf_gradu_tv(mat, basisdata, nr_total, array_in, array_out)
    
end subroutine mf_thmchcoupled_3d

subroutine solver_linearsteady_heat_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, table, &
                            invJ, detJ, prop, Fext, iterations, threshold, linprecond, x, residual)
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
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    logical, intent(in) :: table
    dimension :: table(dimen, 2)

    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: linprecond
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv    
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
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
    call setup_conductivityprop(mat, nc_total, prop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then
        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if
        
        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen, mat%Kmean(1:dimen))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        if (solv%applyfd) call space_eigendecomposition(redsyst)
        call initialize_solver(solv, globsyst, redsyst)
        call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
    else 
        stop 'Unknown method'  
    end if

end subroutine solver_linearsteady_heat_2d

subroutine solver_linearsteady_heat_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            table, invJ, detJ, prop, Fext, iterations, threshold, linprecond, x, residual)
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
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
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

    double precision, intent(in) :: invJ, detJ, prop
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), prop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: linprecond
    integer, intent(in) :: iterations
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv    
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
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
    call setup_conductivityprop(mat, nc_total, prop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v, nc_w/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if
        
        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen, mat%Kmean(1:dimen))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        if (solv%applyfd) call space_eigendecomposition(redsyst)
        call initialize_solver(solv, globsyst, redsyst)
        call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)

    else 
        stop 'Unknown method' 
    end if

end subroutine solver_linearsteady_heat_3d

subroutine solver_lineartransient_heat_2d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, isLumped, table, &
                                invJ, detJ, Cprop, Kprop, tsfactor, Fext, iterations, threshold, linprecond, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    logical, intent(in) :: isLumped
    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ, Cprop, Kprop, tsfactor
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Cprop(nc_total), Kprop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: linprecond
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv    
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
    call init_2basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_u, indj_u, indi_v, indj_v, data_B_u, data_B_v, &
                        data_W_u, data_W_v)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsr2csc(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    call getcsr2csc(globsyst)

    call setup_geometry(mat, dimen, nc_total, invJ, detJ)
    call setup_capacityprop(mat, nc_total, Cprop)
    call setup_conductivityprop(mat, nc_total, Kprop)
    mat%isLumped = isLumped
    mat%scalars = (/1.d0, tsfactor/)
    solv%matrixfreetype = 3
    nc_list = (/nc_u, nc_v/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen, mat%Kmean)
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if
        
        if (solv%applyfd) call space_eigendecomposition(redsyst)
        if (solv%applyfd) redsyst%diageigval_sp = mat%Cmean + tsfactor*redsyst%diageigval_sp
        call initialize_solver(solv, globsyst, redsyst)
        call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_lineartransient_heat_2d

subroutine solver_lineartransient_heat_3d(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, isLumped, table, &
                                invJ, detJ, Cprop, Kprop, tsfactor, Fext, iterations, threshold, linprecond, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreeheat
    use heatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    logical, intent(in) :: isLumped
    logical, intent(in) :: table
    dimension :: table(dimen, 2) 

    double precision, intent(in) :: invJ, detJ, Cprop, Kprop, tsfactor
    dimension :: invJ(dimen, dimen, nc_total), detJ(nc_total), Cprop(nc_total), Kprop(dimen, dimen, nc_total)
    character(len=10), intent(in) :: linprecond
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(thermomat) :: mat
    type(cgsolver) :: solv    
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen)
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
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
    call setup_capacityprop(mat, nc_total, Cprop)
    call setup_conductivityprop(mat, nc_total, Kprop)
    mat%isLumped = isLumped
    mat%scalars = (/1.d0, tsfactor/)
    solv%matrixfreetype = 3
    nc_list = (/nc_u, nc_v, nc_w/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen, mat%Kmean)
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if
        
        if (solv%applyfd) call space_eigendecomposition(redsyst)
        if (solv%applyfd) redsyst%diageigval_sp = mat%Cmean + tsfactor*redsyst%diageigval_sp
        call initialize_solver(solv, globsyst, redsyst)
        call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_lineartransient_heat_3d
